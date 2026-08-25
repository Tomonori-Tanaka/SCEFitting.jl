using Test
using SCEFitting
using LinearAlgebra
using StaticArrays

const _INPUT_FULL = """
[structure]
lattice = [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]]
positions = [[0.2, 0.0, 0.0], [0.8, 0.0, 0.0]]
species = [1, 1]
species_labels = ["Fe"]
pbc = [true, true, true]

[interaction]
nbody = 2
cutoff = 1.5
lmax = [2]
isotropy = false

[symmetry]
backend = "none"
tol = 1.0e-5
"""

# minimal: no [symmetry], no isotropy, no pbc → defaults apply
const _INPUT_MINIMAL = """
[structure]
lattice = [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]]
positions = [[0.0, 0.0, 0.0]]
species = [1]
species_labels = ["Fe"]

[interaction]
nbody = 1
cutoff = 1.5
lmax = [2]
"""

_writetoml(s) = (p = tempname() * ".toml"; write(p, s); p)

@testset "TOML input files" begin
    @testset "read_setup parses structure / interaction / symmetry" begin
        inp = read_setup(_writetoml(_INPUT_FULL))
        @test n_atoms(inp.crystal) == 2
        @test inp.crystal.species_labels == ["Fe"]
        @test inp.crystal.species == [1, 1]
        @test inp.crystal.lattice.vectors == SMatrix{3,3,Float64}(3.0 * I)
        @test inp.crystal.frac_positions[:, 1] ≈ [0.2, 0.0, 0.0]
        @test inp.crystal.frac_positions[:, 2] ≈ [0.8, 0.0, 0.0]
        @test inp.crystal.lattice.pbc == SVector{3,Bool}(true, true, true)
        @test inp.spec.nbody == 2
        @test inp.spec.cutoff[1][1, 1] == 1.5
        @test inp.spec.lmax == [2]
        @test inp.spec.isotropy == false
        @test inp.backend isa NoSymmetry
        @test inp.tol == 1.0e-5
    end

    @testset "SCEBasis(path) == building from the same Crystal/BasisSpec" begin
        path = _writetoml(_INPUT_FULL)
        b_file = SCEBasis(path)
        inp = read_setup(path)
        b_manual = SCEBasis(inp.crystal, inp.spec)
        @test b_file.salc_basis.keys == b_manual.salc_basis.keys
        @test n_salcs(b_file) == n_salcs(b_manual)
    end

    @testset "defaults for omitted keys" begin
        inp = read_setup(_writetoml(_INPUT_MINIMAL))
        @test inp.backend isa NoSymmetry          # no [symmetry] → NoSymmetry
        @test inp.tol == 1.0e-5                    # default tol
        @test inp.spec.isotropy == false    # default isotropy
        @test inp.crystal.lattice.pbc == SVector{3,Bool}(true, true, true)  # default pbc
        @test inp.tie_tol == 1.0e-8                # default same-distance band
    end

    # `tie_tol` changes the emitted basis (a widened band merges near-tie shells), so
    # a setup that needed one must round-trip through its own file — same rule as
    # `images`. The keyword overrides the file, and the file value reaches the
    # constructor's validation (the cap refuses).
    @testset "[interaction].tie_tol is carried and overridable" begin
        s = replace(_INPUT_FULL, "nbody = 2" => "nbody = 2\ntie_tol = 1e-5")
        inp = read_setup(_writetoml(s))
        @test inp.tie_tol == 1.0e-5
        @test SCEBasis(_writetoml(s)) isa SCEBasis                  # builds with it
        @test SCEBasis(_writetoml(s); tie_tol = 1e-7) isa SCEBasis  # override accepted
        bad = replace(_INPUT_FULL, "nbody = 2" => "nbody = 2\ntie_tol = 0.5")
        @test_throws ArgumentError SCEBasis(_writetoml(bad))        # cap enforced
    end

    @testset "keyword arguments override the file's [symmetry]" begin
        path = _writetoml(_INPUT_FULL)
        @test SCEBasis(path; tol = 1e-3).spacegroup.tol == 1e-3
        @test SCEBasis(path).spacegroup.tol == 1e-5
    end

    @testset "backend name maps to the requested backend type" begin
        s = replace(_INPUT_FULL, "backend = \"none\"" => "backend = \"spglib\"")
        inp = read_setup(_writetoml(s))           # mapping only; building would need `using Spglib`
        @test inp.backend isa SpglibBackend
    end

    @testset "error paths" begin
        only_interaction = "[interaction]\nnbody = 1\ncutoff = 1.5\nlmax = [2]\n"
        @test_throws ArgumentError read_setup(_writetoml(only_interaction))   # no [structure]

        only_structure = """
        [structure]
        lattice = [[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]]
        positions = [[0.0,0.0,0.0]]
        species = [1]
        species_labels = ["Fe"]
        """
        @test_throws ArgumentError read_setup(_writetoml(only_structure))     # no [interaction]

        # missing required key inside a section
        no_lattice = replace(_INPUT_FULL,
            "lattice = [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]]\n" => "")
        @test_throws ArgumentError read_setup(_writetoml(no_lattice))

        # lmax length ≠ number of species
        two_species = replace(_INPUT_FULL, "species_labels = [\"Fe\"]" => "species_labels = [\"Fe\", \"Pt\"]")
        @test_throws ArgumentError read_setup(_writetoml(two_species))        # lmax=[2] for 2 species

        # unrecognized backend
        bad_backend = replace(_INPUT_FULL, "backend = \"none\"" => "backend = \"xml\"")
        @test_throws ArgumentError read_setup(_writetoml(bad_backend))

        # malformed lattice (only 2 vectors)
        bad_lattice = replace(_INPUT_FULL,
            "lattice = [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]]" =>
                "lattice = [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0]]")
        @test_throws ArgumentError read_setup(_writetoml(bad_lattice))
    end
end

# ── [moment] ───────────────────────────────────────────────────────────────────────
# Oracles: a hand-written `MomentSpec(; ...)` keyword call (never the parser) and the
# documented defaults; the basis gate compares against the existing Julia path.

# 2-species cell: Fe chain + one Rh; [interaction] is the minimum a setup file needs.
const _INPUT_FERH = """
[structure]
lattice = [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 6.0]]
positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.5], [0.5, 0.5, 0.25]]
species = [1, 1, 2]
species_labels = ["Fe", "Rh"]

[interaction]
nbody = 2
cutoff = 3.1
lmax = [2, 0]
"""

const _MOMENT_FULL = """
[moment]
nbody       = 3
lmax_mark   = 2
lmax_env    = [2, 0]
sampled     = ["Fe"]
marked      = ["Fe", "Rh"]
cutoff_star = 3.1
lsum        = 4
isotropy    = true

[moment.cutoff_pair]
"Fe-Fe" = 4.1
"*-*"   = 3.0
"""

# Same spec in the label-table / boolean spellings.
const _MOMENT_FULL_ALT = """
[moment]
nbody       = 3
lmax_mark   = 2
sampled     = [true, false]
marked      = ["*"]
cutoff_star = 3.1
lsum        = 4
isotropy    = true

[moment.lmax_env]
"*" = 0
Fe  = 2

[moment.cutoff_pair]
"Fe-Fe" = 4.1
"Fe-Rh" = 3.0
"Rh-Rh" = 3.0
"""

const _MOMENT_MINIMAL = """
[moment]
lmax_env    = [1, 0]
sampled     = ["Fe"]
cutoff_pair = 3.1
"""

_ferh_toml(moment::String) = _writetoml(_INPUT_FERH * "\n" * moment)

@testset "[moment] section" begin
    @testset "fields equal the hand-written MomentSpec" begin
        expected = MomentSpec(; lmax_env = [2, 0], sampled = [true, false],
                              lmax_mark = 2, marked = [true, true], nbody = 3,
                              cutoff_pair = [4.1 3.0; 3.0 3.0],
                              cutoff_star = 3.1, lsum = 4, isotropy = true)
        for text in (_MOMENT_FULL, _MOMENT_FULL_ALT)
            got = read_setup(_ferh_toml(text)).moment
            @test got isa MomentSpec
            for f in fieldnames(MomentSpec)
                @test getfield(got, f) == getfield(expected, f)
            end
        end
    end

    @testset "documented defaults for omitted keys" begin
        m = read_setup(_ferh_toml(_MOMENT_MINIMAL)).moment
        @test m.nbody == 3
        @test m.lmax_mark == 2
        @test m.marked == [true, true]
        @test m.cutoff_pair == fill(3.1, 2, 2)
        @test m.cutoff_star == [m.cutoff_pair]   # one entry per star order (nbody = 3)
        @test m.lsum == fill(typemax(Int), 3)     # uncapped, one entry per body order
        @test m.isotropy == true                  # NOT the [interaction] default (false)
        @test read_setup(_ferh_toml(_MOMENT_MINIMAL)).spec.isotropy == false
    end

    @testset "files without [moment] are unchanged" begin
        @test read_setup(_writetoml(_INPUT_FULL)).moment === nothing
        @test read_setup(_writetoml(_INPUT_MINIMAL)).moment === nothing
        @test_throws ArgumentError MomentBasis(_writetoml(_INPUT_FULL))
        try
            MomentBasis(_writetoml(_INPUT_FULL))
        catch err
            @test occursin("[moment]", err.msg)
        end
    end

    @testset "MomentBasis(path) == building from the same Crystal/MomentSpec" begin
        path = _ferh_toml(_MOMENT_FULL)
        inp = read_setup(path)
        mb_file = MomentBasis(path)
        # the manual path takes the hand-written spec, not the parser's
        expected = MomentSpec(; lmax_env = [2, 0], sampled = [true, false],
                              lmax_mark = 2, marked = [true, true], nbody = 3,
                              cutoff_pair = [4.1 3.0; 3.0 3.0],
                              cutoff_star = 3.1, lsum = 4, isotropy = true)
        mb_manual = MomentBasis(inp.crystal, expected)
        @test mb_file.salc_basis.keys == mb_manual.salc_basis.keys
        @test mb_file.marked_atoms == mb_manual.marked_atoms
        @test n_salcs(mb_file) == n_salcs(mb_manual)
        @test n_salcs(mb_file) > 0
        # keywords override the file; [interaction].tie_tol is shared with the moment basis
        @test mb_file.tie_tol == SCEFitting._SAME_DIST_RTOL
        @test MomentBasis(path; tie_tol = 1e-6).tie_tol == 1e-6
        shared = _writetoml(replace(_INPUT_FERH, "nbody = 2" =>
                                    "nbody = 2\ntie_tol = 1e-5") * "\n" * _MOMENT_FULL)
        @test MomentBasis(shared).tie_tol == 1e-5
        @test MomentBasis(path; tol = 1e-3).spacegroup.tol == 1e-3
    end

    @testset "refusals" begin
        bad(text) = @test_throws ArgumentError read_setup(_ferh_toml(text))
        msg(text) = try
            read_setup(_ferh_toml(text)); ""
        catch err
            err isa ArgumentError ? err.msg : rethrow()
        end
        # required keys
        bad(replace(_MOMENT_MINIMAL, "sampled     = [\"Fe\"]\n" => ""))
        bad(replace(_MOMENT_MINIMAL, "cutoff_pair = 3.1\n" => ""))
        bad(replace(_MOMENT_MINIMAL, "lmax_env    = [1, 0]\n" => ""))
        # unknown key, upstream spelling
        @test occursin("unknown key", msg(_MOMENT_MINIMAL * "lmax_enviroment = 2\n"))
        @test occursin("upstream spelling", msg(_MOMENT_MINIMAL * "soc = false\n"))
        # species lists: unknown label, mixed kinds, integers, wrong length, duplicates
        @test occursin("unknown species label",
                       msg(replace(_MOMENT_MINIMAL, "[\"Fe\"]" => "[\"Co\"]")))
        bad(replace(_MOMENT_MINIMAL, "[\"Fe\"]" => "[\"Fe\", true]"))
        bad(replace(_MOMENT_MINIMAL, "[\"Fe\"]" => "[1, 0]"))
        bad(replace(_MOMENT_MINIMAL, "[\"Fe\"]" => "[]"))
        @test occursin("entries for", msg(replace(_MOMENT_MINIMAL, "[\"Fe\"]" => "[true]")))
        @test occursin("duplicate entry",
                       msg(replace(_MOMENT_MINIMAL, "[\"Fe\"]" => "[\"Fe\", \"Fe\"]")))
        # MomentSpec's own consistency rule surfaces unchanged (environment spins on an
        # unsampled species)
        @test occursin("not sampled", msg(replace(_MOMENT_MINIMAL, "[1, 0]" => "[1, 1]")))
        # bare scalar lmax_env, body-order tables, duplicate unordered pair key
        bad(replace(_MOMENT_MINIMAL, "[1, 0]" => "1"))
        @test occursin("not a body order",
                       msg(_MOMENT_MINIMAL * "lsum = { \"Fe\" = 4 }\n"))
        @test occursin("must be an integer",
                       msg(_MOMENT_MINIMAL * "lsum = { 2 = 4.5 }\n"))
        @test occursin("body-order", msg(replace(_MOMENT_MINIMAL, "cutoff_pair = 3.1" =>
                                                 "cutoff_pair = { 2 = 3.1 }")))
        @test occursin("duplicate", msg(replace(_MOMENT_MINIMAL, "cutoff_pair = 3.1" =>
                    "cutoff_pair = { \"Fe-Rh\" = 3.0, \"Rh-Fe\" = 3.0, \"*-*\" = 3.0 }")))
        # `cutoff_star` is the one pointed cutoff that IS per body order. The table must
        # cover exactly 3:nbody — a partial table would leave an order at the default.
        let base = _MOMENT_MINIMAL * "nbody = 4\n"
            got = read_setup(_ferh_toml(base *
                "[moment.cutoff_star]\n3 = 3.1\n[moment.cutoff_star.4]\n" *
                "\"Fe-Fe\" = 2.0\n\"*-*\" = 1.0\n")).moment
            @test got.cutoff_star == [fill(3.1, 2, 2), [2.0 1.0; 1.0 1.0]]
            # a scalar still broadcasts to every star order
            @test read_setup(_ferh_toml(base * "cutoff_star = 2.5\n")).moment.cutoff_star ==
                  [fill(2.5, 2, 2), fill(2.5, 2, 2)]
            @test occursin("must be exactly",
                           msg(base * "[moment.cutoff_star]\n3 = 3.1\n"))      # 4 missing
            @test occursin("must be exactly",
                           msg(base * "[moment.cutoff_star]\n2 = 3.1\n4 = 1.0\n"))
            @test occursin("mixes",
                           msg(base * "[moment.cutoff_star]\n3 = 3.1\n\"*-*\" = 1.0\n"))
            # `lsum` is body-keyed too, but from order 1 and PARTIAL: unnamed orders
            # stay uncapped (the [interaction].lsum spelling, not cutoff_star's).
            let m = read_setup(_ferh_toml(base * "[moment.lsum]\n3 = 6\n4 = 4\n")).moment
                @test m.lsum == [typemax(Int), typemax(Int), 6, 4]
            end
            @test read_setup(_ferh_toml(base * "lsum = 4\n")).moment.lsum == fill(4, 4)
            @test occursin("outside 1:4", msg(base * "[moment.lsum]\n5 = 4\n"))
        end
        # wrong value kinds — booleans are Integers in Julia, so they are refused by name
        bad(_MOMENT_MINIMAL * "nbody = \"three\"\n")
        bad(_MOMENT_MINIMAL * "isotropy = \"yes\"\n")
        bad(_MOMENT_MINIMAL * "isotropy = 1\n")
        bad(_MOMENT_MINIMAL * "nbody = true\n")
        bad(_MOMENT_MINIMAL * "lsum = 4.0\n")
        bad(replace(_MOMENT_MINIMAL, "[1, 0]" => "[true, false]"))
        bad(replace(_MOMENT_MINIMAL, "cutoff_pair = 3.1" => "cutoff_pair = true"))
        bad(_MOMENT_MINIMAL * "nbody = 5\n")        # MomentSpec range rule
        # a setup with all_images builds (minimum-image) but warns
        ai_setup = replace(_INPUT_FERH, "nbody = 2" => "nbody = 2\nimages = \"all_images\"")
        ai = _writetoml(ai_setup * "\n" * _MOMENT_MINIMAL)
        @test (@test_logs (:warn, r"all_images") MomentBasis(ai)) isa MomentBasis
    end
end
