using Test
using MagestyRebuild
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
pair_cutoff = 1.5
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
pair_cutoff = 1.5
lmax = [2]
"""

_writetoml(s) = (p = tempname() * ".toml"; write(p, s); p)

@testset "TOML input files" begin
    @testset "read_input parses structure / interaction / symmetry" begin
        inp = read_input(_writetoml(_INPUT_FULL))
        @test num_atoms(inp.crystal) == 2
        @test inp.crystal.species_labels == ["Fe"]
        @test inp.crystal.species == [1, 1]
        @test inp.crystal.lattice.vectors == SMatrix{3,3,Float64}(3.0 * I)
        @test inp.crystal.frac_positions[:, 1] ≈ [0.2, 0.0, 0.0]
        @test inp.crystal.frac_positions[:, 2] ≈ [0.8, 0.0, 0.0]
        @test inp.crystal.lattice.pbc == SVector{3,Bool}(true, true, true)
        @test inp.interaction.nbody == 2
        @test inp.interaction.pair_cutoff == 1.5
        @test inp.interaction.lmax == [2]
        @test inp.interaction.isotropy == false
        @test inp.backend isa NoSymmetry
        @test inp.tol == 1.0e-5
    end

    @testset "SCEBasis(path) == building from the same Crystal/Interaction" begin
        path = _writetoml(_INPUT_FULL)
        b_file = SCEBasis(path)
        inp = read_input(path)
        b_manual = SCEBasis(inp.crystal, inp.interaction)
        @test b_file.salcs.keys == b_manual.salcs.keys
        @test nsalc(b_file) == nsalc(b_manual)
    end

    @testset "defaults for omitted keys" begin
        inp = read_input(_writetoml(_INPUT_MINIMAL))
        @test inp.backend isa NoSymmetry          # no [symmetry] → NoSymmetry
        @test inp.tol == 1.0e-5                    # default tol
        @test inp.interaction.isotropy == false    # default isotropy
        @test inp.crystal.lattice.pbc == SVector{3,Bool}(true, true, true)  # default pbc
    end

    @testset "keyword arguments override the file's [symmetry]" begin
        path = _writetoml(_INPUT_FULL)
        @test SCEBasis(path; tol = 1e-3).spacegroup.tol == 1e-3
        @test SCEBasis(path).spacegroup.tol == 1e-5
    end

    @testset "backend name maps to the requested backend type" begin
        s = replace(_INPUT_FULL, "backend = \"none\"" => "backend = \"spglib\"")
        inp = read_input(_writetoml(s))           # mapping only; building would need `using Spglib`
        @test inp.backend isa SpglibBackend
    end

    @testset "error paths" begin
        only_interaction = "[interaction]\nnbody = 1\npair_cutoff = 1.5\nlmax = [2]\n"
        @test_throws ArgumentError read_input(_writetoml(only_interaction))   # no [structure]

        only_structure = """
        [structure]
        lattice = [[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]]
        positions = [[0.0,0.0,0.0]]
        species = [1]
        species_labels = ["Fe"]
        """
        @test_throws ArgumentError read_input(_writetoml(only_structure))     # no [interaction]

        # missing required key inside a section
        no_lattice = replace(_INPUT_FULL,
            "lattice = [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]]\n" => "")
        @test_throws ArgumentError read_input(_writetoml(no_lattice))

        # lmax length ≠ number of species
        two_species = replace(_INPUT_FULL, "species_labels = [\"Fe\"]" => "species_labels = [\"Fe\", \"Pt\"]")
        @test_throws ArgumentError read_input(_writetoml(two_species))        # lmax=[2] for 2 species

        # unrecognized backend
        bad_backend = replace(_INPUT_FULL, "backend = \"none\"" => "backend = \"xml\"")
        @test_throws ArgumentError read_input(_writetoml(bad_backend))

        # malformed lattice (only 2 vectors)
        bad_lattice = replace(_INPUT_FULL,
            "lattice = [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]]" =>
                "lattice = [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0]]")
        @test_throws ArgumentError read_input(_writetoml(bad_lattice))
    end
end
