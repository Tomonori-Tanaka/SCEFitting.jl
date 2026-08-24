# Real-data parity tier (S6): SCEFitting's moment channel vs upstream SLCE.jl on
# FeGe B20 (2×2×2, 64 atoms) and FeRh B2 (4×4×4, 128 atoms). See README.md for
# the three kinds of statement made here (acceptance numbers at 1 %, census
# integers as change detectors, absolute column parity) and the data locations.
#
# Both packages are imported QUALIFIED: they export the same names.

# SCOPE NOTE on body order: the acceptance cases pin `nbody = 3` (energy cases 2 or
# 3), because their reference numbers were taken there. The pointed body-order DOOR is
# covered separately by the `nbody = 4` case at the end of the FeGe block. Without it
# the two engines can diverge at N = 4 with this gate green.


using Test
using LinearAlgebra
using Random
using Statistics
using StaticArrays
using Printf
import SCEFitting
import SLCE

const FEGE_DIR = get(ENV, "SCE_PARITY_FEGE_DIR",
    expanduser("~/Packages/_brain_storming/adiabatic-moment-sce/step2_assets/fege"))
const FEGE_POSCAR = get(ENV, "SCE_PARITY_FEGE_POSCAR",
    expanduser("~/jijs/magesty/fege/2x2x2/lc_exp/ge_include/202601/mfa-tau01/300k5/" *
               "nelect_0/input/POSCAR"))
const FERH_DIR = get(ENV, "SCE_PARITY_FERH_DIR",
    expanduser("~/jijs/magesty/ferh/4x4x4/afm/tau_0to0.6/350k3"))

# ── provenance ─────────────────────────────────────────────────────────────────────

function _rev(pkgdir)
    rev = try
        strip(read(`git -C $pkgdir rev-parse --short HEAD`, String))
    catch
        "unknown"
    end
    dirty = try
        isempty(read(`git -C $pkgdir status --porcelain`, String)) ? "" : " +dirty"
    catch
        ""
    end
    return rev * dirty
end
println("SCEFitting ", _rev(dirname(dirname(pathof(SCEFitting)))),
        "  vs  SLCE ", _rev(dirname(dirname(pathof(SLCE)))),
        "  (Julia $VERSION, $(Threads.nthreads()) threads)")

# ── shared helpers ─────────────────────────────────────────────────────────────────

wrap01(x) = (y = mod(x, 1.0); y >= 1.0 - 1e-12 ? 0.0 : y)

function read_poscar(path)
    L = readlines(path)
    sc = parse(Float64, strip(L[2]))
    A = sc .* hcat([parse.(Float64, split(L[k]))[1:3] for k = 3:5]...)
    off = 0
    cnt = try parse.(Int, split(L[7])) catch; off = 1; parse.(Int, split(L[8])) end
    startswith(lowercase(strip(L[8+off])), "d") || error("not direct")
    nat = sum(cnt)
    frac = hcat([parse.(Float64, split(L[8+off+k])[1:3]) for k = 1:nat]...)
    return A, frac, cnt
end

# fixed-group backends, one per package
struct FixedA <: SCEFitting.AbstractSymmetryBackend
    sg::SCEFitting.SpaceGroup
end
SCEFitting.analyze_symmetry(b::FixedA, ::SCEFitting.Crystal; tol::Real = 1e-5) = b.sg
struct FixedB <: SLCE.AbstractSymmetryBackend
    sg::SLCE.SpaceGroup
end
SLCE.analyze_symmetry(b::FixedB, ::SLCE.Crystal; tol::Real = 1e-5) = b.sg

# the same crystal + group in both packages from one (A, frac, species, ops) set
function both_crystals(A, frac, species, labels, rots, trs, symbol, number; tol)
    xa = SCEFitting.Crystal(SCEFitting.Lattice(A), frac, species, labels)
    xb = SLCE.Crystal(SLCE.Lattice(A), frac, species, labels)
    sga = SCEFitting._assemble_spacegroup(xa, rots, trs, symbol, number; tol)
    sgb = SLCE._assemble_spacegroup(xb, rots, trs, symbol, number; tol)
    return xa, sga, xb, sgb
end

# a SALCKey as a package-neutral tuple (decors by their integer fields)
keytuple(k) = (k.body, k.orbit_id, [(d.spin_l, d.disp_k, d.disp_l) for d in k.decors],
               k.L_S, k.Lf, k.block)

# absolute column parity: match by key, then ‖a − b‖ ≤ rtol‖b‖ per column and
# elementwise atol relative to the column's max — no free scale anywhere
function column_parity(Xa, keysa, Xb, keysb; rtol = 1e-10)
    ta = keytuple.(keysa)
    tb = keytuple.(keysb)
    @test length(ta) == length(tb)
    @test Set(ta) == Set(tb)
    pos = Dict(t => j for (j, t) in enumerate(tb))
    worst = 0.0
    for (ja, t) in enumerate(ta)
        haskey(pos, t) || continue
        a = @view Xa[:, ja]
        b = @view Xb[:, pos[t]]
        nb = norm(b)
        d = norm(a .- b)
        worst = max(worst, nb == 0 ? d : d / nb)
        @test d <= rtol * nb
        @test maximum(abs, a .- b) <= 1e-12 * max(maximum(abs, b), 1e-300)
    end
    return worst
end

_skip(msg) = (@warn msg; @test_skip false)

# ══════════════════════════════════════════════════════════════════════════════════
@testset "FeGe B20 2×2×2 — acceptance numbers + column parity" begin
    files = joinpath.(FEGE_DIR, ["tau01.extxyz", "tau05.extxyz", "tau03.extxyz"])
    if !(all(isfile, files) && isfile(FEGE_POSCAR))
        _skip("FeGe data absent (SCE_PARITY_FEGE_DIR / SCE_PARITY_FEGE_POSCAR)")
    else
        # P2₁3 point operations (handwritten) × 2×2×2 translations
        WS = [[1 0 0; 0 1 0; 0 0 1], [-1 0 0; 0 -1 0; 0 0 1], [-1 0 0; 0 1 0; 0 0 -1],
              [1 0 0; 0 -1 0; 0 0 -1],
              [0 0 1; 1 0 0; 0 1 0], [0 0 1; -1 0 0; 0 -1 0], [0 0 -1; -1 0 0; 0 1 0],
              [0 0 -1; 1 0 0; 0 -1 0],
              [0 1 0; 0 0 1; 1 0 0], [0 -1 0; 0 0 1; -1 0 0], [0 1 0; 0 0 -1; -1 0 0],
              [0 -1 0; 0 0 -1; 1 0 0]]
        TS = [[0, 0, 0], [0.5, 0, 0.5], [0, 0.5, 0.5], [0.5, 0.5, 0],
              [0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
              [0, 0, 0], [0, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5]]
        A, frac, cnt = read_poscar(FEGE_POSCAR)
        species = vcat(fill(1, cnt[1]), fill(2, cnt[2]))
        rots = SMatrix{3,3,Float64}[]; trs = SVector{3,Float64}[]
        for (W, t) in zip(WS, TS), i = 0:1, j = 0:1, k = 0:1
            push!(rots, SMatrix{3,3,Float64}(Float64.(W)))
            push!(trs, SVector{3,Float64}(wrap01.((Float64.(t) .+ [i, j, k]) ./ 2)))
        end
        xa, sga, xb, sgb = both_crystals(A, wrap01.(frac), species, ["Fe", "Ge"],
                                         rots, trs, "P2_13(2x2x2)", 198; tol = 1e-4)
        nat = length(species)
        @test nat == 64 && SCEFitting.n_ops(sga) == 96 && SLCE.n_ops(sgb) == 96

        # readers: SCEFitting (reference-checked) and SLCE on the same files
        d1 = SCEFitting.read_extxyz(files[1]; reference = xa)
        d5 = SCEFitting.read_extxyz(files[2]; reference = xa)
        d3 = SCEFitting.read_extxyz(files[3]; reference = xa)
        train = vcat(d1, d5)
        @test length(train) == 400 && length(d3) == 120
        u1 = SLCE.read_extxyz(files[1])
        @test length(u1) == length(d1)
        # reader parity (bitwise on the fields the moment channel reads)
        for c in (1, length(d1))
            @test u1[c].directions == d1[c].directions
            @test u1[c].moments_bare == d1[c].moments_bare
            @test u1[c].magmoms == d1[c].magmoms
            @test u1[c].constraint_mode == d1[c].constraint_mode == 4
        end

        EPS = 2.2e-4
        # reference values: SLCE.jl on this protocol (design record §3.3, the
        # upstream fit_fege.jl log) — an implementation-agreement statement
        REF = Dict("lsum2" => (p = 39, sigma = 0.0398, cv = 0.0403, ho = 0.0308),
                   "full" => (p = 181, sigma = 0.0334, cv = 0.0358, ho = 0.0295))
        specs = [
            ("lsum2",
             SCEFitting.MomentSpec(; lmax_env = [2, 2], sampled = [true, true],
                                   lmax_mark = 2, nbody = 3, cutoff_pair = 4.6,
                                   cutoff_star = 3.0, lsum = 2, marked = [true, false]),
             SLCE.MomentSpec(; lmax_env = [2, 2], sampled = [true, true],
                             lmax_mark = 2, nbody = 3, cutoff_pair = 4.6,
                             cutoff_star = 3.0, lsum = 2, marked = [true, false])),
            ("full",
             SCEFitting.MomentSpec(; lmax_env = [2, 2], sampled = [true, true],
                                   lmax_mark = 2, nbody = 3, cutoff_pair = 4.6,
                                   cutoff_star = 3.0, marked = [true, false]),
             SLCE.MomentSpec(; lmax_env = [2, 2], sampled = [true, true],
                             lmax_mark = 2, nbody = 3, cutoff_pair = 4.6,
                             cutoff_star = 3.0, marked = [true, false])),
        ]
        for (lbl, spa, spb) in specs
            @testset "$lbl" begin
                ref = REF[lbl]
                mba = SCEFitting.MomentBasis(xa, spa; backend = FixedA(sga))
                mbb = SLCE.MomentBasis(xb, spb; backend = FixedB(sgb))
                @test SCEFitting.n_salcs(mba) == SLCE.n_salcs(mbb) == ref.p
                @test mba.marked_atoms == mbb.marked_atoms == collect(1:32)

                ds = SCEFitting.MomentDataset(mba, train; gate_eps = EPS)
                # census (change detector: what the archived data contain)
                @test (count(ds.keep), length(ds.keep)) == (12757, 12800)
                f = SCEFitting.fit(SCEFitting.MomentFit, ds)
                sigma = std(SCEFitting.residuals(f))
                # config-level 5-fold CV, the upstream protocol verbatim
                Random.seed!(1)
                perm = randperm(length(train))
                cvres = Float64[]
                for k = 1:5
                    tc = Set(perm[k:5:end])
                    te = [i for i in eachindex(ds.y) if (ds.row_config[i] in tc) && ds.keep[i]]
                    tr = [i for i in eachindex(ds.y) if !(ds.row_config[i] in tc) && ds.keep[i]]
                    bk = ds.X[tr, :] \ ds.y[tr]
                    append!(cvres, ds.y[te] - ds.X[te, :] * bk)
                end
                cv = std(cvres)
                ds3 = SCEFitting.MomentDataset(mba, d3; gate_eps = EPS)
                @test (count(ds3.keep), length(ds3.keep)) == (3838, 3840)
                rho = ds3.y[ds3.keep] - ds3.X[ds3.keep, :] * SCEFitting.coef(f)
                ho = std(rho)
                @printf("  %-5s p=%3d  sigma %.4f (ref %.4f, %+.2f %%)  CV %.4f (ref %.4f, %+.2f %%)  held-out %.4f (ref %.4f, %+.2f %%)\n",
                        lbl, ref.p, sigma, ref.sigma, 100(sigma / ref.sigma - 1),
                        cv, ref.cv, 100(cv / ref.cv - 1), ho, ref.ho, 100(ho / ref.ho - 1))
                # acceptance at relative 1 % (the reference is rounded to 4 digits,
                # so the tolerance must also cover ±0.5e-4 of rounding)
                tol(x) = 0.01 * x + 0.5e-4
                @test abs(sigma - ref.sigma) <= tol(ref.sigma)
                @test abs(cv - ref.cv) <= tol(ref.cv)
                @test abs(ho - ref.ho) <= tol(ref.ho)

                # column parity on 20 training + 10 held-out configurations
                sel = vcat(1:10, 391:400)
                cfgs = [train[c].directions for c in sel]
                Xa = SCEFitting._design_moment(mba, cfgs, [copy(e) for e in cfgs])
                Xb = SLCE._design_moment(mbb, cfgs, [copy(e) for e in cfgs])
                w = column_parity(Xa, mba.salc_basis.keys, Xb, mbb.salc_basis.keys)
                @printf("  %-5s column parity: worst relative column deviation %.2e\n", lbl, w)
                # dataset parity through SLCE's own dataset door on the same data
                # (targets, gates and masks bitwise; the arithmetic is identical)
                ub = SLCE.read_extxyz(files[3])
                dsb3 = SLCE.MomentDataset(mbb, ub; gate_eps = EPS)
                @test dsb3.y == ds3.y
                @test dsb3.gate == ds3.gate
                @test dsb3.keep == ds3.keep
                @test dsb3.row_atom == ds3.row_atom
                # held-out predictions from each package's own fit of the same rows
                dsb = SLCE.MomentDataset(mbb, vcat(SLCE.read_extxyz(files[1]),
                                                  SLCE.read_extxyz(files[2]));
                                         gate_eps = EPS)
                fb = SLCE.fit(SLCE.MomentFit, dsb)
                pa = ds3.X[ds3.keep, :] * SCEFitting.coef(f)
                pb = dsb3.X[dsb3.keep, :] * SLCE.coef(fb)
                @test norm(pa - pb) <= 1e-10 * norm(pb)
            end
        end

        # ---- the pointed body-order door (both engines general in N) -------------
        # No acceptance number here: the reference protocol never ran at N = 4. What
        # this states is the thing the rest of the harness cannot — that the two
        # packages build the SAME 4-body sector. `cutoff_star = 2.6` keeps it to the
        # nearest-neighbour shell (at 3.0 the same spec reaches 727 columns, 572 of
        # them 4-body, which is a benchmark rather than a gate).
        @testset "nbody = 4" begin
            sp4a = SCEFitting.MomentSpec(; lmax_env = [2, 2], sampled = [true, true],
                                         lmax_mark = 2, nbody = 4, cutoff_pair = 4.6,
                                         cutoff_star = 2.6, lsum = 4,
                                         marked = [true, false])
            sp4b = SLCE.MomentSpec(; lmax_env = [2, 2], sampled = [true, true],
                                   lmax_mark = 2, nbody = 4, cutoff_pair = 4.6,
                                   cutoff_star = 2.6, lsum = 4, marked = [true, false])
            mb4a = SCEFitting.MomentBasis(xa, sp4a; backend = FixedA(sga))
            mb4b = SLCE.MomentBasis(xb, sp4b; backend = FixedB(sgb))
            # the 4-body sector is REALLY there — a spec that silently returned the
            # 3-body basis would satisfy every parity assertion below
            n4a = count(k -> k.body == 4, mb4a.salc_basis.keys)
            @test n4a == count(k -> k.body == 4, mb4b.salc_basis.keys) > 0
            # The cross-package equalities are the oracle; the literals `43` and
            # `[1, 24, 10, 8]` are a REGRESSION PIN — a change detector, not evidence
            # of correctness (what makes the sector correct is the absolute
            # normalization oracle in test/unit/test_momentbasis.jl). Captured
            # 2026-08-25 on this fixture (FeGe B20 2x2x2, the spec above). Recapture
            # only when the truncation deliberately changes; a move with the spec
            # unchanged is the bug the pin exists to catch.
            @test SCEFitting.n_salcs(mb4a) == SLCE.n_salcs(mb4b) == 43
            @test [count(k -> k.body == b, mb4a.salc_basis.keys) for b = 1:4] ==
                  [count(k -> k.body == b, mb4b.salc_basis.keys) for b = 1:4] ==
                  [1, 24, 10, 8]
            sel = vcat(1:6, 396:400)
            cfgs = [train[c].directions for c in sel]
            X4a = SCEFitting._design_moment(mb4a, cfgs, [copy(e) for e in cfgs])
            X4b = SLCE._design_moment(mb4b, cfgs, [copy(e) for e in cfgs])
            w4 = column_parity(X4a, mb4a.salc_basis.keys, X4b, mb4b.salc_basis.keys)
            @printf("  nbody4 p=%3d (%d four-body)  column parity: worst %.2e\n",
                    SCEFitting.n_salcs(mb4a), n4a, w4)
        end
    end
end

# ══════════════════════════════════════════════════════════════════════════════════
@testset "FeRh B2 4×4×4 — mode-1 census + column parity" begin
    pos = joinpath(FERH_DIR, "input", "rh_random", "POSCAR")
    mwp = joinpath(FERH_DIR, "output", "rh_random", "lambda100", "EMBSET")
    mip = joinpath(FERH_DIR, "output", "rh_random", "lambda100", "EMBSET_mint")
    inc(c) = joinpath(FERH_DIR, "input", "rh_random", @sprintf("sample-%03d.INCAR", c))
    if !(isfile(pos) && isfile(mwp) && isfile(mip) && isfile(inc(1)))
        _skip("FeRh data absent (SCE_PARITY_FERH_DIR)")
    else
        NAT = 128
        A, frac, cnt = read_poscar(pos)
        @test sum(cnt) == NAT && length(cnt) == 2
        species = vcat(fill(1, cnt[1]), fill(2, cnt[2]))
        # O_h (48 signed permutations) × 4×4×4 translations
        perms3 = [(1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)]
        rots = SMatrix{3,3,Float64}[]; trs = SVector{3,Float64}[]
        for p in perms3, s1 in (1.0, -1.0), s2 in (1.0, -1.0), s3 in (1.0, -1.0)
            W = zeros(3, 3); W[1, p[1]] = s1; W[2, p[2]] = s2; W[3, p[3]] = s3
            for i = 0:3, j = 0:3, k = 0:3
                push!(rots, SMatrix{3,3,Float64}(W))
                push!(trs, SVector{3,Float64}(i / 4, j / 4, k / 4))
            end
        end
        xa, sga, xb, sgb = both_crystals(A, wrap01.(frac), species, ["Fe", "Rh"],
                                         rots, trs, "Pm-3m(4x4x4)", 221; tol = 1e-4)
        @test SCEFitting.n_ops(sga) == 3072

        # per-config constraint axes from the sample INCARs: an INDEPENDENT
        # M_CONSTR parser (join `\` continuations, take the tag value, normalize
        # each triple, zero stays zero; collinear run ⇒ SAXIS rotation = I)
        function mconstr_axes(path, nat)
            joined = String[]; buf = ""
            for raw in split(read(path, String), '\n')
                line = rstrip(raw)
                if endswith(line, '\\')
                    buf *= chop(line) * " "
                else
                    push!(joined, buf * line); buf = ""
                end
            end
            isempty(buf) || push!(joined, buf)
            val = nothing
            for line in joined
                bare = strip(first(split(line, '#')))
                m = match(r"^\s*M_CONSTR\s*=\s*(.*)$", bare)
                m === nothing || (val = m.captures[1])
            end
            val === nothing && error("no M_CONSTR in $path")
            v = parse.(Float64, split(val))
            length(v) == 3nat || error("M_CONSTR has $(length(v)) numbers in $path")
            axes = zeros(3, nat)
            for a = 1:nat
                m = v[3a-2:3a]
                n = norm(m)
                n <= 1e-12 && continue
                axes[:, a] .= m ./ n
            end
            return axes
        end
        axes = [mconstr_axes(inc(c), NAT) for c = 1:121]
        # the archive's finite-λ physics: the default p99 < 5° angle gate refuses
        # (upstream measured p99 6.53°, max 22.16°), so the disclosed 10° bound
        @test_throws ArgumentError SCEFitting.read_embset_pair(
            mwp, mip; n_atoms = NAT, constraint_mode = 1, constraint_axes = axes)
        data = SCEFitting.read_embset_pair(mwp, mip; n_atoms = NAT, constraint_mode = 1,
                                           constraint_axes = axes,
                                           axis_angle_p99_max = 10.0)
        @test length(data) == 121

        spa = SCEFitting.MomentSpec(; lmax_env = [2, 2], sampled = [true, true],
                                    lmax_mark = 2, nbody = 2, cutoff_pair = 3.2,
                                    marked = [true, true])
        spb = SLCE.MomentSpec(; lmax_env = [2, 2], sampled = [true, true],
                              lmax_mark = 2, nbody = 2, cutoff_pair = 3.2,
                              marked = [true, true])
        mba = SCEFitting.MomentBasis(xa, spa; backend = FixedA(sga))
        mbb = SLCE.MomentBasis(xb, spb; backend = FixedB(sgb))
        @test SCEFitting.n_salcs(mba) == SLCE.n_salcs(mbb) == 10
        @test length(mba.marked_atoms) == 128

        ds = SCEFitting.MomentDataset(mba, data; gate_eps = 2.2e-4)
        rep = Dict(t.orbit => t for t in ds.orbit_report)
        @test sort(collect(keys(rep))) == [1, 65]
        # CHANGE DETECTORS (census of the archived data under the mode rule;
        # captured from the upstream fit_ferh.jl log, 2026-08-20): recapture with
        # the reason recorded if the data or the antiparallel definition change
        @test (rep[1].n_kept, rep[1].n_rows, rep[1].n_anti) == (7744, 7744, 0)
        @test (rep[65].n_kept, rep[65].n_rows, rep[65].n_anti) == (6866, 7744, 3833)
        # hand recount of the census from the raw data: ê·e_MW < 0 over Rh atoms
        nanti = count(dot(d.constraint_axes[:, a], d.directions[:, a]) < 0
                      for d in data for a = 65:128)
        @test nanti == 3833
        f = SCEFitting.fit(SCEFitting.MomentFit, ds)
        sigma = std(SCEFitting.residuals(f))
        @printf("  FeRh p=10  sigma gated %.4f (ref 0.0086)  rows %d\n", sigma, count(ds.keep))
        @test abs(sigma - 0.0086) <= 0.01 * 0.0086 + 0.5e-4

        # column parity on 12 configurations with the mode-1 axes substituted
        sel = 1:10:111
        cfgs = [data[c].directions for c in sel]
        axs = [data[c].constraint_axes for c in sel]
        Xa = SCEFitting._design_moment(mba, cfgs, axs)
        Xb = SLCE._design_moment(mbb, cfgs, axs)
        w = column_parity(Xa, mba.salc_basis.keys, Xb, mbb.salc_basis.keys)
        @printf("  FeRh column parity: worst relative column deviation %.2e\n", w)
    end
end
