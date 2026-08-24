# The moment channel's dataset + fit (src/fitting/momentfit.jl): target assembly
# under the mode rule (hand-arithmetic oracle), mode-1 ≡ mode-4 identity
# equivalence, the decomposability gate (|M| sin²θ, the |M| = 0 pass, gated vs
# ungated disclosure, the antiparallel census), mode-1 zero-axis exclusion,
# per-orbit survival (single- and multi-orbit) + the coverage-floor refusal, the
# loud requirement checks, dataset-level time-reversal covariance (bitwise, both
# modes), the per-orbit μ₀ intercept (residual invariance under a per-orbit
# constant shift), the vanishing-column freeze, the predict_moment door, and the
# persistence refusal. Door tests added after the M4 review panel: the
# `directions` unit door and the zero-moment placeholder door at `MomentDataset`,
# a literal hand-worked row (ê = ẑ / x̂), the multi-orbit two-constant μ₀ shift, the
# normal-equation residual oracle, the Ge-only marked fixture with mixed modes, the
# exact survival ratio + tight floor refusal, and a `groups` spy estimator.
# M5 (upstream bd03517 / a6a6386 / 3d54abb): the band-profile diagnostic, the
# mark-class `salc_groups(mb)` + `GroupAdaptiveRidge(mb)` + `_reduce_to_active`
# freeze reduction, the real face-(a) vanishing fixture end-to-end (closes the
# M4 saboteur's "vanishing wiring never exercised" item), and the local-field /
# coverage / nested simple-floor diagnostics on a hand-derivable tied cell.
#
# NOTE: the planted-model "recovery" gates the mask/solve/predict PLUMBING — the
# data are generated through _design_moment itself. The design matrix has its own
# independent oracles in test_momentbasis.jl.
#
# Reuses the FeGe B20 fixture (_mb_fege / _MBFixedSG / _mb_unit) from
# test_momentbasis.jl — runtests.jl includes that file first.

using Test
using SCEFitting
using SCEFitting: _design_moment, n_salcs
using LinearAlgebra
using Random
using Statistics: quantile, std
using StaticArrays
using SCEFitting: _assemble_spacegroup

isdefined(@__MODULE__, :_mb_fege) || include("test_momentbasis.jl")

# A datum at the reference with exact directions (the direct form, so `e` is
# stored bit for bit — the moments form would re-normalize by an ulp)
function _mf_datum(e::Matrix{Float64}; M::Matrix{Float64}, mode::Int,
                   axes::Union{Matrix{Float64},Nothing} = nothing,
                   mag::Vector{Float64} = ones(size(e, 2)))
    nat = size(e, 2)
    return SpinDatum(0.0, e, mag, zeros(3, nat), zeros(3, nat), M, axes, mode)
end

# An estimator that records the `groups` it is handed and solves by plain
# least squares — gates that `fit(MomentFit, …)` forwards the per-row
# configuration labels of exactly the rows it solves on.
mutable struct _GroupSpy <: SCEFitting.AbstractEstimator
    seen::Vector{Vector{Int}}
end
_GroupSpy() = _GroupSpy(Vector{Int}[])
function SCEFitting.solve_coefficients(est::_GroupSpy, X::AbstractMatrix,
                                       y::AbstractVector; groups = nothing)
    groups === nothing && error("no groups forwarded")
    push!(est.seen, collect(Int, groups))
    return X \ y
end

# Refusal paths: the resolvability gate runs FIRST in the ctor (structural, basis
# only), so a datum-door refusal on this fixture is preceded by its legitimate
# tie-dependency warning; silence it where only the throw is under test.
_quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
_mf_refused(mbx, data; kw...) =
    @test_throws ArgumentError _quiet(() -> MomentDataset(mbx, data; kw...))

# The FeGe primitive pointed design at the fixture spec carries one structural
# tie dependency (columns 2–3, named by moment_resolvability and warned by the
# MomentDataset ctor), so every construction warns "structurally dependent" and
# every OLS solve warns rank-deficient. The helpers assert exactly that (two
# solves per fit: gated + ungated) and keep the suite quiet; match_mode :any lets
# the drop @info pass through where it also fires.
_mf_ds(mbx, data; kw...) =
    @test_logs (:warn, r"structurally dependent") match_mode = :any MomentDataset(
        mbx, data; kw...)
_mf_fit(ds) = @test_logs (:warn, r"rank deficient") (:warn, r"rank deficient") fit(
    MomentFit, ds)

@testset "moment dataset + fit" begin
    xt, sg = _mb_fege()
    bk = _MBFixedSG(sg)
    nat = 8
    spec = MomentSpec(; lmax_env = [1, 1], sampled = [true, true], lmax_mark = 1,
                      nbody = 2, cutoff_pair = 3.0, marked = [true, false])
    mb = MomentBasis(xt, spec; backend = bk)
    p = n_salcs(mb)
    marked = mb.marked_atoms
    nm = length(marked)
    @test nm == 4                      # the Fe 4a orbit; Ge is unmarked
    rng = MersenneTwister(2026)
    beta = randn(rng, p)               # the planted pointed model

    # Decomposable synthetic data from the planted model: y = X·β exactly, M built
    # back as y·ê (+ optional transverse remainder t and target corruption δ on the
    # `dirty` configs — a dirty row has g ≈ |t| and a target the model does NOT fit,
    # so the gate's keep/reject decision is observable in the residuals).
    function _mf_data(ncfg; mode::Int = 4, ndirty::Int = 0, tperp::Float64 = 0.3,
                      delta::Float64 = 0.5)
        data = SpinDatum[]
        for c = 1:ncfg
            e = _mb_unit(rng, nat)
            yv = _design_moment(mb, [e], [e]) * beta
            M = zeros(3, nat)
            for (ai, a) in enumerate(marked)
                y0 = yv[ai]
                if c <= ndirty
                    y0 += delta * randn(rng)
                    t = randn(rng, 3)
                    t .-= dot(t, e[:, a]) .* e[:, a]
                    t .*= tperp / norm(t)
                    M[:, a] .= y0 .* e[:, a] .+ t
                else
                    M[:, a] .= y0 .* e[:, a]
                end
            end
            push!(data, _mf_datum(e; M, mode, axes = mode == 1 ? copy(e) : nothing))
        end
        return data
    end

    @testset "target assembly (hand oracle) + mode rule" begin
        e = _mb_unit(rng, nat)
        M = randn(rng, 3, nat)
        ax = _mb_unit(rng, nat)
        d4 = _mf_datum(e; M, mode = 4)
        d1 = _mf_datum(e; M, mode = 1, axes = copy(ax))
        ds = _mf_ds(mb, [d4, d1]; gate_eps = 1e6, coverage_floor = 0.0)
        @test size(ds.X) == (2nm, p)
        # `order` reads the DIRECTIONS, never the axes: both data share `e`
        @test ds.order[1] == ds.order[2]
        @test ds.order[1] ≈ norm(sum(e[:, a] for a in marked)) / nm rtol = 1e-12
        for (ai, a) in enumerate(marked)
            # mode 4: ê = the datum's direction; mode 1: ê = the constraint axis —
            # both targets and gates recomputed here by hand arithmetic
            y4 = e[1, a] * M[1, a] + e[2, a] * M[2, a] + e[3, a] * M[3, a]
            y1 = ax[1, a] * M[1, a] + ax[2, a] * M[2, a] + ax[3, a] * M[3, a]
            mm = norm(M[:, a])
            @test ds.y[ai] == y4
            @test ds.y[nm + ai] == y1
            # independent gate oracle: |M| sin²θ via the angle, not the code's form
            th4 = acos(clamp(y4 / mm, -1, 1))
            th1 = acos(clamp(y1 / mm, -1, 1))
            @test ds.gate[ai] ≈ mm * sin(th4)^2 rtol = 1e-9
            @test ds.gate[nm + ai] ≈ mm * sin(th1)^2 rtol = 1e-9
            @test ds.row_config[ai] == 1 && ds.row_atom[ai] == a
            @test ds.row_config[nm + ai] == 2 && ds.row_atom[nm + ai] == a
        end
        # the design rows are exactly the marked-column-substituted evaluations
        @test ds.X == vcat(_design_moment(mb, [e], [e]),
                           _design_moment(mb, [e], [ax]))
        # everything decomposes at a huge eps; nothing is undefined
        @test all(ds.defined) && all(ds.keep)
        # literal rows worked by hand (mode 1, axis = a coordinate unit vector):
        #   ê = ẑ, M = (0.1, 0.2, 3.0): y = 3.0, |M⊥|² = 0.05, |M| = √9.05,
        #                                g = 0.05/√9.05
        #   ê = x̂, M = (−1.5, 0.4, 0.0): y = −1.5, |M⊥|² = 0.16, |M| = √2.41,
        #                                g = 0.16/√2.41
        axl = copy(ax)
        Ml = copy(M)
        axl[:, marked[1]] .= [0.0, 0.0, 1.0]; Ml[:, marked[1]] .= [0.1, 0.2, 3.0]
        axl[:, marked[2]] .= [1.0, 0.0, 0.0]; Ml[:, marked[2]] .= [-1.5, 0.4, 0.0]
        dl = _mf_ds(mb, [_mf_datum(e; M = Ml, mode = 1, axes = axl)];
                    gate_eps = 1e6, coverage_floor = 0.0)
        @test dl.y[1] == 3.0
        @test dl.gate[1] ≈ 0.05 / sqrt(9.05) rtol = 1e-12
        @test dl.y[2] == -1.5
        @test dl.gate[2] ≈ 0.16 / sqrt(2.41) rtol = 1e-12
    end

    @testset "the directions door + the zero-moment placeholder door" begin
        e = _mb_unit(rng, nat)
        M = randn(rng, 3, nat)
        # a SpinDatum built field-by-field carries no direction check: the dataset
        # constructor is the door. Off-norm column (both modes), a near-pole column
        # inside no band, an exactly-zero mode-4 column (never "undefined-axis" —
        # refused), and a non-finite column.
        for (what, col) in (("scaled", 1.5 .* e[:, marked[1]]),
                            ("pole", [0.0, 0.0, 1.0 + 5e-7]),
                            ("zero", [0.0, 0.0, 0.0]),
                            ("nan", [NaN, 0.0, 0.0]))
            ebad = copy(e); ebad[:, marked[1]] .= col
            _mf_refused(mb, [_mf_datum(ebad; M, mode = 4)]; gate_eps = 1e6,
                        coverage_floor = 0.0)
            @test_throws ArgumentError _quiet(() -> MomentDataset(
                mb, [_mf_datum(ebad; M, mode = 1, axes = copy(e))]; gate_eps = 1e6,
                coverage_floor = 0.0))
            # an off-norm ENVIRONMENT (unmarked, sampled Ge) column is refused too —
            # it feeds every neighbouring row's harmonics
            ge = first(setdiff(1:nat, marked))
            ebad2 = copy(e); ebad2[:, ge] .= col
            _mf_refused(mb, [_mf_datum(ebad2; M, mode = 4)]; gate_eps = 1e6,
                        coverage_floor = 0.0)
        end
        # the zero-moment door reads ‖MW‖ (magmoms), never M_int: a referenced atom
        # (marked Fe, or sampled-environment Ge) at ‖MW‖ ≤ atol is refused by name,
        # an unreferenced atom is not, and the atol is the caller's
        ge = first(setdiff(1:nat, marked))
        for a in (marked[2], ge)
            mag = ones(nat); mag[a] = 1e-11
            err = try
                _quiet(() -> MomentDataset(mb, [_mf_datum(e; M, mode = 4, mag)];
                                           gate_eps = 1e6, coverage_floor = 0.0))
                nothing
            catch ex
                ex
            end
            @test err isa ArgumentError && occursin("atom $a", err.msg) &&
                  occursin("zero magnetic moment", err.msg)
            # a tighter atol admits the same datum
            ds = _mf_ds(mb, [_mf_datum(e; M, mode = 4, mag)]; gate_eps = 1e6,
                        coverage_floor = 0.0, zero_moment_atol = 0.0)
            @test all(ds.keep)
        end
        # Ge unsampled and unmarked → unreferenced: a quenched Ge moment passes
        spec_fe = MomentSpec(; lmax_env = [1, 0], sampled = [true, false],
                             lmax_mark = 1, nbody = 2, cutoff_pair = 3.0,
                             marked = [true, false])
        mb_fe = MomentBasis(xt, spec_fe; backend = bk)
        @test !SCEFitting._referenced_atoms(mb_fe)[ge]
        @test all(SCEFitting._referenced_atoms(mb)[[marked; ge]])
        mag = ones(nat); mag[ge] = 0.0
        ds = _mf_ds(mb_fe, [_mf_datum(e; M, mode = 4, mag)]; gate_eps = 1e6,
                    coverage_floor = 0.0)
        @test all(ds.keep)
        _mf_refused(mb, [_mf_datum(e; M, mode = 4)]; gate_eps = 1e6,
                    zero_moment_atol = -1.0)
        # an EXACT placeholder (‖MW‖ = 0) is refused even at zero_moment_atol = 0:
        # the rule is ‖MW‖ > atol, so the boundary is a refusal
        mag0 = ones(nat); mag0[ge] = 0.0
        _mf_refused(mb, [_mf_datum(e; M, mode = 4, mag = mag0)]; gate_eps = 1e6,
                    coverage_floor = 0.0, zero_moment_atol = 0.0)
        # M_int = 0 on a marked atom is NOT the placeholder case (‖MW‖ = 1 here):
        # the quenched-passes convention stands
        Mq = copy(M); Mq[:, marked[1]] .= 0.0
        dq = _mf_ds(mb, [_mf_datum(e; M = Mq, mode = 4)]; gate_eps = 0.0,
                    coverage_floor = 0.0)
        @test dq.keep[1] && dq.y[1] == 0.0
    end

    @testset "mode-1 with axes = directions ≡ mode-4 (identity substitution)" begin
        e = _mb_unit(rng, nat)
        M = randn(rng, 3, nat)
        ds4 = _mf_ds(mb, [_mf_datum(e; M, mode = 4)];
                     gate_eps = 1e6, coverage_floor = 0.0)
        ds1 = _mf_ds(mb, [_mf_datum(e; M, mode = 1, axes = copy(e))];
                     gate_eps = 1e6, coverage_floor = 0.0)
        @test ds1.X == ds4.X
        @test ds1.y == ds4.y
    end

    @testset "planted-model recovery (prediction equivalence)" begin
        data = _mf_data(40)
        ds = _mf_ds(mb, data; gate_eps = 1e-8)
        @test all(ds.keep)
        # the OLS rank warning asserted by _mf_fit has a name: the structural
        # resolvability rank of this cell is 4 of 5 columns
        mr = moment_resolvability(mb)
        @test mr.rank == 4
        @test mr.rank + length(mr.null_combinations) == length(mr.kept)
        f = _mf_fit(ds)
        @test norm(residuals(f)) < 1e-7
        @test rmse_moment(f) < 1e-8
        # `groups` forwarded to the estimator: the per-row configuration labels of
        # exactly the solved rows, gated solve first then the ungated one
        spy = _GroupSpy()
        fspy = fit(MomentFit, ds, spy)
        @test length(spy.seen) == 2
        @test spy.seen[1] == ds.row_config[ds.keep]
        @test spy.seen[2] == ds.row_config[ds.defined]
        @test fspy.estimator === spy
        # held-out configuration: the fitted model reproduces the planted model's
        # predictions (rank-robust — any flat direction of the training design is
        # structural and annihilates the held-out rows too)
        eh = _mb_unit(rng, nat)
        model = MomentModel(f)
        @test predict_moment(model, eh) ≈ _design_moment(mb, [eh], [eh]) * beta atol =
            1e-7
        @test predict_moment(f, eh) == predict_moment(model, eh)
        @test coef(f) === f.coeffs
        @test occursin("MomentFit(", sprint(show, f))
        @test occursin("MomentModel(", sprint(show, model))
        @test occursin("MomentDataset(", sprint(show, ds))
    end

    @testset "decomposability gate: keep/reject + disclosure + |M| = 0 pass" begin
        data = _mf_data(40; ndirty = 10)
        ds = @test_logs (:warn, r"structurally dependent") (:info, r"excluded") match_mode =
            :any MomentDataset(mb, data; gate_eps = 1e-6)
        # exactly the dirty configs' rows are rejected (g ≈ tperp ≫ eps on dirty,
        # ≈ 0 on clean), and none are undefined
        @test all(ds.defined)
        @test count(ds.keep) == 30 * nm
        @test all(r -> ds.keep[r] == (ds.row_config[r] > 10), 1:length(ds.keep))
        t = only(ds.orbit_report)
        @test t.n_rows == 40 * nm && t.n_kept == 30 * nm && t.survival == 0.75
        @test t.atoms == marked
        @test t.n_anti == 0                      # mode 4: axis IS the direction
        # hand oracle for the disclosed transverse rms: 40 of 160 defined rows
        # carry |M⊥| = 0.3 exactly (the dirty construction), the rest ≈ 0
        @test t.mperp_rms ≈ sqrt(40 * 0.3^2 / 160) rtol = 1e-9
        # the gate is doing physics: the gated solve fits exactly, the ungated one
        # is polluted by the corrupted dirty targets
        f = _mf_fit(ds)
        @test rmse_moment(f) < 1e-8
        @test rmse_moment(f; gated = false) > 1e-2
        @test f.coeffs != f.coeffs_ungated
        # the ungated solve IS a least-squares solve on the defined rows: its
        # residual is orthogonal to the design (normal equations), and the rmse
        # sits inside the residual magnitudes
        ru = residuals(f; gated = false)
        @test length(ru) == 40 * nm
        @test norm(ds.X[ds.defined, :]' * ru) < 1e-9 * norm(ds.y[ds.defined])
        @test minimum(abs, ru) <= rmse_moment(f; gated = false) <= maximum(abs, ru)
        @test ru[1:10nm] == ds.y[1:10nm] .- ds.X[1:10nm, :] * f.coeffs_ungated
        # a quenched moment (M = 0) has g = 0 and passes with target 0
        e = _mb_unit(rng, nat)
        M = randn(rng, 3, nat)
        M[:, marked[1]] .= 0.0
        dq = @test_logs (:warn, r"structurally dependent") (:info, r"excluded") match_mode =
            :any MomentDataset(mb, [_mf_datum(e; M, mode = 4)];
                               gate_eps = 0.0, coverage_floor = 0.0)
        @test dq.keep[1] && dq.y[1] == 0.0 && dq.gate[1] == 0.0
    end

    @testset "mode-1 zero-axis rows are excluded, not fitted" begin
        data = _mf_data(20; mode = 1)
        # unconstrain atom marked[2] in the first 4 configs
        for c = 1:4
            data[c].constraint_axes[:, marked[2]] .= 0.0
        end
        ds = @test_logs (:warn, r"structurally dependent") (:info,
            r"undefined-axis") match_mode = :any MomentDataset(mb, data;
                                                               gate_eps = 1e-8)
        bad = findall(r -> !ds.defined[r], 1:length(ds.defined))
        @test length(bad) == 4
        @test all(r -> ds.row_atom[r] == marked[2] && ds.row_config[r] <= 4, bad)
        @test all(isnan, ds.y[bad]) && all(isnan, ds.gate[bad])
        @test count(ds.keep) == 20nm - 4
        t = only(ds.orbit_report)
        @test t.n_defined == 20nm - 4 && t.n_kept == 20nm - 4
        # the rms transverse remainder runs over the DEFINED rows only (a NaN
        # from an undefined row would poison it); planted decomposable ⇒ ≈ 0
        @test t.mperp_rms ≈ 0 atol = 1e-12
        # band profile on a PARTIALLY kept configuration: configs 1–4 keep 3 of
        # 4 rows; shift every marked M by c0·ê and profile with the unshifted
        # model — every per-config mean residual must be c0 (a fixed divisor nm
        # would give 0.75·c0 on configs 1–4)
        c0 = 0.41
        shifted = SpinDatum[]
        for d in data
            M = copy(d.moments_bare)
            for a in marked
                M[:, a] .+= c0 .* d.constraint_axes[:, a]
            end
            push!(shifted, _mf_datum(d.directions; M, mode = 1,
                                     axes = copy(d.constraint_axes)))
        end
        f0 = _mf_fit(ds)
        dss = @test_logs (:warn, r"structurally dependent") (:info,
            r"undefined-axis") match_mode = :any MomentDataset(mb, shifted;
                                                               gate_eps = 1e-8)
        profs = moment_band_profile(MomentModel(f0), dss)
        @test length(profs.mean_residual) == 20
        @test all(isapprox.(profs.mean_residual, c0; atol = 1e-7))
        # survival counts the undefined rows in its denominator …
        @test t.survival == (20nm - 4) / (20nm)
        # … so a floor just above it refuses, naming the undefined-axis count
        err = try
            _quiet(() -> MomentDataset(mb, data; gate_eps = 1e-8, coverage_floor = 0.96))
            nothing
        catch ex
            ex
        end
        @test err isa ArgumentError && occursin("4 undefined-axis", err.msg)
        # the placeholder design row is the identity-substituted evaluation
        e1 = data[1].directions
        @test ds.X[bad[1], :] == vec(_design_moment(mb, [e1], [e1])[
            findfirst(==(marked[2]), marked), :])
        # both fits run on the defined rows only, and still recover the model
        f = _mf_fit(ds)
        @test length(residuals(f)) == 20nm - 4
        @test length(residuals(f; gated = false)) == 20nm - 4
        @test rmse_moment(f) < 1e-8
        # the caller's datum was not mutated by the placeholder write
        @test all(iszero, data[1].constraint_axes[:, marked[2]])
    end

    @testset "coverage floor + loud requirements" begin
        # 60 % dirty at a tight eps → survival 0.4 < floor 0.5 → refuse
        data = _mf_data(20; ndirty = 12)
        @test_throws ArgumentError _quiet(() -> MomentDataset(mb, data; gate_eps = 1e-6))
        # ... but an explicitly lowered floor accepts the same data
        ds = @test_logs (:warn, r"structurally dependent") (:info, r"excluded") match_mode =
            :any MomentDataset(mb, data; gate_eps = 1e-6, coverage_floor = 0.3)
        @test count(ds.keep) == 8 * nm
        # requirement checks
        e = _mb_unit(rng, nat)
        M = randn(rng, 3, nat)
        _mf_refused(mb, SpinDatum[]; gate_eps = 1e-6)
        @test_throws UndefKeywordError MomentDataset(mb, [_mf_datum(e; M, mode = 4)])
        _mf_refused(mb, [_mf_datum(e; M, mode = 4)]; gate_eps = -1.0)
        _mf_refused(mb, [_mf_datum(e; M, mode = 4)]; gate_eps = 1e-6,
                    coverage_floor = 1.5)
        no_m = SpinDatum(0.0, e, ones(nat), zeros(3, nat), zeros(3, nat),
                         nothing, nothing, 4)
        @test_throws ArgumentError _quiet(() -> MomentDataset(mb, [no_m]; gate_eps = 1e-6))
        no_mode = SpinDatum(0.0, e, ones(nat), zeros(3, nat), zeros(3, nat),
                            M, nothing, nothing)
        _mf_refused(mb, [no_mode]; gate_eps = 1e-6)
        # atom-count mismatch dies at the door
        e4 = _mb_unit(rng, 4)
        small = _mf_datum(e4; M = randn(rng, 3, 4), mode = 4)
        @test_throws ArgumentError _quiet(() -> MomentDataset(mb, [small]; gate_eps = 1e-6))
        # gate_eps = 0 with a floor of 0 constructs an empty-keep dataset; the fit
        # door refuses it (the ctor's floor would otherwise have fired first)
        dirty = _mf_data(5; ndirty = 5)
        ds0 = @test_logs (:warn, r"structurally dependent") (:info,
            r"excluded") match_mode = :any MomentDataset(mb, dirty; gate_eps = 0.0,
                                                         coverage_floor = 0.0)
        @test !any(ds0.keep)
        @test_throws ArgumentError fit(MomentFit, ds0)
    end

    @testset "time reversal at the dataset level (bitwise)" begin
        data = _mf_data(6)
        flipped = [_mf_datum(-d.directions; M = -d.moments_bare, mode = 4)
                   for d in data]
        ds = _mf_ds(mb, data; gate_eps = 1e-8)
        dsf = _mf_ds(mb, flipped; gate_eps = 1e-8)
        # even-Σl labels (mark rank counts): every design entry is TR-even exactly,
        # and y = ê·M is invariant under (ê, M) → (−ê, −M)
        @test dsf.X == ds.X
        @test dsf.y == ds.y
        # same statement for mode 1 (axes flip with everything else)
        data1 = _mf_data(4; mode = 1)
        flip1 = [_mf_datum(-d.directions; M = -d.moments_bare, mode = 1,
                           axes = -d.constraint_axes) for d in data1]
        dsa = _mf_ds(mb, data1; gate_eps = 1e-8)
        dsb = _mf_ds(mb, flip1; gate_eps = 1e-8)
        @test dsb.X == dsa.X
        @test dsb.y == dsa.y
    end

    @testset "per-orbit μ₀ intercept absorbs a constant shift" begin
        data = _mf_data(30)
        shifted = SpinDatum[]
        c0 = 0.37
        for d in data
            M = copy(d.moments_bare)
            for a in marked
                M[:, a] .+= c0 .* d.directions[:, a]     # y → y + c0, M⊥ unchanged
            end
            push!(shifted, _mf_datum(d.directions; M, mode = 4))
        end
        f = _mf_fit(_mf_ds(mb, data; gate_eps = 1e-8))
        fs = _mf_fit(_mf_ds(mb, shifted; gate_eps = 1e-8))
        # the l = 0 1-body [MARK] μ₀ column of the (single) Fe orbit absorbs the
        # shift: residuals are invariant (which coefficients move is not asserted)
        @test maximum(abs, residuals(fs) - residuals(f)) < 1e-9
        # and the shifted targets really did move by c0 on every kept row
        @test all(fs.dataset.y .- f.dataset.y .≈ c0)
    end

    @testset "per-orbit μ₀ (two orbits): distinct shifts, one intercept each" begin
        # Fe and Ge both marked: a shift c0 on every Fe row and c1 ≠ c0 on every Ge
        # row must be absorbed with invariant residuals — a single SHARED intercept
        # could not absorb two different constants
        spec2 = MomentSpec(; lmax_env = [1, 1], sampled = [true, true],
                           lmax_mark = 1, nbody = 2, cutoff_pair = 3.0,
                           marked = [true, true])
        mb2 = MomentBasis(xt, spec2; backend = bk)
        m2 = mb2.marked_atoms
        beta2 = randn(rng, n_salcs(mb2))
        data2 = SpinDatum[]
        shifted2 = SpinDatum[]
        c0, c1 = 0.37, -0.81
        for c = 1:30
            e = _mb_unit(rng, nat)
            yv = _design_moment(mb2, [e], [e]) * beta2
            M = zeros(3, nat)
            Ms = zeros(3, nat)
            for (ai, a) in enumerate(m2)
                M[:, a] .= yv[ai] .* e[:, a]
                Ms[:, a] .= (yv[ai] + (a <= 4 ? c0 : c1)) .* e[:, a]
            end
            push!(data2, _mf_datum(e; M, mode = 4))
            push!(shifted2, _mf_datum(e; M = Ms, mode = 4))
        end
        f = _mf_fit(_mf_ds(mb2, data2; gate_eps = 1e-8))
        fs = _mf_fit(_mf_ds(mb2, shifted2; gate_eps = 1e-8))
        @test maximum(abs, residuals(fs) - residuals(f)) < 1e-9
        dsy = fs.dataset.y .- f.dataset.y
        fe_rows = fs.dataset.orbit_rep .== 1
        @test all(dsy[fe_rows] .≈ c0) && all(dsy[.!fe_rows] .≈ c1)
        # the two l = 0 [MARK] columns are orbit indicators: constant 1 on their
        # own orbit's rows, exactly 0 elsewhere
        mu = findall(k -> k.body == 1 && k.L_S == 0, mb2.salc_basis.keys)
        @test length(mu) == 2
        X = f.dataset.X
        for (j, o) in zip(mu, (1, 5))
            rows = f.dataset.orbit_rep .== o
            @test all(X[rows, j] .== X[findfirst(rows), j]) && X[findfirst(rows), j] != 0
            @test all(iszero, X[.!rows, j])
        end
    end

    @testset "Ge-only marked basis with mixed modes (row/atom layout)" begin
        # marked = [false, true]: rows are the Ge atoms 5:8 only; the Fe atoms are
        # environment. A mode-4 and a mode-1 datum in one dataset; the mode-1
        # datum unconstrains atom 7 → exactly row nm + 3 is undefined.
        specg = MomentSpec(; lmax_env = [1, 1], sampled = [true, true],
                           lmax_mark = 1, nbody = 2, cutoff_pair = 3.0,
                           marked = [false, true])
        mbg = MomentBasis(xt, specg; backend = bk)
        @test mbg.marked_atoms == [5, 6, 7, 8]
        e4 = _mb_unit(rng, nat); M4 = randn(rng, 3, nat)
        e1 = _mb_unit(rng, nat); M1 = randn(rng, 3, nat)
        ax1 = _mb_unit(rng, nat); ax1[:, 7] .= 0.0
        ds = @test_logs (:warn, r"structurally dependent") (:info,
            r"undefined-axis") match_mode = :any MomentDataset(
            mbg, [_mf_datum(e4; M = M4, mode = 4),
                  _mf_datum(e1; M = M1, mode = 1, axes = ax1)];
            gate_eps = 1e6, coverage_floor = 0.0)
        @test ds.row_atom == repeat(5:8, 2)
        @test ds.row_config == repeat(1:2, inner = 4)
        @test findall(!, ds.defined) == [4 + 3]
        @test isnan(ds.y[7])
        @test ds.y[3] == dot(e4[:, 7], M4[:, 7])
        @test ds.y[8] == dot(ax1[:, 8], M1[:, 8])
        @test only(ds.orbit_report).orbit == 5
        @test only(ds.orbit_report).atoms == [5, 6, 7, 8]
        @test only(ds.orbit_report).n_defined == 7
        # the Fe columns of `directions` are environment coordinates on every row:
        # perturb one Fe direction and every Ge row of that config moves (the
        # Fe sites are within the 3.0 Å pair cutoff of each Ge), no other row does
        e4b = copy(e4); e4b[:, 1] .= _mb_unit(rng, 1)[:, 1]
        dsb = _mf_ds(mbg, [_mf_datum(e4b; M = M4, mode = 4),
                           _mf_datum(e1; M = M1, mode = 1, axes = ax1)];
                     gate_eps = 1e6, coverage_floor = 0.0)
        @test dsb.X[5:8, :] == ds.X[5:8, :]
        @test all(r -> dsb.X[r, :] != ds.X[r, :], 1:4)
    end

    @testset "predict_moment plumbing + the door" begin
        f = _mf_fit(_mf_ds(mb, _mf_data(12); gate_eps = 1e-8))
        model = MomentModel(f)
        e = _mb_unit(rng, nat)
        ax = _mb_unit(rng, nat)
        # explicit axes reproduce the training-row evaluation; default is axes = e
        @test predict_moment(model, e; axes = ax) ==
              vec(_design_moment(mb, [e], [ax]) * model.coeffs)
        @test predict_moment(model, e) == predict_moment(model, e; axes = e)
        @test_throws DimensionMismatch predict_moment(model, e;
                                                      axes = _mb_unit(rng, 4))
        @test_throws DimensionMismatch MomentModel(mb, randn(rng, p + 1))
        # the DOOR: off-norm e is refused; off-norm MARKED axes column is refused;
        # a marked axes column with a component > 1 inside the norm band is refused;
        # unmarked axes columns are free and provably never read
        ebad = copy(e); ebad[:, 1] .*= 1.5
        @test_throws ArgumentError predict_moment(model, ebad)
        axbad = copy(e); axbad[:, marked[1]] .*= 1.5
        @test_throws ArgumentError predict_moment(model, e; axes = axbad)
        axpole = copy(e); axpole[:, marked[1]] .= [0.0, 0.0, 1.0 + 5e-7]
        @test_throws ArgumentError predict_moment(model, e; axes = axpole)
        ge = first(setdiff(1:nat, marked))
        ax2 = copy(e); ax2[:, ge] .= [10.0, -3.0, 0.5]
        @test predict_moment(model, e; axes = ax2) == predict_moment(model, e)
        # a bad `e` is refused even when perfectly good axes are supplied (the
        # environment is read from `e`, never from `axes`)
        @test_throws ArgumentError predict_moment(model, ebad; axes = e)
        # closed-form component readout: ê = ẑ on every marked column returns the
        # z-component of the model moment, i.e. the same design evaluated with
        # substituted ẑ columns
        axz = copy(e); for a in marked; axz[:, a] .= [0.0, 0.0, 1.0]; end
        @test predict_moment(model, e; axes = axz) ==
              vec(_design_moment(mb, [e], [axz]) * model.coeffs)
    end

    @testset "resolvability wiring + vanishing-column freeze" begin
        ds = _mf_ds(mb, _mf_data(10); gate_eps = 1e-8)
        mr = moment_resolvability(mb)
        @test ds.vanishing == mr.vanishing == Int[]
        @test ds.dependent == mr.null_combinations
        @test length(only(ds.dependent)) == 2      # the columns-2/3 tie
        # freeze MECHANISM (no vanishing column exists on this fixture, so drive
        # the branch directly through the inner ctor): fit must return exact 0.0
        # for the frozen column and the reduced solve for the rest
        ds2 = MomentDataset(ds.basis, ds.X, ds.y, ds.defined, ds.keep, ds.gate,
                            ds.row_config, ds.row_atom, ds.orbit_rep,
                            ds.orbit_report, ds.order, [2], ds.dependent,
                            ds.gate_eps, ds.coverage_floor)
        f2 = fit(MomentFit, ds2)     # dropping col 2 removes the tie: no warning
        @test f2.coeffs[2] == 0.0 && f2.coeffs_ungated[2] == 0.0
        act = [1, 3, 4, 5]
        @test f2.coeffs[act] ≈ ds.X[ds.keep, act] \ ds.y[ds.keep] rtol = 1e-10
    end

    @testset "antiparallel census (mode-1 gauge rows)" begin
        data = _mf_data(6; mode = 1)
        for c = 1:3
            data[c].constraint_axes[:, marked[1]] .*= -1.0   # re-gauged axis
        end
        ds = @test_logs (:warn, r"structurally dependent") (:info,
            r"antiparallel") match_mode = :any MomentDataset(mb, data; gate_eps = 1e-8)
        t = only(ds.orbit_report)
        @test t.n_anti == 3
        @test all(ds.keep)          # the gate is even in y: gauge rows still pass
        # `order` is mode-independent: the same directions as mode 4 give the
        # same coordinate (the re-gauged axes do not enter)
        ds4 = _mf_ds(mb, [_mf_datum(d.directions; M = d.moments_bare, mode = 4)
                          for d in data]; gate_eps = 1e-8)
        @test ds4.order == ds.order
    end

    @testset "multi-orbit bookkeeping (Fe + Ge marked)" begin
        spec2 = MomentSpec(; lmax_env = [1, 1], sampled = [true, true],
                           lmax_mark = 1, nbody = 2, cutoff_pair = 3.0,
                           marked = [true, true])
        mb2 = MomentBasis(xt, spec2; backend = bk)
        m2 = mb2.marked_atoms
        @test m2 == collect(1:8)
        e = _mb_unit(rng, nat)
        # Fe rows decomposable, Ge rows carry a large transverse remainder
        M = zeros(3, nat)
        for a = 1:4
            M[:, a] .= 1.3 .* e[:, a]
        end
        for a = 5:8
            t = randn(rng, 3); t .-= dot(t, e[:, a]) .* e[:, a]
            M[:, a] .= 0.2 .* e[:, a] .+ 0.4 .* t ./ norm(t)
        end
        d = _mf_datum(e; M, mode = 4)
        # per-orbit disaggregation: Fe orbit (rep 1) clean, Ge orbit (rep 5) dirty
        ds = @test_logs (:warn, r"structurally dependent") (:info, r"excluded") match_mode =
            :any MomentDataset(mb2, [d]; gate_eps = 1e-6, coverage_floor = 0.0)
        @test [t.orbit for t in ds.orbit_report] == [1, 5]
        @test ds.orbit_report[1].atoms == [1, 2, 3, 4]
        @test ds.orbit_report[2].atoms == [5, 6, 7, 8]
        @test ds.orbit_rep == [1, 1, 1, 1, 5, 5, 5, 5]
        @test ds.orbit_report[1].n_kept == 4 && ds.orbit_report[2].n_kept == 0
        @test ds.orbit_report[2].mperp_rms ≈ 0.4 rtol = 1e-9
        @test length(ds.dependent) == 2
        # the per-orbit coverage refusal names the failing (Ge) orbit
        err = try
            _quiet(() -> MomentDataset(mb2, [d]; gate_eps = 1e-6))
            nothing
        catch ex
            ex
        end
        @test err isa ArgumentError && occursin("orbit 5", err.msg)
    end

    @testset "the dataset door is hard on an unclassifiable basis" begin
        # the primitive-cell STAR basis keeps repeated-image environment members
        # (upstream parity); moment_resolvability refuses it, and MomentDataset
        # propagates that refusal instead of building a design it cannot vouch for
        sstar = MomentSpec(; lmax_env = [2, 2], sampled = [true, true], lmax_mark = 2,
                           nbody = 3, cutoff_pair = 3.3, cutoff_star = 3.3)
        mstar = MomentBasis(xt, sstar; backend = bk)
        e = _mb_unit(rng, nat)
        @test_throws SCEFitting.UnclassifiableBasis _quiet(() -> MomentDataset(
            mstar, [_mf_datum(e; M = randn(rng, 3, nat), mode = 4)]; gate_eps = 1e6,
            coverage_floor = 0.0))
    end

    @testset "band-profile diagnostic" begin
        data = _mf_data(24)
        ds = _mf_ds(mb, data; gate_eps = 1e-8)
        f = _mf_fit(ds)
        # the order coordinate is the marked-sublattice |⟨e⟩| — hand oracle
        for c in (1, 24)
            m = sum(data[c].directions[:, a] for a in marked)
            @test ds.order[c] ≈ norm(m) / nm rtol = 1e-12
        end
        prof = moment_band_profile(f)
        # re-computation of the profile from public fields (row → config
        # mapping; valid here because every row is kept — the residual
        # DEFINITION is gated by the planted-slope recovery below)
        pred = ds.X * coef(f)
        mres = [sum((ds.y .- pred)[(c-1)*nm+1:c*nm]) / nm for c = 1:24]
        @test prof.mean_residual ≈ mres atol = 1e-12
        @test sum(b.n for b in prof.bands) == 24
        @test prof.bands[1].hi <= prof.bands[end].lo   # ordered bins
        # bin-free line: independent least-squares via the normal equations
        A = [ones(24) ds.order]
        ab = (A' * A) \ (A' * mres)
        @test prof.intercept ≈ ab[1] atol = 1e-10
        @test prof.slope ≈ ab[2] atol = 1e-10
        @test abs(prof.r) <= 1 + 1e-12
        # a planted linear-in-order corruption is recovered by the slope
        data2 = SpinDatum[]
        for (c, d) in enumerate(data)
            M = copy(d.moments_bare)
            for a in marked
                M[:, a] .+= (0.5 * ds.order[c]) .* d.directions[:, a]
            end
            push!(data2, _mf_datum(d.directions; M, mode = 4))
        end
        ds2 = _mf_ds(mb, data2; gate_eps = 1e-8)
        prof2 = moment_band_profile(MomentModel(f), ds2)
        @test prof2.slope ≈ prof.slope + 0.5 atol = 1e-6
        # gated rows are excluded from the profile: a config whose rows are all
        # gate-rejected disappears
        data3 = _mf_data(6; ndirty = 1)
        ds3 = @test_logs (:warn, r"structurally dependent") (:info,
            r"excluded") match_mode = :any MomentDataset(mb, data3;
                                                         gate_eps = 1e-6,
                                                         coverage_floor = 0.5)
        prof3 = moment_band_profile(MomentModel(f), ds3)
        @test length(prof3.mean_residual) == 5
        # equal-count binning on n = 5 < … not divisible by 4: every config in a
        # band, counts differ by ≤ 1, bands ordered and covering the order range
        @test sum(b.n for b in prof3.bands) == 5
        @test maximum(b.n for b in prof3.bands) - minimum(b.n for b in prof3.bands) <= 1
        @test all(b.lo <= b.hi for b in prof3.bands)
        @test all(prof3.bands[i].hi <= prof3.bands[i + 1].lo
                  for i = 1:length(prof3.bands) - 1)
        @test prof3.bands[1].lo == minimum(prof3.order)
        @test prof3.bands[end].hi == maximum(prof3.order)
        # Pearson r carries the slope's sign: a perfect linear plant has r → ±1
        @test prof2.r > 0.99
        data2m = SpinDatum[]
        for (c, d) in enumerate(data)
            M = copy(d.moments_bare)
            for a in marked
                M[:, a] .-= (0.5 * ds.order[c]) .* d.directions[:, a]
            end
            push!(data2m, _mf_datum(d.directions; M, mode = 4))
        end
        prof2m = moment_band_profile(MomentModel(f), _mf_ds(mb, data2m; gate_eps = 1e-8))
        @test prof2m.r < -0.99
        @test prof2m.slope ≈ prof.slope - 0.5 atol = 1e-6
        # the MomentFit method profiles the GATED coefficients: on a dirty
        # dataset the gated model reproduces the clean configs exactly while
        # the ungated one is polluted (test "decomposability gate")
        datad = _mf_data(40; ndirty = 10)
        dsd = @test_logs (:warn, r"structurally dependent") (:info,
            r"excluded") match_mode = :any MomentDataset(mb, datad; gate_eps = 1e-6)
        fd = _mf_fit(dsd)
        profd = moment_band_profile(fd)
        @test length(profd.mean_residual) == 30         # the 10 dirty configs are gone
        @test all(abs.(profd.mean_residual) .< 1e-7)
        @test maximum(abs, fd.coeffs_ungated - fd.coeffs) > 1e-3   # the two differ
        # planted slope recovered through a gated-out configuration too
        datads = SpinDatum[]
        for (c, d) in enumerate(datad)
            M = copy(d.moments_bare)
            for a in marked
                M[:, a] .+= (0.5 * dsd.order[c]) .* d.directions[:, a]
            end
            push!(datads, _mf_datum(d.directions; M, mode = 4))
        end
        dsds = @test_logs (:warn, r"structurally dependent") (:info,
            r"excluded") match_mode = :any MomentDataset(mb, datads; gate_eps = 1e-6)
        @test moment_band_profile(MomentModel(fd), dsds).slope ≈ 0.5 atol = 1e-6
        # degenerate order spread (the same directions twice, bitwise): slope and
        # r are NaN, the intercept is the plain mean residual — by hand
        # mres = [0, 0.2] (the second datum's targets shifted by 0.2)
        ed = data[1].directions
        M0 = data[1].moments_bare
        M2 = copy(M0)
        for a in marked
            M2[:, a] .+= 0.2 .* ed[:, a]
        end
        dsdeg = _mf_ds(mb, [_mf_datum(ed; M = M0, mode = 4),
                            _mf_datum(ed; M = M2, mode = 4)]; gate_eps = 1e-8)
        pdeg = moment_band_profile(MomentModel(f), dsdeg)
        @test isnan(pdeg.slope) && isnan(pdeg.r)
        @test pdeg.intercept ≈ 0.1 atol = 1e-7
        # equal-COUNT bins: two configs give two singleton bands at one value
        @test sum(b.n for b in pdeg.bands) == 2
        @test all(b.lo == b.hi == pdeg.order[1] for b in pdeg.bands)
        @test_throws ArgumentError moment_band_profile(f; nbins = 0)
        # fewer configs than bins: every config still lands in a band (the naive
        # div-by-nbins binning dropped the high-|⟨e⟩| end — upstream review)
        ds4 = _mf_ds(mb, _mf_data(3); gate_eps = 1e-8)
        prof4 = moment_band_profile(MomentModel(f), ds4; nbins = 4)
        @test sum(b.n for b in prof4.bands) == 3
        @test length(prof4.bands) == 3
        # a model on a different basis object is refused (same-width nonsense door)
        spec_alt = MomentSpec(; lmax_env = [1, 1], sampled = [true, true],
                              lmax_mark = 1, nbody = 2, cutoff_pair = 3.0,
                              marked = [true, false])
        mb_alt = MomentBasis(xt, spec_alt; backend = bk)
        malt = MomentModel(mb_alt, coef(f))
        @test_throws ArgumentError moment_band_profile(malt, ds3)
        # the dataset's dependency record is a COPY of the cached resolvability
        # record, never an alias into the basis cache
        @test ds.dependent == moment_resolvability(mb).null_combinations
        @test all(a !== b for (a, b) in zip(ds.dependent,
                                             moment_resolvability(mb).null_combinations))
    end

    @testset "salc_groups(mb): mark classes split, gauge blocks fold" begin
        # pair basis, both species marked: every Fe–Ge pair orbit carries one
        # Fe-marked and one Ge-marked class sharing (body, orbit_id, decors) —
        # the energy-side key folds them, the pointed key must not
        spg = MomentSpec(; lmax_env = [1, 1], sampled = [true, true],
                         lmax_mark = 1, nbody = 2, cutoff_pair = 3.0,
                         marked = [true, true])
        mbg = MomentBasis(xt, spg; backend = bk)
        gl = SCEFitting.salc_groups(mbg)
        n = n_salcs(mbg)
        ks = mbg.salc_basis.keys
        @test length(gl) == n
        @test sort(unique(gl)) == 1:maximum(gl)     # contiguous, no gaps
        # the pointed key REFINES the energy-side key
        ekey(j) = (ks[j].body, ks[j].orbit_id, ks[j].decors)
        for j = 1:n, k = (j + 1):n
            gl[j] == gl[k] && @test ekey(j) == ekey(k)
        end
        # independent (design-side) oracle: a column writes rows of exactly its
        # mark class's atoms, so two mark classes of a MIXED-species pair orbit
        # have DISJOINT row support while sharing the energy-side key — the
        # recorded defect the pointed key exists to fix
        cfgs = [_mb_unit(rng, nat) for _ = 1:3]
        X = _design_moment(mbg, cfgs, [copy(c) for c in cfgs])
        nmg = length(mbg.marked_atoms)
        sup(j) = Set(mbg.marked_atoms[a] for c = 0:2 for a = 1:nmg
                     if X[c*nmg+a, j] != 0.0)
        sups = [sup(j) for j = 1:n]
        split_disjoint = false
        for j = 1:n, k = (j + 1):n
            if ekey(j) == ekey(k) && gl[j] != gl[k] &&
               isempty(intersect(sups[j], sups[k]))
                split_disjoint = true
            end
        end
        @test split_disjoint
        # every energy-side pair-orbit key on this cell splits into ≥ 2 groups
        # (the resolvability census promises ≥ 2 mark classes per pair orbit)
        for u in unique(ekey(j) for j = 1:n if ks[j].body == 2)
            @test length(unique(gl[j] for j = 1:n if ekey(j) == u)) >= 2
        end
        # GAR convenience: labels = salc_groups, unit weights, one per group
        gar = GroupAdaptiveRidge(mbg; lambda = 1e-3)
        @test gar.column_groups == gl
        @test gar.group_weights == ones(maximum(gl))

        # gauge blocks of ONE mark class fold into ONE group (star basis:
        # multi-column groups at these caps) — same group ⇒ same row support
        sps = MomentSpec(; lmax_env = [2, 2], sampled = [true, true],
                         lmax_mark = 2, nbody = 3, cutoff_pair = 3.0,
                         cutoff_star = 3.0, marked = [true, true])
        mbs = MomentBasis(xt, sps; backend = bk)
        gls = SCEFitting.salc_groups(mbs)
        Gs = maximum(gls)
        @test any(g -> count(==(g), gls) > 1, 1:Gs)
        Xs = _design_moment(mbs, cfgs, [copy(c) for c in cfgs])
        nms = length(mbs.marked_atoms)
        sups2 = [Set(mbs.marked_atoms[a] for c = 0:2 for a = 1:nms
                     if Xs[c*nms+a, j] != 0.0) for j = 1:n_salcs(mbs)]
        for g = 1:Gs
            js = findall(==(g), gls)
            length(js) > 1 &&
                @test all(sups2[j] == sups2[js[1]] for j in js[2:end])
        end
    end

    @testset "GroupAdaptiveRidge follows the vanishing-column freeze" begin
        # inner-ctor injection (as in the freeze testset): estimator column
        # metadata must shrink with the active mask — full-basis labels, one
        # column frozen, and the fit must run with the frozen coef EXACTLY 0.0
        ds = _mf_ds(mb, _mf_data(10); gate_eps = 1e-8)
        ds2 = MomentDataset(ds.basis, ds.X, ds.y, ds.defined, ds.keep, ds.gate,
                            ds.row_config, ds.row_atom, ds.orbit_rep,
                            ds.orbit_report, ds.order, [2], ds.dependent,
                            ds.gate_eps, ds.coverage_floor)
        gar = GroupAdaptiveRidge(mb; lambda = 0.0)
        f = fit(MomentFit, ds2, gar)
        @test f.coeffs[2] == 0.0 && f.coeffs_ungated[2] == 0.0
        @test all(isfinite, f.coeffs)
        # lambda = 0: the group-adaptive solve degenerates to the unpenalized
        # solve on the active columns — same answer as the OLS freeze path
        act = [1, 3, 4, 5]
        @test f.coeffs[act] ≈ ds.X[ds.keep, act] \ ds.y[ds.keep] rtol = 1e-8
        # the stored estimator is the caller's (un-reduced) object
        @test f.estimator === gar
        # mismatched labels (built on a DIFFERENT basis) refuse loudly
        bad = GroupAdaptiveRidge([1, 2, 3], ones(3); lambda = 0.0)
        @test_throws DimensionMismatch fit(MomentFit, ds2, bad)
    end

    @testset "vanishing columns end-to-end (real face-(a) fixture)" begin
        # The pointed analog of the energy side's face (a): a four-fold WS tie
        # (Δf = (.5, .5, .25) on a 3×3×6 cell) partially fused by m_y, which
        # fixes both atoms and permutes the ±y images. With isotropy = false the
        # Lf ≠ 0 blocks pair geometry-odd member weights with spin-component
        # content, and on CELL-PERIODIC data (both images read the same atom)
        # the odd content cancels identically — 16 of 38 columns vanish. Under
        # isotropy = true no case is known: for Lf = 0 the transport matrix is
        # D⁰ = 1, so every image of an assignment carries an IDENTICAL folded
        # weight — nothing of opposite sign exists to cancel under the periodic
        # fold (an argument for single-assignment Lf = 0 blocks, not a general
        # theorem; probes found no isotropic vanishing on this cell). That is
        # why the freeze-mechanism test above needs ctor injection and THIS one
        # does not. P1 on the same crystal is the face-(b) control:
        # dependencies, never per-column cancellation.
        crv = Crystal(Lattice(Matrix(Diagonal([3.0, 3.0, 6.0]))),
                      [0.1 0.6; 0.0 0.5; 0.03 0.28], [1, 1], ["Fe"])
        I3 = SMatrix{3,3,Float64}(Matrix(1.0 * I(3)))
        MY = SMatrix{3,3,Float64}(Diagonal([1.0, -1.0, 1.0]))
        t0 = SVector{3,Float64}(0, 0, 0)
        sgv = _assemble_spacegroup(crv, [I3, MY], [t0, t0], "Pm-tie", 6; tol = 1e-5)
        spv = MomentSpec(; lmax_env = [1], sampled = [true], lmax_mark = 1,
                         nbody = 2, cutoff_pair = 2.7, isotropy = false)
        mbv = MomentBasis(crv, spv; backend = _MBFixedSG(sgv))
        mrv = moment_resolvability(mbv)
        @test !isempty(mrv.vanishing)
        @test length(mrv.vanishing) < n_salcs(mbv)
        # upstream measured 16 of 38 on this cell — a change-detecting PIN
        # (captured from SLCE.jl 3d54abb's test comment; recapture if the
        # enumeration or the screen changes), not a derivation
        @test (length(mrv.vanishing), n_salcs(mbv)) == (16, 38)

        # independent numerical oracle: on random cell-periodic configurations
        # the vanishing columns are EXACTLY zero (opposite member weights on
        # identical factors — cancellation in floating point is exact), and
        # every kept column is not
        ev = [_mb_unit(rng, 2) for _ = 1:8]
        Xv = _design_moment(mbv, ev, [copy(e) for e in ev])
        # the physical claim, tolerance form:
        @test all(j -> norm(Xv[:, j]) < 1e-12, mrv.vanishing)
        # ... and a change-detecting pin on the EXACT cancellation (paired
        # member weights are bitwise opposite and the per-column accumulation
        # order is deterministic — an implementation property, not an FP
        # theorem; recapture if the design accumulation order changes)
        @test all(j -> norm(Xv[:, j]) == 0.0, mrv.vanishing)
        @test all(j -> norm(Xv[:, j]) > 1e-8, setdiff(1:n_salcs(mbv), mrv.vanishing))

        # the null report is COMPLETE: one combination per flat direction
        @test length(mrv.null_combinations) ==
              length(mrv.kept) - mrv.rank

        # P1 control: face (b) only — dependencies, no vanishing column. This
        # control found a real defect upstream: with a WIDE signature block
        # (more kept columns than rows) the economy SVD's V lists only
        # min(r, c) directions and the report used to come back EMPTY.
        sgp = _assemble_spacegroup(crv, [I3], [t0], "P1-tie", 1; tol = 1e-5)
        mbp = MomentBasis(crv, spv; backend = _MBFixedSG(sgp))
        mrp = moment_resolvability(mbp)
        @test isempty(mrp.vanishing)
        @test length(mrp.null_combinations) == length(mrp.kept) - mrp.rank > 0
        # independent oracle: every reported combination annihilates the
        # numerical design on random cell-periodic data
        Xp = _design_moment(mbp, ev, [copy(e) for e in ev])
        for comb in mrp.null_combinations
            v = zeros(n_salcs(mbp))
            for (j, w) in comb
                v[j] = w
            end
            @test norm(Xp * v) < 1e-8 * norm(Xp) * norm(v)
        end
        # the wide block's dependency log is summarized past 8 combinations
        @test length(mrp.null_combinations) > 8
        datap = [_mf_datum(e; M = 1.2 .* e, mode = 4) for e in ev]
        dsp = @test_logs (:warn, r"structurally dependent") MomentDataset(
            mbp, datap; gate_eps = 1e-8)
        @test length(dsp.dependent) == length(mrp.null_combinations)

        # end-to-end: the dataset door warns and stores, fit freezes to EXACT
        # zero, and the model predicts — no ctor injection anywhere
        datav = [_mf_datum(e; M = 1.2 .* e, mode = 4) for e in ev]
        dsv = @test_logs (:warn, r"vanish identically") (:warn,
                          r"structurally dependent") match_mode = :any MomentDataset(
            mbv, datav; gate_eps = 1e-8)
        @test dsv.vanishing == mrv.vanishing
        fv = @test_logs (:warn, r"rank deficient") (:warn, r"rank deficient") fit(
            MomentFit, dsv)
        @test all(fv.coeffs[mrv.vanishing] .== 0.0)
        @test all(fv.coeffs_ungated[mrv.vanishing] .== 0.0)
        @test all(isfinite, residuals(fv))
        modelv = MomentModel(fv)
        @test all(isfinite, predict_moment(modelv, ev[1]))
        # the GAR reduction rides the same freeze on a REAL vanishing set
        fg = fit(MomentFit, dsv, GroupAdaptiveRidge(mbv; lambda = 1e-6))
        @test all(fg.coeffs[mrv.vanishing] .== 0.0)
    end

    @testset "salc_groups: same mark atoms, different mark sites split" begin
        # Review-found merge case (upstream): a canonical member carrying two
        # periodic images of one atom projects two distinct mark placements onto
        # the SAME reference-cell atom set — an atom-set-only key folds them. On
        # this cell the two columns also alias on periodic data (|cos| = 1; the
        # resolvability layer discloses that dependency separately), but the
        # group key is STRUCTURAL: mark class, never column values.
        crm = Crystal(Lattice(Matrix(Diagonal([3.0, 4.0, 5.0]))),
                      [0.0 0.5; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
        sgm = _assemble_spacegroup(crm, [SMatrix{3,3,Float64}(Matrix(1.0 * I(3)))],
                                   [SVector{3,Float64}(0, 0, 0)], "P1-img", 1;
                                   tol = 1e-5)
        spm = MomentSpec(; lmax_env = [1], sampled = [true], lmax_mark = 1,
                         nbody = 3, cutoff_pair = 3.3, cutoff_star = 3.3)
        mbm = MomentBasis(crm, spm; backend = _MBFixedSG(sgm))
        salm = SCEFitting.salcs(mbm)
        glm = SCEFitting.salc_groups(mbm)
        # NOTE: `_marksets` re-derives the key the same way `salc_groups` does,
        # so this is a SELF-CONSISTENCY check of the split (the key must see the
        # site set), not an independent oracle — on this cell the two columns
        # alias on periodic data (|cos| = 1), so no design-side row-support
        # oracle exists; the independent oracle is the FeGe testset above
        function _marksets(s)
            mem = s.members[1]
            at = Int[]; st = Int[]
            for t in mem.terms
                ms = findfirst(sl -> sl.factor.channel == SCEFitting.DISP, t.slots)
                push!(st, t.slots[ms].site)
                push!(at, mem.atoms[t.slots[ms].site])
            end
            return (sort!(unique!(at)), sort!(unique!(st)))
        end
        nmm = n_salcs(mbm)
        nfound = 0
        for j = 1:nmm, k = (j + 1):nmm
            (salm[j].key.body, salm[j].key.orbit_id, salm[j].key.decors) ==
                (salm[k].key.body, salm[k].key.orbit_id, salm[k].key.decors) ||
                continue
            aj, sj = _marksets(salm[j])
            ak, sk = _marksets(salm[k])
            if aj == ak && sj != sk
                nfound += 1
                @test glm[j] != glm[k]      # the pre-fix key merged these
            end
        end
        @test nfound >= 2                   # the fixture actually has the case
    end

    @testset "penalty metric: the μ₀ intercepts are exempt, not shrunk" begin
        # `_intercept_columns` against an INDEPENDENT read of the design: a μ₀ column
        # is the indicator of its mark class, so its entries are exactly 0 or 1 and it
        # is constant over every configuration. Nothing else can be.
        data = _mf_data(24)
        ds = _mf_ds(mb, data; gate_eps = 1e-8)
        ic = SCEFitting._intercept_columns(mb)
        @test !isempty(ic)
        for j = 1:p
            col = ds.X[:, j]
            isind = all(v -> v == 0.0 || v == 1.0, col) && any(isone, col)
            @test (j in ic) == isind
        end
        # every μ₀ column belongs to a different mark class, and they partition the
        # rows of the marked atoms
        @test sum(ds.X[:, ic]; dims = 2) == ones(size(ds.X, 1), 1)

        m = penalty_metric(mb; nconfig = 512, seed = 5)
        @test length(m) == p
        @test all(iszero, m[ic])                       # exempt, exactly
        rest = setdiff(1:p, vcat(ic, ds.vanishing))
        @test all(>(0), m[rest])
        # `false` is the deliberate comparison: the intercepts are then penalized
        mp = penalty_metric(mb; free_intercepts = false, nconfig = 512, seed = 5)
        @test all(>(0), mp[ic])
        @test mp[rest] == m[rest]
        # a μ₀ column is the constant 1, so its reference second moment is exactly 1
        @test all(x -> x ≈ 1.0, mp[ic])
    end

    @testset "penalty metric: μ₀ at λ → ∞ is the intercept-only least squares" begin
        # With the metric, the penalty never touches the intercepts, so as λ → ∞ the
        # penalized columns are crushed and μ₀ converges to the least-squares fit of
        # y on the intercept columns ALONE — which, those being class indicators, is
        # each class's mean target. Without the metric μ₀ is shrunk to zero instead.
        # The λ → ∞ limit is the oracle: at finite λ the intercepts still move,
        # `β_F = (X_F'X_F)⁻¹X_F'(y − X_P β_P(λ))`.
        data = _mf_data(24)
        ds = _mf_ds(mb, data; gate_eps = 1e-8)
        ic = SCEFitting._intercept_columns(mb)
        rows = ds.keep
        means = [sum(ds.y[rows] .* ds.X[rows, j]) / sum(ds.X[rows, j]) for j in ic]
        for est in (Ridge(mb; lambda = 1e12, metric_nconfig = 256),
                    AdaptiveRidge(mb; lambda = 1e12, metric_nconfig = 256),
                    GroupAdaptiveRidge(mb; lambda = 1e12, metric_nconfig = 256))
            f = fit(MomentFit, ds, est)
            @test coef(f)[ic] ≈ means rtol = 1e-6
            @test maximum(abs, coef(f)[setdiff(1:p, ic)]) < 1e-6 * maximum(abs, means)
        end
        # the control: no metric ⇒ the same λ crushes μ₀ too
        fu = fit(MomentFit, ds, Ridge(mb; lambda = 1e12, metric = nothing))
        @test maximum(abs, coef(fu)[ic]) < 1e-6 * maximum(abs, means)
    end

    @testset "penalty metric: the freeze reduction covers every estimator" begin
        # A metric makes Ridge/AdaptiveRidge column-structured too, so the
        # vanishing-column freeze has to cut it down with the design.
        m = collect(1.0:p)
        active = trues(p)
        active[[2, 5]] .= false
        for est in (Ridge(; lambda = 0.5, metric = m),
                    AdaptiveRidge(; lambda = 0.5, metric = m),
                    GroupAdaptiveRidge(salc_groups(mb), ones(maximum(salc_groups(mb)));
                                       lambda = 0.5, metric = m))
            red = SCEFitting._reduce_to_active(est, active)
            @test red.metric == m[active]
            @test SCEFitting._reduce_to_active(est, trues(p)) === est
        end
        # length mismatches are refused by name, not by a deep DimensionMismatch
        @test_throws DimensionMismatch SCEFitting._reduce_to_active(
            Ridge(; lambda = 0.5, metric = ones(p + 1)), active)
    end

    @testset "_reduce_to_active: non-uniform weights follow the relabeling" begin
        # the unit-weight convenience ctor cannot see a weight-permutation bug
        gar = GroupAdaptiveRidge([3, 1, 2, 3, 2, 1], [10.0, 20.0, 30.0];
                                 lambda = 0.5, epsilon = 1e-7, max_iter = 7,
                                 tol = 1e-5)
        active = BitVector([true, false, true, true, false, true])
        red = SCEFitting._reduce_to_active(gar, active)
        @test red.column_groups == [1, 2, 1, 3]     # 3,2,3,1 → first-appearance
        @test red.group_weights == [30.0, 20.0, 10.0]
        @test red.group_sizes == [2, 1, 1]
        @test red.lambda == 0.5 && red.epsilon == 1e-7 &&
              red.max_iter == 7 && red.tol == 1e-5
        # all-active passes the estimator through untouched (identity)
        @test SCEFitting._reduce_to_active(gar, trues(6)) === gar
        # an emptied group is relabeled away, weights follow
        red2 = SCEFitting._reduce_to_active(gar, BitVector([true, false, true, true,
                                                            false, false]))
        @test red2.column_groups == [1, 2, 1]
        @test red2.group_weights == [30.0, 20.0]
        # an unstructured estimator is the identity
        @test SCEFitting._reduce_to_active(OLS(), active) === OLS()
    end

    @testset "local field + coverage + simple floor" begin
        # 2-atom Fe chain-like cell with an EXACT ±x tie (frac x = 0, 0.5 on a
        # 3 Å axis): each atom sees the other through TWO tied images, so the
        # pair-consistent field is h₁ = 2ê_other — hand-derivable throughout.
        crf = Crystal(Lattice(Matrix(Diagonal([3.0, 4.0, 5.0]))),
                      [0.0 0.5; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
        sgf = _assemble_spacegroup(crf, [SMatrix{3,3,Float64}(Matrix(1.0 * I(3)))],
                                   [SVector{3,Float64}(0, 0, 0)], "P1-lf", 1;
                                   tol = 1e-5)
        spf = MomentSpec(; lmax_env = [1], sampled = [true], lmax_mark = 1,
                         nbody = 2, cutoff_pair = 1.6)
        mbf = MomentBasis(crf, spf; backend = _MBFixedSG(sgf))

        # hand oracle: ê₁ = ẑ, ê₂ = (sinθ, 0, cosθ) ⇒ h₁(1) = 2ê₂ and
        # h₁(2) = 2ê₁; both norms 2 (tie multiplicity), both alignments cosθ
        th = 0.7
        ef = zeros(3, 2)
        ef[:, 1] = [0.0, 0.0, 1.0]
        ef[:, 2] = [sin(th), 0.0, cos(th)]
        lf = moment_local_field(mbf, [ef])
        @test lf.h1 ≈ [2.0, 2.0]
        @test lf.edoth ≈ [cos(th), cos(th)]
        @test lf.row_config == [1, 1] && lf.row_atom == [1, 2]
        # the SpinDatum method resolves mode 4 to the identity — bitwise
        dm4 = _mf_datum(ef; M = 1.3 .* ef, mode = 4)
        lfd = moment_local_field(mbf, [dm4])
        @test lfd.h1 == lf.h1 && lfd.edoth == lf.edoth
        # mode-1 zero axis: h₁ still computed, alignment undefined
        ax1 = copy(ef)
        ax1[:, 1] .= 0.0
        lf1 = moment_local_field(mbf, [_mf_datum(ef; M = 1.3 .* ef, mode = 1,
                                                 axes = ax1)])
        @test lf1.h1 == lf.h1
        @test isnan(lf1.edoth[1]) && lf1.edoth[2] ≈ cos(th)
        # mode 1 with a genuinely different axis reads THAT axis: ê₁ = x̂ against
        # h₁(1) = 2ê₂ gives sinθ
        axx = copy(ef)
        axx[:, 1] .= [1.0, 0.0, 0.0]
        lfx = moment_local_field(mbf, [_mf_datum(ef; M = 1.3 .* ef, mode = 1,
                                                 axes = axx)])
        @test lfx.edoth[1] ≈ sin(th) && lfx.edoth[2] ≈ cos(th)
        # a datum without a mode is refused (the axis rule needs the class)
        nomode = SpinDatum(0.0, ef, ones(2), zeros(3, 2), zeros(3, 2),
                           1.3 .* ef, nothing, nothing)
        @test_throws ArgumentError moment_local_field(mbf, [nomode])
        # doors: non-unit config column; non-unit nonzero marked axis
        bad = copy(ef)
        bad[:, 2] .*= 1.5
        @test_throws ArgumentError moment_local_field(mbf, [bad])
        @test_throws ArgumentError moment_local_field(mbf, [ef]; axes = [bad])
        @test_throws ArgumentError moment_local_field(mbf, [ef]; axes = [ef, ef])
        # a species the basis never reads (lmax_env = 0) contributes nothing:
        # Fe–Ge cell, Ge environment off ⇒ the Fe row sees h₁ = 0 (edoth NaN),
        # the Ge row still reads its Fe neighbors
        crg = Crystal(Lattice(Matrix(Diagonal([3.0, 4.0, 5.0]))),
                      [0.0 0.5; 0.0 0.0; 0.0 0.0], [1, 2], ["Fe", "Ge"])
        sgg = _assemble_spacegroup(crg, [SMatrix{3,3,Float64}(Matrix(1.0 * I(3)))],
                                   [SVector{3,Float64}(0, 0, 0)], "P1-lg", 1;
                                   tol = 1e-5)
        spg = MomentSpec(; lmax_env = [1, 0], sampled = [true, false],
                         lmax_mark = 1, nbody = 2, cutoff_pair = 1.6,
                         marked = [true, true])
        mbg2 = MomentBasis(crg, spg; backend = _MBFixedSG(sgg))
        lfg = moment_local_field(mbg2, [ef])
        @test lfg.h1[1] == 0.0 && isnan(lfg.edoth[1])       # Fe: Ge env is off
        @test lfg.h1[2] ≈ 2.0 && lfg.edoth[2] ≈ cos(th)     # Ge: reads Fe

        # atoms ≠ slots: the Ge-only FeGe basis (marked 5:8) must report ATOM
        # numbers and read the ATOM's neighbors — a slot-indexed implementation
        # would read atoms 1:4's environment
        specg8 = MomentSpec(; lmax_env = [1, 1], sampled = [true, true],
                            lmax_mark = 1, nbody = 2, cutoff_pair = 3.0,
                            marked = [false, true])
        mbg8 = MomentBasis(xt, specg8; backend = bk)
        specall = MomentSpec(; lmax_env = [1, 1], sampled = [true, true],
                             lmax_mark = 1, nbody = 2, cutoff_pair = 3.0,
                             marked = [true, true])
        mball = MomentBasis(xt, specall; backend = bk)
        e8 = _mb_unit(rng, nat)
        lf8 = moment_local_field(mbg8, [e8])
        lfall = moment_local_field(mball, [e8])
        @test lf8.row_atom == [5, 6, 7, 8]
        # a row's field depends on the atom, not on which atoms are marked
        @test lf8.h1 == lfall.h1[5:8] && lf8.edoth == lfall.edoth[5:8]
        # perturbing Fe atom 1 (an environment atom of every Ge within 3.0 Å)
        # moves every Ge row's field; re-axing atom 7 moves only edoth[3]
        e8b = copy(e8); e8b[:, 1] .= _mb_unit(rng, 1)[:, 1]
        lf8b = moment_local_field(mbg8, [e8b])
        @test all(lf8b.h1 .!= lf8.h1)
        ax8 = copy(e8); ax8[:, 7] .= _mb_unit(rng, 1)[:, 1]
        lf8c = moment_local_field(mbg8, [e8]; axes = [ax8])
        @test lf8c.h1 == lf8.h1
        @test lf8c.edoth[[1, 2, 4]] == lf8.edoth[[1, 2, 4]] && lf8c.edoth[3] != lf8.edoth[3]

        # a 1-body basis reads no environment: the field is empty by construction
        # (cutoff_pair is a required spec field the model never uses there)
        sp1 = MomentSpec(; lmax_env = [1], sampled = [true], lmax_mark = 1,
                         nbody = 1, cutoff_pair = 1.6)
        mb1 = MomentBasis(crf, sp1; backend = _MBFixedSG(sgf))
        @test SCEFitting._pair_neighbors(mb1) == Dict(1 => Int[], 2 => Int[])
        lf1b = moment_local_field(mb1, [ef])
        @test lf1b.h1 == [0.0, 0.0] && all(isnan, lf1b.edoth)

        # the neighbor tie band is the BASIS's: a near tie (Δd/d ≈ 4e-5) is one
        # image at the default band and two at tie_tol = 1e-4, and h₁ follows
        crn = Crystal(Lattice(Matrix(Diagonal([3.0, 4.0, 5.0]))),
                      [0.0 0.5 + 1e-5; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
        sgn = _assemble_spacegroup(crn, [SMatrix{3,3,Float64}(Matrix(1.0 * I(3)))],
                                   [SVector{3,Float64}(0, 0, 0)], "P1-near", 1;
                                   tol = 1e-5)
        mbn_default = MomentBasis(crn, spf; backend = _MBFixedSG(sgn))
        mbn_wide = MomentBasis(crn, spf; backend = _MBFixedSG(sgn), tie_tol = 1e-4)
        @test SCEFitting._pair_neighbors(mbn_default)[1] == [2]
        @test SCEFitting._pair_neighbors(mbn_wide)[1] == [2, 2]
        @test moment_local_field(mbn_default, [ef]).h1 ≈ [1.0, 1.0]
        @test moment_local_field(mbn_wide, [ef]).h1 ≈ [2.0, 2.0]

        # coverage: hand-checkable quantile + fractions
        tr = (; h1 = collect(0.1:0.1:10.0), edoth = fill(0.9, 100))
        nw = (; h1 = [0.5, 9.98, 11.0, 12.0], edoth = [0.8, -0.2, NaN, -0.3])
        cov = moment_coverage(tr, nw; q = 0.99)
        @test cov.threshold ≈ quantile(tr.h1, 0.99)
        @test cov.threshold ≈ 9.901                # hand: 0.1·(1 + 0.99·99)
        @test cov.frac_beyond ≈ 3 / 4          # 9.98, 11, 12 exceed
        @test cov.frac_anti ≈ 2 / 3            # NaN row excluded from the denominator
        @test cov.n_rows == 4 && cov.n_defined == 3
        @test_throws ArgumentError moment_coverage(tr, nw; q = 1.0)
        @test_throws ArgumentError moment_coverage((; h1 = Float64[], edoth = Float64[]),
                                                   nw)
        @test isnan(moment_coverage(tr, (; h1 = [1.0], edoth = [NaN])).frac_anti)

        # ---- the nested simple-feature floor (M2-5) ----
        # plant y EXACTLY in the feature space: y = a₀ + b₀·Σ_j ê_i·ê_j, moments
        # decomposable by construction (M ∥ ê ⇒ gate 0). On this pair basis the
        # P₁ shell sum lies in the SALC span, so the SALC fit is nested above
        # the floor. The tied cell carries 2 structural dependencies (warned)
        # and the OLS solves warn rank-deficient — both expected and asserted.
        nbf = SCEFitting._pair_neighbors(mbf)
        @test nbf[1] == [2, 2] && nbf[2] == [1, 1]     # the ±x tie, by hand
        a0, b0 = 1.1, 0.07
        dataf = SpinDatum[]
        for c = 1:12
            ec = _mb_unit(rng, 2)
            M = zeros(3, 2)
            for a = 1:2
                f1 = sum(ec[:, a]' * ec[:, j] for j in nbf[a])
                M[:, a] = (a0 + b0 * f1) .* ec[:, a]
            end
            push!(dataf, _mf_datum(ec; M = M, mode = 4))
        end
        dsf = @test_logs (:warn, r"structurally dependent") match_mode = :any MomentDataset(
            mbf, dataf; gate_eps = 1e-8)
        @test all(dsf.keep)
        ff = @test_logs (:warn, r"rank deficient") (:warn, r"rank deficient") fit(
            MomentFit, dsf)
        fl = moment_simple_floor(ff, dataf; lmax = 1)
        @test fl.sigma_floor < 1e-12            # the floor model IS the truth
        @test fl.sigma_model <= fl.sigma_floor + 1e-12   # nested (P₁ ∈ span)
        @test all(<(1e-10), fl.inclusion)       # every feature representable
        @test fl.n_features == 4                # 2 orbits × (intercept + P₁)
        # recover the planted parameters from the floor coefficients
        @test all(isapprox.(fl.coef[1:2:end], a0; atol = 1e-10))
        @test all(isapprox.(fl.coef[2:2:end], b0; atol = 1e-10))
        # lmax = 2 on an lmax_env = 1 basis: the P₂ column is NOT representable
        # — the inclusion report says so instead of silently claiming a bound
        fl2 = moment_simple_floor(ff, dataf; lmax = 2)
        p2cols = [k for (k, lb) in enumerate(fl2.feature_labels)
                  if occursin("P2", lb)]
        # (3x² − 1 against span{1, x} with x = ê₁·ê₂ uniform on [−1, 1]: the
        # relative residual is ≈ 1 in expectation; 0.1 leaves room for a 12-sample
        # draw while staying 9 decades above the representable columns)
        @test all(>(0.1), fl2.inclusion[p2cols])
        @test all(<(1e-10), fl2.inclusion[setdiff(1:fl2.n_features, p2cols)])
        # the bound flag and the disclosure counts
        @test fl.nested_bound === true          # OLS fit: the bound applies
        @test fl.n_rows == count(dsf.keep) && 1 <= fl.design_rank <= 6
        fgar = fit(MomentFit, dsf, GroupAdaptiveRidge(mbf; lambda = 1e-3))
        flg = moment_simple_floor(fgar, dataf; lmax = 1)
        @test flg.nested_bound === false        # shrinkage: no nested bound
        # MODE 1: the floor's feature is P_l(ê_i·e_j) with ê the ROW AXIS, not
        # the direction — plant it with random axes (M ∥ ê ⇒ gate 0) and recover
        a1, b1 = -0.4, 0.13
        data1 = SpinDatum[]
        for c = 1:12
            ec = _mb_unit(rng, 2)
            axc = _mb_unit(rng, 2)
            M = zeros(3, 2)
            for a = 1:2
                f1 = sum(axc[:, a]' * ec[:, j] for j in nbf[a])
                M[:, a] = (a1 + b1 * f1) .* axc[:, a]
            end
            push!(data1, _mf_datum(ec; M = M, mode = 1, axes = axc))
        end
        ds1 = @test_logs (:warn, r"structurally dependent") match_mode = :any MomentDataset(
            mbf, data1; gate_eps = 1e-8)
        @test all(ds1.keep)
        ff1 = @test_logs (:warn, r"rank deficient") (:warn, r"rank deficient") fit(
            MomentFit, ds1)
        fl1 = moment_simple_floor(ff1, data1; lmax = 1)
        @test fl1.sigma_floor < 1e-12
        @test all(isapprox.(fl1.coef[1:2:end], a1; atol = 1e-10))
        @test all(isapprox.(fl1.coef[2:2:end], b1; atol = 1e-10))
        # the replay config is the first one whose MARKED AXES are all nonzero:
        # zero atom 2's axis in configs 1–4 → config 5 replays, no refusal, and
        # the floor runs on the kept rows only
        data1z = copy(data1)
        for c = 1:4
            axz = copy(data1[c].constraint_axes); axz[:, 2] .= 0.0
            data1z[c] = _mf_datum(data1[c].directions; M = data1[c].moments_bare,
                                  mode = 1, axes = axz)
        end
        ds1z = @test_logs (:warn, r"structurally dependent") (:info,
            r"undefined-axis") match_mode = :any MomentDataset(mbf, data1z;
                                                               gate_eps = 1e-8)
        @test count(ds1z.keep) == 24 - 4
        ff1z = @test_logs (:warn, r"rank deficient") (:warn, r"rank deficient") fit(
            MomentFit, ds1z)
        fl1z = moment_simple_floor(ff1z, data1z; lmax = 1)
        @test fl1z.n_rows == 20 && fl1z.sigma_floor < 1e-12
        # every config carries a zero axis → the replay door falls back to the
        # target-only door (disclosed in the docstring) and still returns
        data1a = [_mf_datum(d.directions; M = d.moments_bare, mode = 1,
                            axes = (ax = copy(d.constraint_axes); ax[:, 2] .= 0.0; ax))
                  for d in data1]
        ds1a = @test_logs (:warn, r"structurally dependent") (:info,
            r"undefined-axis") match_mode = :any MomentDataset(mbf, data1a;
                                                               gate_eps = 1e-8,
                                                               coverage_floor = 0.0)
        ff1a = @test_logs (:warn, r"rank deficient") (:warn, r"rank deficient") fit(
            MomentFit, ds1a)
        fl1a = moment_simple_floor(ff1a, data1a; lmax = 1)
        @test fl1a.n_rows == 12 && fl1a.n_features == 2     # atom-1 orbit only
        # the floor reads the GATED rows and the GATED residuals: on a dataset
        # with gate-rejected rows, n_rows = kept count and sigma_model is the
        # gated (exact) fit, not the polluted ungated one
        datag = SpinDatum[]
        for c = 1:12
            ec = _mb_unit(rng, 2)
            M = zeros(3, 2)
            for a = 1:2
                f1 = sum(ec[:, a]' * ec[:, j] for j in nbf[a])
                M[:, a] = (a0 + b0 * f1) .* ec[:, a]
            end
            if c <= 3                          # dirty: transverse + corrupted
                t = [0.0, 1.0, 0.0]; t .-= dot(t, ec[:, 1]) .* ec[:, 1]
                M[:, 1] .+= 0.3 .* t ./ norm(t) .+ 0.5 .* ec[:, 1]
            end
            push!(datag, _mf_datum(ec; M = M, mode = 4))
        end
        dsg = @test_logs (:warn, r"structurally dependent") (:info,
            r"excluded") match_mode = :any MomentDataset(mbf, datag; gate_eps = 1e-6)
        @test count(dsg.keep) == 24 - 3
        ffg = @test_logs (:warn, r"rank deficient") (:warn, r"rank deficient") fit(
            MomentFit, dsg)
        flgd = moment_simple_floor(ffg, datag; lmax = 1)
        @test flgd.n_rows == 21
        @test flgd.sigma_model < 1e-8 && flgd.sigma_floor < 1e-12
        @test std(residuals(ffg; gated = false)) > 1e-2
        # _design_moment's shape door is a plain ArgumentError (not a
        # TaskFailedException from inside the threaded loop)
        @test_throws ArgumentError _design_moment(mbf, [randn(3, 4)], [randn(3, 4)])
        @test_throws ArgumentError _design_moment(mbf, [ef], [randn(3, 3)])
        # the basis records its tie band; diagnostics read it back
        @test mbf.tie_tol == SCEFitting._SAME_DIST_RTOL
        @test_throws ArgumentError moment_local_field(mbf, Matrix{Float64}[])
        # pairing doors: wrong length; silently re-paired data (targets move)
        @test_throws ArgumentError moment_simple_floor(ff, dataf[1:5]; lmax = 1)
        swapped = copy(dataf)
        swapped[1], swapped[2] = swapped[2], swapped[1]
        @test_throws ArgumentError moment_simple_floor(ff, swapped; lmax = 1)
        @test_throws ArgumentError moment_simple_floor(ff, dataf; lmax = 4)
        # environment door: same TARGETS, moved environment — the target check
        # alone cannot see it (y = ê·M is invariant under co-rotating M with a
        # new ê), the design-row replay must
        moved = copy(dataf)
        eold = dataf[1].directions
        yold = [dot(eold[:, a], dataf[1].moments_bare[:, a]) for a = 1:2]
        enew = _mb_unit(rng, 2)
        Mnew = hcat((yold[a] .* enew[:, a] for a = 1:2)...)
        moved[1] = _mf_datum(enew; M = Mnew, mode = 4)
        @test_throws ArgumentError moment_simple_floor(ff, moved; lmax = 1)
        # Legendre kernel by hand at x = 0.3
        @test SCEFitting._legendre(1, 0.3) == 0.3
        @test SCEFitting._legendre(2, 0.3) ≈ (3 * 0.09 - 1) / 2
        @test SCEFitting._legendre(3, 0.3) ≈ (5 * 0.027 - 0.9) / 2
        @test_throws ArgumentError SCEFitting._legendre(4, 0.3)
    end

    @testset "persistence is refused by name (design record §4.2)" begin
        f = _mf_fit(_mf_ds(mb, _mf_data(6); gate_eps = 1e-8))
        path = tempname() * ".toml"
        for x in (mb, f, MomentModel(f))
            err = try; SCEFitting.save(path, x); nothing; catch e; e; end
            @test err isa ArgumentError && occursin("not persisted", err.msg)
            @test_throws ArgumentError SCEFitting.save(IOBuffer(), x)
        end
        @test !isfile(path)
    end
end
