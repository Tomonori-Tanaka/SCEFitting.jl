using Test
using SCEFitting
using SCEFitting: solve_coefficients, salc_groups, group_costs, cost_weights
using LinearAlgebra
using Statistics
using Random
import Tables

# A centered synthetic regression fixture: X column-centered, y centered, so the
# estimator contract (`solve_coefficients` receives a centered problem) holds.
function _centered_problem(rng, n, btrue; noise = 0.01)
    p = length(btrue)
    X = randn(rng, n, p)
    X .-= mean(X; dims = 1)
    y = X * btrue .+ noise .* randn(rng, n)
    y .-= mean(y)
    return X, y
end

@testset "cost-weighted group selection" begin
    @testset "GroupAdaptiveRidge construction and validation" begin
        est = GroupAdaptiveRidge([1, 1, 2], [1.0, 2.0]; lambda = 0.5)
        @test est.lambda === 0.5
        @test est.column_groups == [1, 1, 2]
        @test est.group_weights == [1.0, 2.0]
        @test est.group_sizes == [2, 1]
        @test est.epsilon === 1e-8 && est.max_iter === 50 && est.tol === 1e-6
        @test SCEFitting.islinear(est)

        # invalid scalars
        @test_throws ArgumentError GroupAdaptiveRidge([1], [1.0]; lambda = -1.0)
        @test_throws ArgumentError GroupAdaptiveRidge([1], [1.0]; lambda = Inf)
        @test_throws ArgumentError GroupAdaptiveRidge([1], [1.0];
                                                      lambda = 1.0, epsilon = 0.0)
        @test_throws ArgumentError GroupAdaptiveRidge([1], [1.0];
                                                      lambda = 1.0, max_iter = 0)
        @test_throws ArgumentError GroupAdaptiveRidge([1], [1.0]; lambda = 1.0, tol = 0.0)
        # invalid labels / weights
        @test_throws ArgumentError GroupAdaptiveRidge(Int[], Float64[]; lambda = 1.0)
        @test_throws ArgumentError GroupAdaptiveRidge([1, 3, 3], [1.0, 1.0, 1.0];
                                                      lambda = 1.0)
        @test_throws ArgumentError GroupAdaptiveRidge([0, 1], [1.0]; lambda = 1.0)
        @test_throws ArgumentError GroupAdaptiveRidge([1, 2], [1.0]; lambda = 1.0)
        @test_throws ArgumentError GroupAdaptiveRidge([1, 2], [1.0, -1.0]; lambda = 1.0)
        @test_throws ArgumentError GroupAdaptiveRidge([1, 2], [1.0, 0.0]; lambda = 1.0)
        @test_throws ArgumentError GroupAdaptiveRidge([1, 2], [1.0, NaN]; lambda = 1.0)

        # construction copies: caller mutation cannot leak in
        cg = [1, 1, 2]
        wv = [1.0, 2.0]
        est2 = GroupAdaptiveRidge(cg, wv; lambda = 1.0)
        cg[1] = 99
        wv[1] = -5.0
        @test est2.column_groups == [1, 1, 2]
        @test est2.group_weights == [1.0, 2.0]

        # label length must match the design-matrix column count
        rng = MersenneTwister(1)
        X, y = _centered_problem(rng, 20, zeros(4))
        @test_throws DimensionMismatch solve_coefficients(est, X, y)   # 3 labels, 4 cols

        # compact show (never dumps the label vectors)
        s = sprint(show, est)
        @test occursin("GroupAdaptiveRidge(3 columns in 2 groups", s)
        @test occursin("lambda=0.5", s)
    end

    @testset "singleton groups + unit weights reproduce AdaptiveRidge" begin
        rng = MersenneTwister(42)
        btrue = zeros(10)
        btrue[2] = 1.5
        btrue[7] = -2.0
        X, y = _centered_problem(rng, 80, btrue)
        lam = 1e-3
        b_ar = solve_coefficients(AdaptiveRidge(; lambda = lam), X, y)
        b_gar = solve_coefficients(
            GroupAdaptiveRidge(collect(1:10), ones(10); lambda = lam), X, y)
        # mathematically identical iterations (wⱼ = 1/(βⱼ² + ε), same init, same
        # convergence test) through separate accumulation code — pin tightly, not bitwise
        @test isapprox(b_gar, b_ar; rtol = 1e-10)
        # exact lambda == 0 routes to the same QR OLS path
        est0 = GroupAdaptiveRidge(collect(1:10), ones(10); lambda = 0.0)
        @test solve_coefficients(est0, X, y) == solve_coefficients(OLS(), X, y)
    end

    @testset "group-sparse recovery and weight monotonicity" begin
        rng = MersenneTwister(7)
        labels = repeat(1:4; inner = 3)          # p = 12, four groups of three
        btrue = zeros(12)
        btrue[4:6] .= [1.0, -0.8, 0.6]           # signal on group 2 only
        X, y = _centered_problem(rng, 120, btrue)
        b = solve_coefficients(GroupAdaptiveRidge(labels, ones(4); lambda = 1e-2), X, y)
        active = b[4:6]
        inactive = b[[1:3; 7:12]]
        @test minimum(abs, active) > 50 * maximum(abs, inactive)   # off-groups crushed
        @test isapprox(active, btrue[4:6]; rtol = 0.05)

        # a heavier fixed weight on a group strictly shrinks that group's norm
        btrue2 = copy(btrue)
        btrue2[7:9] .= [0.5, 0.4, -0.3]          # add signal on group 3
        y2 = X * btrue2 .+ 0.01 .* randn(rng, 120)
        y2 .-= mean(y2)
        w_hi = [1.0, 1.0, 2.0, 1.0]
        b_lo = solve_coefficients(GroupAdaptiveRidge(labels, ones(4); lambda = 0.05), X, y2)
        b_hi = solve_coefficients(GroupAdaptiveRidge(labels, w_hi; lambda = 0.05), X, y2)
        @test norm(b_hi[7:9]) < norm(b_lo[7:9])
    end

    # --- basis-driven helpers -------------------------------------------------------
    lat = Lattice(Matrix(3.0 * I(3)))
    crystal = Crystal(lat, [0.2 -0.2; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
    # NoSymmetry (P1): anisotropic basis with several (orbit, ls) groups split into
    # Lf/block channels, plus a small isotropic basis for by-hand entry counting
    basis = SCEBasis(crystal, BasisSpec(; nbody = 2, cutoff = 1.5, lmax = [2],
                                        isotropy = false))
    small = SCEBasis(crystal, BasisSpec(; nbody = 2, cutoff = 1.5, lmax = [1],
                                        isotropy = true))

    @testset "salc_groups: contiguous labels matching (body, orbit_id, ls) runs" begin
        ks = basis.salc_basis.keys
        labels = salc_groups(basis)
        m = length(ks)
        @test length(labels) == m == n_salcs(basis)
        G = maximum(labels)
        @test sort(unique(labels)) == 1:G
        tup = [(k.body, k.orbit_id, k.decors) for k in ks]
        @test G == length(unique(tup))
        @test all((labels[i] == labels[j]) == (tup[i] == tup[j])
                  for i = 1:m for j = 1:m)
        # at least one group carries multiple Lf/block channels (else the fixture is
        # too weak to distinguish group from column selection)
        @test G < m
    end

    @testset "group_costs: brute-force union, additivity, validation" begin
        # both fixtures: `small` (single group) and `basis` (many (orbit, ls) groups,
        # so the additivity claim is exercised with G > 1)
        for b in (small, basis)
            labels = salc_groups(b)
            costs = group_costs(b, labels)
            G = maximum(labels)
            @test length(costs) == G

            # independent brute force with a structurally different key encoding:
            # nested tuples rather than the vectors `_EntryKey` uses, and the slot
            # labels spelled out here rather than routed through `_slotkey` or
            # `_term_spin_ls`.
            bysets = [Set{Any}() for _ = 1:G]
            for (j, s) in enumerate(salcs(b))
                for mem in s.members, t in mem.terms
                    slotkey = Tuple((Int(sl.factor.channel), sl.site,
                                     sl.factor.k, sl.factor.l) for sl in t.slots)
                    for idx in CartesianIndices(t.folded)
                        t.folded[idx] == 0.0 && continue
                        key = (Tuple(mem.atoms), Tuple(Tuple.(mem.shifts)),
                               slotkey, Tuple(idx))
                        push!(bysets[labels[j]], key)
                    end
                end
            end
            # The identity the slot price rests on, asserted directly: on a
            # pure-spin basis every term's slot list is the identity site map, so
            # the slot count IS the body order. A regression that made a term
            # carry two slots on one site would leave every stored number
            # unchanged and only show up here.
            @test all(length(t.slots) == s.body == length(mem.atoms) &&
                      all(sl.factor.channel === SCEFitting.SPIN for sl in t.slots) &&
                      [sl.site for sl in t.slots] == collect(eachindex(mem.atoms))
                      for s in salcs(b) for mem in s.members for t in mem.terms)

            # sweep pricing: each distinct entry costs one site-program slot per
            # SLOT of the term (k[3] of the brute-force key is the slot tuple).
            # On a pure-spin basis every slot is a distinct site, so this is still
            # the member's site count = body order — the equality below is what
            # pins that, and it is NOT caught by byte-equality of the basis.
            # [Backported from SLCE.jl a596ea3.]
            slotsum(S) = sum(k -> length(k[3]), S; init = 0)
            @test costs == slotsum.(bysets)
            # the sweep price is bounded below by the bare entry count (a 1-body
            # entry costs 1 slot, an N-body one N slots)
            @test all(costs .>= length.(bysets))

            # additivity: groups never share entries, so the sum is the global
            # count — equal to the single-group partition's cost
            total = group_costs(b, ones(Int, n_salcs(b)))
            @test length(total) == 1
            @test sum(costs) == total[1]
            @test total[1] == slotsum(union(bysets...))

            # default labels argument = salc_groups partition
            @test group_costs(b) == costs
        end
        @test maximum(salc_groups(basis)) > 1     # the multi-group case was real

        # malformed labels
        @test_throws ArgumentError group_costs(small, ones(Int, n_salcs(small) + 1))
        badgap = ones(Int, n_salcs(small))
        badgap[1] = 3                          # skips label 2
        @test_throws ArgumentError group_costs(small, badgap)
    end

    @testset "cost_weights and the basis convenience constructor" begin
        labels = salc_groups(basis)
        c = group_costs(basis, labels)
        G = maximum(labels)
        p = [count(==(g), labels) for g = 1:G]

        lw0 = cost_weights(basis; theta = 0.0)
        @test lw0.labels == labels
        @test lw0.weights == sqrt.(p)                       # exactly √p_g at theta = 0

        lw1 = cost_weights(basis; theta = 1.0)
        @test lw1.weights ≈ sqrt.(p) .* (c ./ mean(c))

        lw5 = cost_weights(basis; theta = 0.5)
        @test all(min.(lw0.weights, lw1.weights) .<= lw5.weights .+ 1e-12) &&
              all(lw5.weights .<= max.(lw0.weights, lw1.weights) .+ 1e-12)

        @test_throws ArgumentError cost_weights(basis; theta = -0.1)
        @test_throws ArgumentError cost_weights(basis; theta = 1.1)

        est = GroupAdaptiveRidge(basis; lambda = 1e-3, theta = 1.0)
        @test est.column_groups == labels
        @test est.group_weights == lw1.weights
        @test est.lambda === 1e-3
    end

    # --- GCV / effective degrees of freedom -----------------------------------------
    @testset "effective_dof / gcv match the dense hat matrix (n > p)" begin
        rng = MersenneTwister(23)
        nconf = 30
        configs = [Matrix(randcfg(rng, 2)) for _ = 1:nconf]
        energies = randn(rng, nconf)
        ds = SCEDataset(small, configs, energies)
        X, y, _, _, _ = SCEFitting._assemble_problem(ds, 0.0)
        n = size(X, 1)
        @test n > size(X, 2)                       # overdetermined fixture

        lam = 0.3
        fr = fit(SCEFit, ds, Ridge(; lambda = lam))
        Hr = X * ((X' * X + lam * I) \ X')
        @test effective_dof(fr) ≈ tr(Hr) + 1 rtol = 1e-8
        rss_r = sum(abs2, y .- X * coef(fr))
        @test gcv(fr) ≈ n * rss_r / (n - (tr(Hr) + 1))^2 rtol = 1e-8

        # adaptive members: converged weights recomputed from the fitted coefficients
        fa = fit(SCEFit, ds, AdaptiveRidge(; lambda = 1e-3))
        Da = Diagonal(1.0 ./ (coef(fa) .^ 2 .+ 1e-8))
        Ha = X * ((X' * X + 1e-3 * Da) \ X')
        @test effective_dof(fa) ≈ tr(Ha) + 1 rtol = 1e-8

        glabels = salc_groups(small)
        fg = fit(SCEFit, ds, GroupAdaptiveRidge(glabels, ones(maximum(glabels));
                                                lambda = 1e-3))
        beta = coef(fg)
        wg = [1.0 / (sum(abs2, beta[glabels .== glabels[j]]) +
                     count(==(glabels[j]), glabels) * 1e-8) for j in eachindex(beta)]
        Hg = X * ((X' * X + 1e-3 * Diagonal(wg)) \ X')
        @test effective_dof(fg) ≈ tr(Hg) + 1 rtol = 1e-8
        @test gcv(fg) ≈ n * sum(abs2, y .- X * beta) / (n - (tr(Hg) + 1))^2 rtol = 1e-8

        fo = fit(SCEFit, ds, OLS())
        @test effective_dof(fo) ≈ rank(X) + 1
        @test isfinite(gcv(fo))
    end

    @testset "_edof with unpenalized columns (dense hat reference)" begin
        rng = MersenneTwister(4127)
        # Independent oracle: the dense hat matrix of the penalized normal equations,
        # `tr(X(X'X + λD)⁻¹X')`, formed here without any of `_edof`'s machinery.
        dense_df(X, lam, d) = tr(X * ((X' * X + lam * Diagonal(d)) \ X'))

        lam = 0.37
        # (a) overdetermined: no / one / several unpenalized columns, and the cached
        #     Gram keyword (whose `D^{-1/2}` form is the one that divides by zero)
        n, p = 40, 12
        X = randn(rng, n, p)
        base = 0.3 .+ rand(rng, p)
        for free in (Int[], [3], [1, 7, 12])
            d = copy(base)
            d[free] .= 0.0
            @test SCEFitting._edof(X, lam, d) ≈ dense_df(X, lam, d) rtol = 1e-9
            @test SCEFitting._edof(X, lam, d; XtX = X' * X) ≈
                  dense_df(X, lam, d) rtol = 1e-9
        end

        # (b) underdetermined: the dual side, with the unpenalized block still thin
        n2, p2 = 9, 25
        X2 = randn(rng, n2, p2)
        d2 = 0.3 .+ rand(rng, p2)
        d2[[2, 5]] .= 0.0
        @test SCEFitting._edof(X2, lam, d2) ≈ dense_df(X2, lam, d2) rtol = 1e-9

        # (c) λ → ∞ leaves exactly the unpenalized degrees of freedom, and an
        #     all-unpenalized diagonal is the plain design rank
        d3 = copy(base)
        d3[[2, 4, 9]] .= 0.0
        @test SCEFitting._edof(X, 1e12, d3) ≈ 3.0 rtol = 1e-6
        @test SCEFitting._edof(X, lam, zeros(p)) == Float64(p)

        # (d) a rank-deficient unpenalized block leaves the fit unidentified: refuse
        #     by name rather than report a finite dof for it
        Xd = copy(X)
        Xd[:, 5] = Xd[:, 3]
        dd = copy(base)
        dd[[3, 5]] .= 0.0
        @test_throws ArgumentError SCEFitting._edof(Xd, lam, dd)
    end

    @testset "gcv in the underdetermined regime and guards" begin
        rng = MersenneTwister(31)
        nconf = 12
        configs = [Matrix(randcfg(rng, 2)) for _ = 1:nconf]
        # an in-span group-sparse target: a pure-noise target would drive the adaptive
        # fit to interpolation (df → n) and the GCV guard to Inf by design
        ds0 = SCEDataset(basis, configs, zeros(nconf))
        btrue = zeros(n_salcs(basis))
        btrue[1:3] .= [0.5, -0.3, 0.2]
        energies = 0.7 .+ ds0.X_E * btrue .+ 0.001 .* randn(rng, nconf)
        ds = SCEDataset(basis, configs, energies)   # p = n_salcs(basis) > 12 rows
        @test n_salcs(basis) > nconf
        f = fit(SCEFit, ds, GroupAdaptiveRidge(basis; lambda = 1e-2))
        @test effective_dof(f) - 1 < nconf
        @test isfinite(gcv(f))
        # near-interpolating: the guard returns Inf instead of a meaningless score
        # (column centering makes the rows sum to zero, so rank ≤ n − 1 and df + 1 → n)
        fni = fit(SCEFit, ds, Ridge(; lambda = 1e-14))
        @test gcv(fni) == Inf
        # non-linear estimators are rejected
        fp = fit(SCEFit, ds, PrecomputedPilot(zeros(n_salcs(basis))))
        @test_throws ArgumentError gcv(fp)
        @test_throws ArgumentError effective_dof(fp)
    end

    @testset "gcv on a torque co-fit (smoke; grouped CV is the ground truth)" begin
        rng = MersenneTwister(37)
        nconf = 10
        configs = [Matrix(randcfg(rng, 2)) for _ = 1:nconf]
        energies = randn(rng, nconf)
        torques = [randn(rng, 3, 2) for _ = 1:nconf]
        ds = SCEDataset(small, configs, energies, torques)
        f = fit(SCEFit, ds, Ridge(; lambda = 0.1); torque_weight = 0.3)
        @test isfinite(gcv(f))
        @test effective_dof(f) > 1

        # M4 accounting at torque-only weight: the energy rows enter with weight
        # exactly zero, so GCV must not count them (n_eff = n_T) and must not
        # charge the intercept to the counted rows. Wiring check: gcv == the
        # uncounted-corrected _gcv_score, ≠ the naive one — both inequalities
        # assert the correction is ACTIVE. [Backported from SLCE.jl 54457ca.]
        fw1 = fit(SCEFit, ds, Ridge(; lambda = 0.1); torque_weight = 1.0)
        X1, y1, _, _, _ = SCEFitting._assemble_problem(ds, 1.0)
        lam1, w1 = SCEFitting._penalty_diagonal(fw1.estimator, fw1.jphi)
        neff1 = SCEFitting._gcv_neff(fw1)
        @test neff1 == length(ds.y_T)
        @test gcv(fw1) == first(SCEFitting._gcv_score(X1, y1, fw1.jphi, lam1, w1;
                                                      n_eff = neff1, intercept = 0.0))
        @test gcv(fw1) != first(SCEFitting._gcv_score(X1, y1, fw1.jphi, lam1, w1))
        @test effective_dof(fw1) ==
              SCEFitting._edof(X1, lam1, w1)              # no intercept charged at w = 1
    end

    @testset "diagnostics refuse a refit result; refit refuses a GroupAdaptiveRidge" begin
        # [Backported from SLCE.jl 54457ca, review M3.] effective_dof/gcv rebuild
        # the FULL design; on a refit (solved on a support) that df/score describes
        # a model the refit did not return — refuse by name.
        rng = MersenneTwister(53)
        nconf = 30
        configs = [Matrix(randcfg(rng, 2)) for _ = 1:nconf]
        ds0 = SCEDataset(basis, configs, zeros(nconf))
        btrue = zeros(n_salcs(basis))
        btrue[1:2] .= [0.5, -0.3]
        energies = 0.2 .+ ds0.X_E * btrue .+ 0.001 .* randn(rng, nconf)
        ds = SCEDataset(basis, configs, energies)
        f = fit(SCEFit, ds, GroupAdaptiveRidge(basis; lambda = 1e-3))
        @test f.support === nothing
        fr = refit(f)
        @test fr.support isa Vector{Int} && issorted(fr.support)
        @test_throws ArgumentError effective_dof(fr)
        @test_throws ArgumentError gcv(fr)
        # the refusal names the situation, not a dimension mismatch
        err = try
            gcv(fr)
        catch e
            e
        end
        @test err isa ArgumentError && occursin("refit", err.msg)
        # refit refuses a GroupAdaptiveRidge up front (its column_groups index the
        # full design, not the support)
        @test_throws ArgumentError refit(f, GroupAdaptiveRidge(basis; lambda = 1e-3))
    end

    @testset "the unpenalized rank cut is min(size), not max(size)" begin
        # `_rank_df` feeds `effective_dof`/`gcv` for every OLS fit and every λ = 0
        # point of a λ-path. Engineer a marginal singular value that sits ABOVE the
        # documented cut (`min(size)·eps·σ₁`) and watch a `max`-based cut drop it as
        # soon as the row count exceeds the column count.
        # [Backported from SLCE.jl 896180e.]
        p = 6
        for n in (10, 1000)
            U = qr(randn(MersenneTwister(4), n, p)).Q[:, 1:p]
            V = qr(randn(MersenneTwister(5), p, p)).Q
            svals = [1.0, 0.5, 0.25, 0.1, 0.05, 5 * p * eps(Float64)]
            Xm = U * Diagonal(svals) * V'
            @test SCEFitting._rank_df(Xm) == Float64(p)
            @test SCEFitting._rank_df(Xm) == Float64(rank(Xm))   # agrees with LinearAlgebra
        end
    end

    # --- λ path + Pareto selection --------------------------------------------------
    @testset "_select_pareto rule" begin
        sp = SCEFitting._select_pareto
        scores = [1.0, 1.02, 1.5]
        costs = [100.0, 40.0, 10.0]
        @test sp(scores, costs, 0.05) == 2      # 1.02 within 5% of 1.0, cheaper
        @test sp(scores, costs, 0.6) == 3       # 1.5 admitted too, cheapest wins
        @test sp(scores, costs, 0.0) == 1       # only the exact minimum eligible
        @test sp([1.0, 1.0], [5.0, 5.0], 0.1) == 1          # cost tie → larger λ
        @test sp([Inf, 1.0, 1.1], [1.0, 9.0, 8.0], 0.2) == 3  # Inf never eligible
        @test_throws ArgumentError sp([Inf, Inf], [1.0, 2.0], 0.1)
        @test_throws ArgumentError sp(scores, costs, -0.1)
        @test_throws DimensionMismatch sp(scores, [1.0], 0.1)
    end

    # group-sparse in-span fixture shared by the path tests
    rngp = MersenneTwister(41)
    nconf_p = 40
    configs_p = [Matrix(randcfg(rngp, 2)) for _ = 1:nconf_p]
    ds0_p = SCEDataset(basis, configs_p, zeros(nconf_p))
    labels_p = salc_groups(basis)
    btrue_p = zeros(n_salcs(basis))
    btrue_p[labels_p .== 1] .= 0.4
    btrue_p[labels_p .== 3] .= -0.25
    energies_p = 0.3 .+ ds0_p.X_E * btrue_p .+ 0.005 .* randn(rngp, nconf_p)
    ds_p = SCEDataset(basis, configs_p, energies_p)
    est_p = GroupAdaptiveRidge(basis; lambda = 1.0)   # template λ is ignored

    @testset "select_fit: warm-started path matches cold fits" begin
        # ÷ nconf_p: the energy-only branch now whitens by 1/√n_E (the objective is
        # the MSE at every weight), so a λ that used to act on the SSE Gram acts
        # n_E times harder — the grid moves with the convention to keep the same
        # regularization regime. [Backported from SLCE.jl 572cbe0.]
        lams = 10.0 .^ range(0, -6; length = 6) ./ nconf_p
        path = select_fit(ds_p, est_p; lambdas = lams, criterion = :gcv)
        @test path.lambda == sort(unique(Float64.(lams)); rev = true)
        # warm-started path entries reproduce cold single-λ fits: both converge to the
        # same fixed point but stop within the IRLS tol (1e-6 on coefficients), so the
        # scores agree to a few multiples of that tolerance, not exactly
        for i in (2, 4)
            est_i = GroupAdaptiveRidge(est_p.column_groups, est_p.group_weights;
                                       lambda = path.lambda[i])
            fc = fit(SCEFit, ds_p, est_i)
            @test isapprox(path.score[i], gcv(fc); rtol = 1e-4)
            @test isapprox(path.edof[i], effective_dof(fc); rtol = 1e-4)
        end
        # the selected fit is the cold fit at the selected λ, byte-for-byte
        est_s = GroupAdaptiveRidge(est_p.column_groups, est_p.group_weights;
                                   lambda = path.lambda[path.selected])
        @test path.fit.jphi == fit(SCEFit, ds_p, est_s).jphi
        # the selected row's score/edof are re-derived from that cold fit, so the
        # displayed row is self-consistent (not the warm path value)
        @test path.score[path.selected] ≈ gcv(path.fit) rtol = 1e-12
        @test path.edof[path.selected] ≈ effective_dof(path.fit) rtol = 1e-12
        # alive count grows (weakly) as λ decreases along the descending path
        @test all(diff(path.n_alive) .>= 0)
    end

    @testset "select_fit: end-to-end selection, tables, and validation" begin
        # ÷ nconf_p: same grid rescale as above (energy-only 1/√n_E whitening)
        lams = 10.0 .^ range(1, -7; length = 10) ./ nconf_p
        path = select_fit(ds_p, est_p; lambdas = lams, criterion = :gcv)
        np = length(path.lambda)
        @test length(path.score) == length(path.edof) == np
        @test length(path.n_alive) == length(path.cost) == np
        @test 1 <= path.selected <= np
        @test isfinite(path.score[path.selected])
        # the relative alive floor discriminates: the cost column varies along the
        # path (a flat column would mean the cost axis is vacuous)
        @test path.cost[1] < path.cost[end]
        @test path.threshold > 0
        @test minimum(path.n_alive) < maximum(path.n_alive)
        # widening the accuracy tolerance can only cheapen the selection
        pwide = select_fit(ds_p, est_p; lambdas = lams, criterion = :gcv, delta = 0.5)
        @test pwide.cost[pwide.selected] <= path.cost[path.selected]
        # de-bias on exactly the reported support and recover the data
        fr = refit(path.fit; threshold = path.threshold)
        @test r2_energy(fr) > 0.95
        # predicted cost at the selected point matches an independent recomputation
        # from the returned fit and threshold
        X, _, _, _, _ = SCEFitting._assemble_problem(ds_p, 0.0)
        cn = [norm(X[:, j]) for j = 1:size(X, 2)]
        cg = group_costs(basis, labels_p)
        aliveg = unique(labels_p[[j for j in eachindex(path.fit.jphi)
                                  if abs(path.fit.jphi[j]) * cn[j] > path.threshold]])
        @test path.cost[path.selected] == sum(cg[aliveg])
        @test path.n_alive[path.selected] == length(aliveg)
        # Tables view
        cols = Tables.columntable(path)
        @test keys(cols) == (:lambda, :score, :edof, :n_alive, :cost, :selected)
        @test count(cols.selected) == 1
        # show smoke
        @test occursin("← selected", sprint(show, MIME"text/plain"(), path))

        # an absolute threshold reproduces refit's rule verbatim; threshold = 0
        # counts every group alive (documented degenerate case)
        p0 = select_fit(ds_p, est_p; lambdas = [1.0, 1e-3], threshold = 0.0)
        @test all(==(maximum(labels_p)), p0.n_alive)
        @test p0.threshold === 0.0

        # grouped-CV criterion runs and selects a finite score
        pcv = select_fit(ds_p, est_p; lambdas = 10.0 .^ range(0, -5; length = 5),
                         criterion = :cv, nfolds = 4)
        @test all(isnan, pcv.edof)
        @test isfinite(pcv.score[pcv.selected])

        # validation errors
        @test_throws ArgumentError select_fit(ds_p, est_p; lambdas = Float64[])
        @test_throws ArgumentError select_fit(ds_p, est_p; lambdas = [-1.0])
        @test_throws ArgumentError select_fit(ds_p, est_p; lambdas = [1.0],
                                              criterion = :bogus)
        @test_throws ArgumentError select_fit(ds_p, est_p; lambdas = [1.0],
                                              delta = -0.1)
        @test_throws ArgumentError select_fit(ds_p, est_p; lambdas = [1.0],
                                              costs = [1.0])
        @test_throws ArgumentError select_fit(ds_p, est_p; lambdas = [1.0], nfolds = 1)
        @test_throws ArgumentError select_fit(ds_p, est_p; lambdas = [1.0],
                                              threshold = -1.0)
        @test_throws ArgumentError select_fit(ds_p, est_p; lambdas = [1.0],
                                              torque_weight = 1.5)
        @test_throws ArgumentError select_fit(ds_p, est_p; lambdas = [1.0],
                                              torque_weight = 0.5)   # no torque data
        # too few resampling units for grouped CV
        ds_tiny = SCEDataset(basis, configs_p[1:5], energies_p[1:5])
        @test_throws ArgumentError select_fit(ds_tiny, est_p; lambdas = [1.0, 0.1],
                                              criterion = :cv)
    end

    @testset "select_support: threshold-swept refit front" begin
        ds_tr2 = ds_p[1:30]
        ds_ev = ds_p[31:40]                       # held-out evaluation slice
        f = fit(SCEFit, ds_tr2, GroupAdaptiveRidge(basis; lambda = 1e-4))
        sp = select_support(f; npoints = 8, evalset = ds_ev)

        nt = length(sp.threshold)
        @test issorted(sp.threshold; rev = true)                  # sparsest first
        @test length(sp.n_alive) == length(sp.cost) == nt
        @test length(sp.score) == length(sp.rmse_energy) == nt
        @test issorted(sp.n_alive) && issorted(sp.cost)           # nested supports
        @test sp.threshold[end] == 0.0                            # full-support anchor
        @test sp.n_alive[end] == maximum(salc_groups(basis))
        @test minimum(sp.n_alive) < maximum(sp.n_alive)           # front discriminates
        @test all(isnan, sp.rmse_torque)                          # energy-only evalset
        @test sp.score ≈ sp.rmse_energy .^ 2                      # w = 0 objective
        @test 1 <= sp.selected <= nt
        @test isfinite(sp.score[sp.selected])
        # the returned fit is the refit at the selected threshold, byte-for-byte
        @test sp.fit.jphi == refit(f; threshold = sp.threshold[sp.selected]).jphi
        # in-span group-sparse target: the selected sparse refit generalizes, and
        # beats the full-support refit (which is an underdetermined min-norm OLS
        # here, n = 30 rows < p = 44 columns — it overfits the held-out slice)
        @test sp.rmse_energy[sp.selected] < 1e-2
        @test sp.rmse_energy[sp.selected] < sp.rmse_energy[end]
        # Pareto: widening delta can only cheapen the selection
        spw = select_support(f; npoints = 8, evalset = ds_ev, delta = 0.5)
        @test spw.cost[spw.selected] <= sp.cost[sp.selected]

        # explicit thresholds vector; single full-support point
        sp1 = select_support(f; thresholds = [0.0])
        @test length(sp1.threshold) == 1 && sp1.selected == 1

        # Tables + show smoke
        cols = Tables.columntable(sp)
        @test keys(cols) == (:threshold, :n_alive, :cost, :score, :rmse_energy,
                             :rmse_torque, :selected)
        @test count(cols.selected) == 1
        @test occursin("← selected", sprint(show, MIME"text/plain"(), sp))

        # coupled-site pin: the selected row re-derives from the returned refit —
        # groups with a nonzero de-biased coefficient are exactly the alive groups
        lab_b = salc_groups(basis)
        cg_b = group_costs(basis, lab_b)
        ag = unique(lab_b[findall(!=(0.0), sp.fit.jphi)])
        @test length(ag) == sp.n_alive[sp.selected]
        @test sum(cg_b[ag]) == sp.cost[sp.selected]

        # the auto grid: midpoints between distinct magnitudes + the 0.0 anchor;
        # exact ties collapse (tied groups die together, no empty-support point)
        st = SCEFitting._support_thresholds
        @test st(5, [3.0, 2.0, 1.0]) == [2.5, 1.5, 0.0]
        @test st(4, [2.0]) == [0.0]                       # G = 1 → anchor only
        @test st(9, [1.0, 1.0, 0.5]) == [0.75, 0.0]       # tie: no t = 1.0 point
        @test_throws ArgumentError st(1, [1.0])

        # torque co-fit: rmse_torque populated; energy-only evalset rejected
        rngs = MersenneTwister(53)
        cfg_s = [Matrix(randcfg(rngs, 2)) for _ = 1:12]
        en_s = randn(rngs, 12)
        tq_s = [randn(rngs, 3, 2) for _ = 1:12]
        ds_cofit = SCEDataset(small, cfg_s, en_s, tq_s)
        glab = salc_groups(small)
        fco = fit(SCEFit, ds_cofit, GroupAdaptiveRidge(glab, ones(maximum(glab));
                                                       lambda = 1e-3);
                  torque_weight = 0.4)
        spc = select_support(fco; npoints = 3)
        @test all(isfinite, spc.rmse_torque)
        @test length(spc.threshold) == 1        # single-group basis: grid collapses
        @test_throws ArgumentError select_support(fco;
                                                  evalset = SCEDataset(small, cfg_s, en_s))

        # validation errors
        @test_throws ArgumentError select_support(f; delta = -0.1)
        @test_throws ArgumentError select_support(f; thresholds = Float64[])
        @test_throws ArgumentError select_support(f; thresholds = [-1.0])
        @test_throws ArgumentError select_support(f; npoints = 1)
        # an Integer `thresholds` used to silently mean "n-point grid"; the count
        # and the explicit vector are separate keywords now, so it is a TypeError
        # [Backported behavior from SLCE.jl 2259c54.]
        @test_throws TypeError select_support(f; thresholds = 8)
        @test_throws ArgumentError select_support(f; costs = [1.0])
        @test_throws ArgumentError select_support(f; estimator = PrecomputedPilot(
            zeros(n_salcs(basis))))
        ds_other = SCEDataset(small, [Matrix(randcfg(rngs, 2)) for _ = 1:4],
                              randn(rngs, 4))
        @test_throws ArgumentError select_support(f; evalset = ds_other)
    end

    @testset "cross_validate" begin
        # energy-only fixture: 12 configs so nfolds = 3 needs no reduction
        rngc = MersenneTwister(67)
        nconf_c = 12
        cfg_c = [Matrix(randcfg(rngc, 2)) for _ = 1:nconf_c]
        ds0_c = SCEDataset(small, cfg_c, zeros(nconf_c))
        en_c = 0.2 .+ ds0_c.X_E * fill(0.3, n_salcs(small)) .+
               0.01 .* randn(rngc, nconf_c)
        ds_c = SCEDataset(small, cfg_c, en_c)

        cv = cross_validate(ds_c, OLS(); nfolds = 3)
        @test cv.nfolds == 3 && cv.seed == 1 && cv.torque_weight == 0.0
        @test sum(cv.n_holdout) == nconf_c
        @test all(isfinite, cv.score) && all(isfinite, cv.rmse_energy)
        @test all(isnan, cv.rmse_torque) && isnan(cv.pooled_rmse_torque)
        # energy-only: score = MSE_E (rmse stores the sqrt, so equality is to 1 ulp)
        @test cv.score ≈ cv.rmse_energy .^ 2 rtol = 1e-14
        @test cv.pooled_score ≈ cv.pooled_rmse_energy^2 rtol = 1e-14

        # manual reconstruction of one fold: same folds, same fit, same residuals
        folds_c = SCEFitting._grouped_folds(collect(1:nconf_c), 3, 1)
        ho1 = findall(==(1), folds_c)
        f1 = fit(SCEFit, ds_c[findall(!=(1), folds_c)], OLS())
        h1 = ds_c[ho1]
        @test cv.n_holdout[1] == length(ho1)
        @test cv.rmse_energy[1] ≈
              sqrt(mean(abs2, h1.y_E .- (f1.j0 .+ h1.X_E * f1.jphi))) rtol = 1e-12

        # pooled MSE = holdout-size-weighted mean of the per-fold MSEs
        @test cv.pooled_rmse_energy^2 ≈
              sum(cv.n_holdout .* cv.rmse_energy .^ 2) / nconf_c rtol = 1e-12

        # deterministic in the seed; a different seed re-deals the folds
        @test cross_validate(ds_c, OLS(); nfolds = 3).rmse_energy == cv.rmse_energy
        @test cross_validate(ds_c, OLS(); nfolds = 3, seed = 2).rmse_energy !=
              cv.rmse_energy

        # torque co-fit fixture: both error axes reported, even at torque_weight = 0
        tq_c = [randn(rngc, 3, 2) for _ = 1:nconf_c]
        ds_ct = SCEDataset(small, cfg_c, en_c, tq_c)
        cvt = cross_validate(ds_ct, OLS(); torque_weight = 0.4, nfolds = 3)
        @test all(isfinite, cvt.rmse_torque) && isfinite(cvt.pooled_rmse_torque)
        @test cvt.score ≈ 0.6 .* cvt.rmse_energy .^ 2 .+ 0.4 .* cvt.rmse_torque .^ 2
        @test cvt.pooled_score ≈ 0.6 * cvt.pooled_rmse_energy^2 +
                                 0.4 * cvt.pooled_rmse_torque^2
        cvt0 = cross_validate(ds_ct, OLS(); torque_weight = 0.0, nfolds = 3)
        @test all(isfinite, cvt0.rmse_torque)           # measured despite w = 0
        @test cvt0.score ≈ cvt0.rmse_energy .^ 2 rtol = 1e-14

        # works with a penalized estimator (the intended comparison use)
        glab_c = salc_groups(small)
        est_c = GroupAdaptiveRidge(glab_c, ones(maximum(glab_c)); lambda = 1e-4)
        cvg = cross_validate(ds_c, est_c; nfolds = 3)
        @test all(isfinite, cvg.score)

        # fold reduction: 7 configs support only 2 folds of >= 3 configs
        ds7 = ds_c[1:7]
        cv7 = @test_logs (:warn, r"reducing CV folds") match_mode = :any begin
            cross_validate(ds7, OLS(); nfolds = 5)
        end
        @test cv7.nfolds == 2 && sum(cv7.n_holdout) == 7
        # minimum valid size: 6 configs at nfolds = 2 succeeds without a warning
        cv6 = @test_logs cross_validate(ds_c[1:6], OLS(); nfolds = 2)
        @test cv6.nfolds == 2 && sum(cv6.n_holdout) == 6

        # Tables + show smoke
        tbl = Tables.columntable(cv)
        @test tbl.fold == [1, 2, 3] && tbl.rmse_energy == cv.rmse_energy
        @test occursin("pooled", sprint(show, MIME"text/plain"(), cvt))
        @test occursin("rmse_T", sprint(show, MIME"text/plain"(), cvt))
        @test !occursin("rmse_T", sprint(show, MIME"text/plain"(), cv))

        # validation errors
        @test_throws ArgumentError cross_validate(ds_c, OLS(); torque_weight = 1.5)
        @test_throws ArgumentError cross_validate(ds_c, OLS(); torque_weight = 0.5)
        @test_throws ArgumentError cross_validate(ds_c, OLS(); nfolds = 1)
        @test_throws ArgumentError cross_validate(ds_c[1:5], OLS())
        @test_throws ArgumentError cross_validate(ds_c, PrecomputedPilot(
            zeros(n_salcs(small))))
        @test_throws ArgumentError cross_validate(ds_c, AdaptiveLasso(
            pilot = PrecomputedPilot(zeros(n_salcs(small)))))
    end
end
