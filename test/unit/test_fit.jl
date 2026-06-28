using Test
using MagestyRebuild
using LinearAlgebra
using Random

struct _DummyEstimator <: MagestyRebuild.AbstractEstimator end

function randcfg(rng, nat)
    M = Matrix{Float64}(undef, 3, nat)
    for a = 1:nat
        v = randn(rng, 3)
        M[:, a] = v / norm(v)
    end
    return M
end

@testset "fit / predict" begin
    lat = Lattice(Matrix(3.0 * I(3)))
    crystal = Crystal(lat, [0.2 -0.2; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
    interaction = Interaction(; nbody = 2, pair_cutoff = 1.5, lmax = [2], isotropy = true)
    basis = SCEBasis(crystal, interaction)          # NoSymmetry backend (no Spglib needed)
    m = length(basis.salcs)
    @test m > 0

    rng = MersenneTwister(1)
    configs = [randcfg(rng, 2) for _ = 1:40]

    # synthetic energies that lie exactly in the SALC span
    ds0 = SCEDataset(basis, configs, zeros(length(configs)))
    true_jphi = randn(rng, m)
    true_j0 = 0.7
    y = true_j0 .+ ds0.X_E * true_jphi
    ds = SCEDataset(basis, configs, y)

    @testset "OLS recovers an in-span target" begin
        f = fit(SCEFit, ds, OLS())
        @test r2_energy(f) ≈ 1.0 atol = 1e-9
        @test rmse_energy(f) < 1e-6
        @test isapprox(predict_energy(f, configs), y; atol = 1e-7)
        @test nobs(f) == 40
        @test length(coef(f)) == m
    end

    @testset "Ridge(0) ≈ OLS; Ridge(λ) shrinks" begin
        f = fit(SCEFit, ds, OLS())
        f0 = fit(SCEFit, ds, Ridge(lambda = 0.0))
        @test isapprox(coef(f0), coef(f); atol = 1e-6)
        fλ = fit(SCEFit, ds, Ridge(lambda = 10.0))
        @test norm(coef(fλ)) <= norm(coef(f)) + 1e-9
    end

    @testset "unknown estimator → friendly error" begin
        @test_throws ErrorException solve_coefficients(_DummyEstimator(), ds.X_E, ds.y_E)
    end

    # The GLMNet-backed estimators are constructed and validated in the core (the
    # type is dependency-free); the actual solve needs `using GLMNet` and is exercised
    # in test/glmnet/. Without the extension loaded, a friendly error points there.
    @testset "ElasticNet / Lasso: construction, validation, deferred backend" begin
        @test Lasso() === ElasticNet(; alpha = 1.0)
        @test Lasso(lambda = 0.1).alpha == 1.0
        en = ElasticNet(; alpha = 0.4, lambda = 0.2, select = :lambda_1se, nfolds = 5, seed = 9)
        @test en.alpha == 0.4 && en.lambda == 0.2 && en.select === :lambda_1se
        @test ElasticNet().lambda === nothing                       # default: choose by CV
        @test !MagestyRebuild.islinear(Lasso())                     # not a closed-form estimator

        @test_throws ArgumentError ElasticNet(; alpha = 1.5)
        @test_throws ArgumentError ElasticNet(; alpha = -0.1)
        @test_throws ArgumentError ElasticNet(; lambda = -1.0)
        @test_throws ArgumentError ElasticNet(; select = :bogus)
        @test_throws ArgumentError ElasticNet(; nfolds = 1)

        # GLMNet is not loaded in this environment ⇒ the fallback fires.
        @test_throws ErrorException solve_coefficients(Lasso(), ds.X_E, ds.y_E)
    end
end
