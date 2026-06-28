# GLMNet-export validation: fit SCE coefficients with the Lasso / ElasticNet
# estimators and confirm the penalized solver honors the column-centered contract,
# recovers a sparse signal, and selects the penalty by cross-validation. Runs in a
# separate environment (heavy GLMNet / Fortran binary dependency), mirroring
# `test/sunny/` and `test/oracle/`. The `ElasticNet` / `Lasso` types and their
# argument validation are covered Estimator-free in `test/unit/test_fit.jl`; here we
# exercise the actual GLMNet solve through the `MagestyRebuildGLMNetExt` extension.

using Test
using Random
using LinearAlgebra
import MagestyRebuild
import GLMNet   # triggers MagestyRebuildGLMNetExt

const MR = MagestyRebuild

_rdir(rng) = (v = randn(rng, 3); v / norm(v))
randcfg(rng, nat) = reduce(hcat, (_rdir(rng) for _ = 1:nat))
avg(x) = sum(x) / length(x)

@testset "GLMNet estimators" begin
    lat = MR.Lattice(Matrix(3.0 * I(3)))
    crystal = MR.Crystal(lat, [0.2 -0.2; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
    inter = MR.Interaction(; nbody = 2, pair_cutoff = 1.5, lmax = [2], isotropy = false)
    basis = MR.SCEBasis(crystal, inter)          # NoSymmetry ⇒ a wide (44-column) basis
    m = MR.nsalc(basis)
    @test m > 10
    keys = basis.salcs.keys
    kiso = findfirst(k -> sort(collect(k.ls)) == [1, 1] && k.Lf == 0, keys)
    @test kiso !== nothing                       # the isotropic (Heisenberg) channel

    rng = MersenneTwister(7)
    configs = [randcfg(rng, 2) for _ = 1:80]

    # A dense in-span target (no noise) and the design matrix it lives in.
    ds0 = MR.SCEDataset(basis, configs, zeros(length(configs)))
    βdense = randn(rng, m)
    ydense = 0.4 .+ ds0.X_E * βdense

    @testset "tiny λ ≈ OLS, and the analytic-j0 contract holds" begin
        ds = MR.SCEDataset(basis, configs, ydense)
        fols = MR.fit(MR.SCEFit, ds, MR.OLS())
        fl = MR.fit(MR.SCEFit, ds, MR.Lasso(lambda = 1e-6))
        @test isapprox(MR.predict_energy(fols, configs), ydense; atol = 1e-8)   # OLS interpolates
        @test isapprox(MR.predict_energy(fl, configs), ydense; rtol = 2e-2)     # near-unpenalized
        # j0 is recovered analytically downstream regardless of the estimator, so the
        # in-sample mean prediction always equals the mean target.
        @test avg(MR.predict_energy(fl, configs)) ≈ avg(ydense)
    end

    # A sparse signal (only the isotropic column) plus noise: the regime the Lasso is for.
    ynoisy = 0.4 .+ 1.0 .* ds0.X_E[:, kiso] .+ 0.01 .* randn(rng, length(configs))
    dssparse = MR.SCEDataset(basis, configs, ynoisy)

    @testset "CV Lasso recovers the support and is sparse" begin
        f = MR.fit(MR.SCEFit, dssparse, MR.Lasso())              # CV, :lambda_min
        @test argmax(abs.(MR.coef(f))) == kiso                    # signal column dominates
        @test MR.r2_energy(f) > 0.95                              # and still fits
    end

    @testset ":lambda_1se shrinks more than :lambda_min and is sparse" begin
        fmin = MR.fit(MR.SCEFit, dssparse, MR.Lasso(select = :lambda_min))
        f1se = MR.fit(MR.SCEFit, dssparse, MR.Lasso(select = :lambda_1se))
        # λ_1se ≥ λ_min ⇒ more shrinkage; ‖β‖₁ is monotone non-increasing in λ.
        @test norm(MR.coef(f1se), 1) <= norm(MR.coef(fmin), 1) + 1e-8
        # With 44 columns and a 1-sparse signal, the parsimonious 1-SE model is sparse.
        @test count(iszero, MR.coef(f1se)) >= 1
    end

    @testset "seeded CV is reproducible" begin
        a = MR.coef(MR.fit(MR.SCEFit, dssparse, MR.Lasso(seed = 3)))
        b = MR.coef(MR.fit(MR.SCEFit, dssparse, MR.Lasso(seed = 3)))
        @test a == b
    end

    @testset "ElasticNet (α<1) differs from the Lasso and is less sparse" begin
        λ = 0.02
        fl = MR.fit(MR.SCEFit, dssparse, MR.Lasso(lambda = λ))
        fe = MR.fit(MR.SCEFit, dssparse, MR.ElasticNet(alpha = 0.3, lambda = λ))
        @test !isapprox(MR.coef(fe), MR.coef(fl); atol = 1e-8)
        @test count(iszero, MR.coef(fe)) <= count(iszero, MR.coef(fl))
    end

    @testset "energy+torque co-fit flows through GLMNet" begin
        truem = MR.SCEModel(basis, 0.4, βdense, keys)
        torqs = MR.predict_torque(truem, configs)
        dst = MR.SCEDataset(basis, configs, ydense, torqs)
        # explicit λ: both observables recovered
        f = MR.fit(MR.SCEFit, dst, MR.Lasso(lambda = 1e-5); torque_weight = 0.3)
        @test MR.r2_energy(f) > 0.9
        @test MR.r2_torque(f) > 0.9
        # CV path: exercises configuration-grouped folds (energy + torque rows of one
        # config stay together) and still fits both observables.
        fcv = MR.fit(MR.SCEFit, dst, MR.Lasso(); torque_weight = 0.3)
        @test MR.r2_energy(fcv) > 0.9
        @test MR.r2_torque(fcv) > 0.9
    end
end
