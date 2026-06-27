using Test
using MagestyRebuild
using LinearAlgebra
using StaticArrays
using Random

# Self-contained config generator (fresh draw per column, so neighboring spins
# differ and odd-Lf / anisotropic channels do not vanish).
function _randcfg(rng, nat)
    M = Matrix{Float64}(undef, 3, nat)
    for a = 1:nat
        v = randn(rng, 3)
        M[:, a] = v / norm(v)
    end
    return M
end

# A random proper rotation (det = +1) from a QR of a Gaussian matrix.
function _rand_rotation(rng)
    Q, R = qr(randn(rng, 3, 3))
    Q = Matrix(Q) * Diagonal(sign.(diag(R)))   # fix the QR sign ambiguity
    det(Q) < 0 && (Q[:, 1] .*= -1)             # make it proper
    return Q
end

# Central-difference torque from the energy surface: τ_a = e_a × ∇_tan E, with the
# perturbation renormalized onto the sphere so the finite difference reconstructs
# the tangent gradient (the radial part drops out of the cross product anyway).
function _torque_fd(model, config; δ = 1e-6)
    nat = size(config, 2)
    T = Matrix{Float64}(undef, 3, nat)
    for a = 1:nat
        e = config[:, a]
        g = zeros(3)
        for d = 1:3
            ep = copy(config)
            em = copy(config)
            step = zeros(3); step[d] = δ
            ep[:, a] = (e .+ step) ./ norm(e .+ step)
            em[:, a] = (e .- step) ./ norm(e .- step)
            g[d] = (predict_energy(model, ep) - predict_energy(model, em)) / (2δ)
        end
        T[:, a] = cross(SVector{3}(e), SVector{3}(g))
    end
    return T
end

@testset "torque" begin
    rng = MersenneTwister(7)
    lat = Lattice(Matrix(3.0 * I(3)))
    crystal = Crystal(lat, [0.2 -0.2; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])

    @testset "predict_torque = e × ∇E (finite differences, anisotropic)" begin
        interaction = Interaction(; nbody = 2, pair_cutoff = 1.5, lmax = [2], isotropy = false)
        basis = SCEBasis(crystal, interaction)   # NoSymmetry (P1): exercises Lf > 0
        m = length(basis.salcs)
        @test m > 0
        model = SCEModel(basis, 0.0, randn(rng, m), basis.salcs.keys)
        for _ = 1:8
            c = _randcfg(rng, 2)
            @test isapprox(predict_torque(model, c), _torque_fd(model, c); atol = 1e-5)
        end
    end

    @testset "global-rotation equivariance (isotropic model)" begin
        interaction = Interaction(; nbody = 2, pair_cutoff = 1.5, lmax = [1], isotropy = true)
        basis = SCEBasis(crystal, interaction)
        m = length(basis.salcs)
        model = SCEModel(basis, 0.3, randn(rng, m), basis.salcs.keys)
        for _ = 1:5
            c = _randcfg(rng, 2)
            R = _rand_rotation(rng)
            τ = predict_torque(model, c)
            τR = predict_torque(model, R * c)
            @test isapprox(τR, R * τ; atol = 1e-10)   # τ is an axial vector under proper R
        end
    end

    @testset "energy+torque co-fit recovers an in-span model" begin
        interaction = Interaction(; nbody = 2, pair_cutoff = 1.5, lmax = [2], isotropy = false)
        basis = SCEBasis(crystal, interaction)
        m = length(basis.salcs)
        configs = [_randcfg(rng, 2) for _ = 1:60]

        true_jphi = randn(rng, m)
        true_j0 = 0.4
        skel = SCEDataset(basis, configs, zeros(length(configs)))   # to read X_E / X_T
        energies = true_j0 .+ skel.X_E * true_jphi
        # per-config torque targets, consistent with the same model
        model0 = SCEModel(basis, true_j0, true_jphi, basis.salcs.keys)
        torques = [predict_torque(model0, c) for c in configs]

        ds = SCEDataset(basis, configs, energies, torques)
        @test has_torque(ds)

        f = fit(SCEFit, ds, OLS(); torque_weight = 0.3)
        @test isapprox(coef(f), true_jphi; atol = 1e-6)
        @test isapprox(intercept(f), true_j0; atol = 1e-6)
        @test r2_energy(f) ≈ 1.0 atol = 1e-9
        @test r2_torque(f) ≈ 1.0 atol = 1e-9
        @test rmse_torque(f) < 1e-6
        # predicted torques match the targets
        @test isapprox(predict_torque(f, configs[1]), torques[1]; atol = 1e-7)
    end

    @testset "pure-torque fit (weight = 1) still pins jϕ and recovers j0" begin
        interaction = Interaction(; nbody = 2, pair_cutoff = 1.5, lmax = [1], isotropy = true)
        basis = SCEBasis(crystal, interaction)
        m = length(basis.salcs)
        configs = [_randcfg(rng, 2) for _ = 1:50]
        true_jphi = randn(rng, m)
        true_j0 = -0.9
        model0 = SCEModel(basis, true_j0, true_jphi, basis.salcs.keys)
        energies = [predict_energy(model0, c) for c in configs]
        torques = [predict_torque(model0, c) for c in configs]
        ds = SCEDataset(basis, configs, energies, torques)
        f = fit(SCEFit, ds, OLS(); torque_weight = 1.0)
        @test isapprox(coef(f), true_jphi; atol = 1e-6)
        @test isapprox(intercept(f), true_j0; atol = 1e-6)   # j0 from energy data alone
    end

    @testset "error paths" begin
        interaction = Interaction(; nbody = 2, pair_cutoff = 1.5, lmax = [1], isotropy = true)
        basis = SCEBasis(crystal, interaction)
        configs = [_randcfg(rng, 2) for _ = 1:10]
        ds_e = SCEDataset(basis, configs, zeros(10))           # energy-only
        @test !has_torque(ds_e)
        @test_throws ArgumentError fit(SCEFit, ds_e, OLS(); torque_weight = 0.5)
        @test_throws ArgumentError r2_torque(fit(SCEFit, ds_e, OLS()))

        torques = [zeros(3, 2) for _ = 1:10]
        ds_t = SCEDataset(basis, configs, zeros(10), torques)
        @test_throws ArgumentError fit(SCEFit, ds_t, OLS(); torque_weight = 1.5)
        @test_throws ArgumentError fit(SCEFit, ds_t, OLS(); torque_weight = -0.1)
        # wrong torque block shape
        @test_throws ArgumentError SCEDataset(basis, configs, zeros(10), [zeros(3, 3) for _ = 1:10])
        # mismatched count
        @test_throws ArgumentError SCEDataset(basis, configs, zeros(10), [zeros(3, 2) for _ = 1:9])
    end
end
