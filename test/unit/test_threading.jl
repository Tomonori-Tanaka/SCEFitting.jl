using Test
using SCEFitting
using LinearAlgebra
using Random

# The design-matrix assembly (`_design_energy` / `_design_torque`) and the vector
# `predict_energy` / `predict_torque` forms are multithreaded over independent
# columns / slots. These checks pin the threaded output to a race-free serial
# reference (and to the independent scalar predict path): they hold at any thread
# count and would flag a data race when run with `julia -t N>1`.

@testset "threaded assembly / prediction equals serial" begin
    @info "running with $(Threads.nthreads()) thread(s)"
    Threads.nthreads() == 1 &&
        @warn "threading tests run serial; launch `julia -t N>1` to exercise the parallel path"

    lat = Lattice(Matrix(3.0 * I(3)))
    crystal = Crystal(lat, [0.2 -0.2; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
    interaction = BasisSpec(; nbody = 2, cutoff = 1.5, lmax = [2], isotropy = true)
    basis = SCEBasis(crystal, interaction)
    salcs = basis.salc_basis.salcs
    m = length(salcs)
    nat = 2

    rng = MersenneTwister(7)
    configs = [randcfg(rng, nat) for _ = 1:64]
    cfgs = [Matrix{Float64}(c) for c in configs]

    @testset "_design_energy: threaded == serial double loop" begin
        X = SCEFitting._design_energy(basis, cfgs)             # threaded
        Xref = Matrix{Float64}(undef, length(cfgs), m)
        for j = 1:m, i in eachindex(cfgs)                      # race-free serial
            Xref[i, j] = SCEFitting.evaluate_salc(salcs[j], cfgs[i])
        end
        @test X == Xref                                        # exact: same kernel, deterministic
        @test SCEFitting._design_energy(basis, cfgs) == X      # idempotent across calls
    end

    @testset "penalty_metric: threaded == serial double loop" begin
        # The metric is built column-parallel inside a constructor, and it enters
        # every penalized coefficient, so it belongs to the same bitwise contract as
        # the design kernels: a schedule-dependent reduction here would move fits.
        K = 128
        cfgs_ref = SCEFitting._reference_configs(nat, K, 1)
        for w in (0.0, 1.0)
            got = penalty_metric(basis; torque_weight = w, nconfig = K, seed = 1)
            ref = Vector{Float64}(undef, m)
            for j = 1:m                                      # race-free serial
                s1 = 0.0
                s2 = 0.0
                st = 0.0
                for c in cfgs_ref
                    if w < 1.0
                        phi = SCEFitting.evaluate_salc(salcs[j], c)
                        s1 += phi
                        s2 += phi * phi
                    end
                    if w > 0.0
                        G = zeros(3, nat)
                        SCEFitting.accumulate_grad!(G, salcs[j], c, 1.0,
                                                    SCEFitting.SALCScratch())
                        for a = 1:nat
                            ea = view(c, :, a)
                            ga = view(G, :, a)
                            t = [ga[2] * ea[3] - ga[3] * ea[2],
                                 ga[3] * ea[1] - ga[1] * ea[3],
                                 ga[1] * ea[2] - ga[2] * ea[1]]
                            st += sum(abs2, t)
                        end
                    end
                end
                varj = max(0.0, s2 / K - (s1 / K)^2)
                ref[j] = (1 - w) * varj + w * (st / K) / (3 * nat)
            end
            @test got ≈ ref rtol = 1e-12
            @test penalty_metric(basis; torque_weight = w, nconfig = K, seed = 1) == got
        end
    end

    # Fit a model so the columns carry real coefficients for the cross-checks.
    y = 0.7 .+ SCEFitting._design_energy(basis, cfgs) * randn(MersenneTwister(3), m)
    ds = SCEDataset(basis, configs, y)
    model = SCEPredictor(fit(SCEFit, ds, OLS()))
    jphi = coef(model)

    @testset "_design_torque: threaded X_T·jϕ == scalar predict_torque" begin
        X_T = SCEFitting._design_torque(basis, cfgs)           # threaded
        @test X_T == SCEFitting._design_torque(basis, cfgs)    # idempotent
        # Independent path: assemble the full-model torque from the scalar kernel,
        # flattened config-major / atom-major / xyz, and compare to X_T·jϕ.
        block = 3 * nat
        ref = Vector{Float64}(undef, length(cfgs) * block)
        for ci in eachindex(cfgs)
            τ = predict_torque(model, cfgs[ci])                # scalar (all SALCs)
            rb = block * (ci - 1)
            for a = 1:nat, d = 1:3
                ref[rb + 3 * (a - 1) + d] = τ[d, a]
            end
        end
        @test isapprox(X_T * jphi, ref; atol = 1e-12, rtol = 0)
    end

    @testset "vector predict_* (threaded) == scalar map (serial)" begin
        @test predict_energy(model, configs) == [predict_energy(model, c) for c in configs]
        @test predict_torque(model, configs) == [predict_torque(model, c) for c in configs]
    end
end
