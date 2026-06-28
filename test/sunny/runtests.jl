# Sunny-export validation: build a Sunny.System from a fitted SCEModel and confirm
# its classical energy reproduces the SCE energy. Runs in a separate environment
# (heavy Sunny dependency), mirroring `test/oracle/`. The conversion math itself is
# already covered Sunny-free in `test/unit/test_sunny.jl`; here we check the actual
# Sunny.System assembly end-to-end.

using Test
using Random
using LinearAlgebra
using StaticArrays
import MagestyRebuild
import Sunny   # triggers MagestyRebuildSunnyExt

const MR = MagestyRebuild

_rdir(rng) = (v = randn(rng, 3); v / norm(v))
_rcfg(rng, n) = reshape(reduce(vcat, (_rdir(rng) for _ = 1:n)), 3, n)

# Max |Sunny energy − (predict_energy − j0)| over random configurations.
function energy_error(model, spins, mode; seed = 5, ntrial = 20)
    sys = MR.to_sunny(model; spins = spins, mode = mode)
    nat = MR.num_atoms(model.basis.crystal)
    rng = MersenneTwister(seed)
    me = 0.0
    for _ = 1:ntrial
        e = _rcfg(rng, nat)
        for i = 1:nat
            Sunny.set_dipole!(sys, e[:, i], (1, 1, 1, i))
        end
        me = max(me, abs(Sunny.energy(sys) - (MR.predict_energy(model, e) - model.j0)))
    end
    return me
end

@testset "Sunny export vs SCE" begin
    lat = MR.Lattice(Matrix(3.0 * I(3)))
    cr2 = MR.Crystal(lat, [0.2 -0.2; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
    rng = MersenneTwister(1)

    @testset "bilinear exchange: energy matches across modes and spins" begin
        b = MR.SCEBasis(cr2, MR.Interaction(; nbody = 2, pair_cutoff = 1.5, lmax = [1],
                                            isotropy = false))
        model = MR.SCEModel(b, 0.5, randn(rng, MR.nsalc(b)), b.salcs.keys)
        @test energy_error(model, 1.5, :dipole_uncorrected) < 1e-10
        @test energy_error(model, 1.0, :dipole) < 1e-10            # exchange is mode-independent
        @test energy_error(model, 2.5, :dipole) < 1e-10           # and spin-length-independent
    end

    @testset "single-ion anisotropy: classical energy matches" begin
        cr1 = MR.Crystal(lat, reshape([0.0, 0, 0], 3, 1), [1], ["Fe"])
        bo = MR.SCEBasis(cr1, MR.Interaction(; nbody = 1, pair_cutoff = 1.5, lmax = [2],
                                             isotropy = false))
        mo = MR.SCEModel(bo, 0.0, randn(rng, MR.nsalc(bo)), bo.salcs.keys)
        @test energy_error(mo, 1.5, :dipole_uncorrected) < 1e-10
    end

    @testset "pair + single-ion together (supported channels) match" begin
        # lmax=[2] also brings unsupported ls=[2,2] pairs; zero those so the exported
        # System and the model carry the same energy, then check consistency.
        b = MR.SCEBasis(cr2, MR.Interaction(; nbody = 2, pair_cutoff = 1.5, lmax = [2],
                                            isotropy = false))
        keys = b.salcs.keys
        jphi = randn(rng, MR.nsalc(b))
        for k in eachindex(keys)
            MR._classify_salc(keys[k].ls) === :unsupported && (jphi[k] = 0.0)
        end
        model = MR.SCEModel(b, 0.2, jphi, keys)
        @test energy_error(model, 1.5, :dipole_uncorrected) < 1e-10
    end

    @testset "unsupported channels trigger a warning" begin
        b = MR.SCEBasis(cr2, MR.Interaction(; nbody = 2, pair_cutoff = 1.5, lmax = [2],
                                            isotropy = false))
        model = MR.SCEModel(b, 0.0, ones(MR.nsalc(b)), b.salcs.keys)
        @test_logs (:warn,) match_mode = :any MR.to_sunny(model; spins = 1.5,
                                                          mode = :dipole_uncorrected)
    end

    @testset "spins as a per-species mapping" begin
        b = MR.SCEBasis(cr2, MR.Interaction(; nbody = 2, pair_cutoff = 1.5, lmax = [1],
                                            isotropy = true))
        model = MR.SCEModel(b, 0.0, [0.01], b.salcs.keys)
        @test energy_error(model, Dict("Fe" => 1.5), :dipole_uncorrected) < 1e-10
    end
end
