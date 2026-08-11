using Test
using SCEFitting
using LinearAlgebra
using Random

# The training-data boundary enforces the spin-config contract the harmonic kernels
# assume (3 × n_atoms, finite, unit-norm columns). These checks pin that malformed
# input fails loudly at construction / prediction instead of silently biasing the fit.

@testset "data-boundary validation" begin
    lat = Lattice(Matrix(3.0 * I(3)))
    crystal = Crystal(lat, [0.2 -0.2; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
    interaction = BasisSpec(; nbody = 2, cutoff = 1.5, lmax = [2], isotropy = true)
    basis = SCEBasis(crystal, interaction)
    nat = n_atoms(crystal)

    rng = MersenneTwister(1)
    good = [randcfg(rng, nat) for _ = 1:8]
    energies = randn(rng, length(good))

    @testset "valid input constructs / predicts" begin
        ds = SCEDataset(basis, good, energies)
        @test nobs(fit(SCEFit, ds, OLS())) == length(good)
        model = SCEPredictor(fit(SCEFit, ds, OLS()))
        @test predict_energy(model, good) == [predict_energy(model, c) for c in good]
        @test length(predict_torque(model, good)) == length(good)
    end

    model = SCEPredictor(fit(SCEFit, SCEDataset(basis, good, energies), OLS()))

    @testset "SCEDataset rejects malformed configs" begin
        # non-unit column
        bad = deepcopy(good); bad[3][:, 1] .*= 1.5
        @test_throws ArgumentError SCEDataset(basis, bad, energies)
        # NaN
        bad = deepcopy(good); bad[2][1, 2] = NaN
        @test_throws ArgumentError SCEDataset(basis, bad, energies)
        # wrong atom count
        @test_throws DimensionMismatch SCEDataset(basis, [randcfg(rng, nat + 1)], [0.0])
        # wrong row count
        @test_throws ArgumentError SCEDataset(basis, [randn(rng, 2, nat)], [0.0])
        # config / energy length mismatch
        @test_throws DimensionMismatch SCEDataset(basis, good, energies[1:end-1])
    end

    @testset "torque dataset rejects malformed configs / lengths" begin
        torques = [zeros(3, nat) for _ in good]
        @test nobs(fit(SCEFit, SCEDataset(basis, good, energies, torques), OLS())) == length(good)
        bad = deepcopy(good); bad[1][:, 2] .= 0.0           # zero-norm column
        @test_throws ArgumentError SCEDataset(basis, bad, energies, torques)
        @test_throws ArgumentError SCEDataset(basis, good, energies, torques[1:end-1])
    end

    @testset "predict rejects malformed configs (scalar + batch)" begin
        @test_throws ArgumentError predict_energy(model, 2.0 .* good[1])
        @test_throws ArgumentError predict_torque(model, 2.0 .* good[1])
        @test_throws DimensionMismatch predict_energy(model, randcfg(rng, nat + 1))
        # batch path validates serially up front (clean error, not a TaskFailedException)
        bad = deepcopy(good); bad[4][:, 1] .*= 3.0
        @test_throws ArgumentError predict_energy(model, bad)
        @test_throws ArgumentError predict_torque(model, bad)
    end

    @testset "atol keyword loosens / tightens the unit-norm check" begin
        nearly = deepcopy(good); nearly[1][:, 1] .*= (1 + 1e-4)
        @test_throws ArgumentError SCEDataset(basis, nearly, energies)         # default 1e-6
        @test nobs(fit(SCEFit, SCEDataset(basis, nearly, energies; atol = 1e-3), OLS())) == length(good)
    end

    @testset "the atol cap binds the dataset door" begin
        # This door validates WITHOUT projecting, so an uncapped atol admits
        # whatever clears the widened band raw into the design: atol = 0.5 lets a
        # moment-scaled column (‖e‖ = 0.6) through and biases the design by tens
        # of percent. [Backported from SLCE.jl 54457ca, review M2.]
        scaled = deepcopy(good); scaled[1][:, 1] .*= 0.6
        err = try
            SCEDataset(basis, scaled, energies; atol = 0.5)
        catch e
            e
        end
        @test err isa ArgumentError && occursin("exceeds the hard cap", err.msg)
        @test_throws ArgumentError SCEDataset(basis, good, energies; atol = 2e-2)
        @test_throws ArgumentError SCEDataset(basis, good, energies; atol = -1.0)
        @test_throws ArgumentError SCEDataset(basis, good, energies; atol = NaN)
        # a capped widening is still accepted (rejects nothing legitimate)
        nearly = deepcopy(good); nearly[1][:, 1] .*= (1 + 1e-4)
        @test SCEDataset(basis, nearly, energies; atol = 1e-3) isa SCEDataset
    end

    @testset "a component outside [-1,1] is refused, and atol cannot license it" begin
        # The kernel's Legendre recursion has the HARD domain |e_z| ≤ 1, so this used
        # to pass the boundary and throw a DomainError from inside the threaded design
        # assembly, naming neither the config nor the atom. Tightening `atol` does not
        # help: ‖e‖ − 1 = 5e-9 clears a 1e-8 band and still throws.
        # [Backported from SLCE.jl 8ee6739.]
        pole = deepcopy(good); pole[1][:, 1] = [0.0, 0.0, 1 + 5e-9]
        @test abs(norm(pole[1][:, 1]) - 1) < 1e-8          # inside any sane atol band
        @test_throws ArgumentError SCEDataset(basis, pole, energies)
        @test_throws ArgumentError SCEDataset(basis, pole, energies; atol = 1e-3)
        @test_throws ArgumentError predict_energy(model, pole[1])
        # ... while an off-norm column whose components all sit inside the domain is
        # still accepted, so this rejects nothing the tolerance was meant to allow
        tilt = deepcopy(good); tilt[1][:, 1] = [0.6, 0.8, 0.0] .* (1 + 1e-7)
        @test maximum(abs, tilt[1][:, 1]) < 1
        @test SCEDataset(basis, tilt, energies) isa SCEDataset
    end
end

# The basis containers carry contracts every downstream consumer reads without
# re-checking: `SALCBasis` addresses design-matrix columns by key, and `SCEBasis`
# indexes `spec.lmax` / `spec.cutoff` by the crystal's species ids. Both used to be
# plain structs, so the field-wise call bypassed the validating outer constructor —
# which persistence actually takes. These pin the inner constructors.
# [Backported from SLCE.jl bde9ded.]
@testset "basis constructor validation" begin
    lat = Lattice(Matrix(3.0 * I(3)))
    cr = Crystal(lat, [0.2 -0.2; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
    spec = BasisSpec(; nbody = 2, cutoff = 1.5, lmax = [2], isotropy = true)
    basis = SCEBasis(cr, spec)
    sb = basis.salc_basis
    good = sb.salcs
    keys = sb.keys
    @test length(good) >= 2                     # teeth: the permutations below are real

    @testset "SALCBasis" begin
        @test SCEFitting.SALCBasis(good, keys).fingerprint == sb.fingerprint
        # the fingerprint is derived, so it cannot be supplied out of sync
        @test_throws MethodError SCEFitting.SALCBasis(good, keys, hash(keys))
        @test_throws MethodError SCEFitting.SALCBasis(good, keys, zero(UInt64))
        @test_throws DimensionMismatch SCEFitting.SALCBasis(good, keys[1:end-1])
        # keys must mirror the SALCs, in order
        shuffled = copy(keys); shuffled[1], shuffled[2] = shuffled[2], shuffled[1]
        @test_throws ArgumentError SCEFitting.SALCBasis(good, shuffled)
        # reversed: mirrors nothing and is not sorted
        @test_throws ArgumentError SCEFitting.SALCBasis(reverse(good), reverse(keys))
        # a strictly-decreasing pair that DOES mirror its SALCs is caught by the sort check
        rs = reverse(good)
        @test_throws ArgumentError SCEFitting.SALCBasis(rs, SCEFitting.SALCKey[s.key for s in rs])
        # duplicate column
        dup = [good[1], good[1]]
        @test_throws ArgumentError SCEFitting.SALCBasis(dup, SCEFitting.SALCKey[s.key for s in dup])
    end

    @testset "SCEBasis field-wise call cannot skip the species check" begin
        two = BasisSpec(; nbody = 2, cutoff = 1.5, lmax = [2, 2], isotropy = true)
        @test_throws ArgumentError SCEBasis(cr, basis.spacegroup, sb, two)
        @test_throws ArgumentError SCEBasis(cr, two)                # outer, fails early
        # the honest field-wise call still works — this rejects nothing legitimate
        @test SCEBasis(cr, basis.spacegroup, sb, spec) isa SCEBasis
    end
end
