# The code-agnostic DFT data boundary (src/io/dftsource.jl): `SpinDatum`,
# `read_configs`, and `SCEDataset(basis, data/src)`. The headline gate is the torque
# convention `τ_a = m_a × B_a` (the physical / Landau–Lifshitz torque) — CLAUDE.md
# designates this file as the convention source every DFT adapter must match, so a
# closed-form sign check lives here, next to the definition. A silent sign flip
# would bias every energy+torque co-fit while all self-consistency tests still pass.

using Test
using SCEFitting
using LinearAlgebra
using Random

# A minimal in-memory source, exercising the `read_configs` extension seam the way a
# real adapter (e.g. SCETools.VASP.Oszicar) does.
struct _MemSource <: AbstractDFTSource
    data::Vector{SpinDatum}
end
SCEFitting.read_configs(s::_MemSource) = s.data

struct _EmptySource <: AbstractDFTSource end   # no read_configs method on purpose

@testset "DFT data boundary (SpinDatum / read_configs)" begin
    @testset "torque convention: τ = m × B, closed form" begin
        # m = m·x̂, B = B·ŷ ⇒ τ = m×B = mB·ẑ (and NOT −mB·ẑ: the old, pre-LL sign)
        moments = reshape([2.0, 0.0, 0.0], 3, 1)
        field = reshape([0.0, 0.5, 0.0], 3, 1)
        d = SpinDatum(-1.0, moments, field)
        @test d.torques[:, 1] ≈ [0.0, 0.0, 1.0]
        @test d.directions[:, 1] ≈ [1.0, 0.0, 0.0]
        @test d.magmoms[1] ≈ 2.0
        @test d.energy == -1.0
        # generic cross-product check on random data
        rng = MersenneTwister(11)
        mo = randn(rng, 3, 4)
        Bf = randn(rng, 3, 4)
        dr = SpinDatum(0.0, mo, Bf)
        for a = 1:4
            @test dr.torques[:, a] ≈ cross(mo[:, a], Bf[:, a]) atol = 1e-14
            @test dr.directions[:, a] ≈ mo[:, a] ./ norm(mo[:, a]) atol = 1e-14
        end
    end

    @testset "zero-moment placeholder: ẑ direction, zero torque" begin
        moments = [1.0 0.0; 0.0 0.0; 0.0 0.0]         # atom 2 is quenched
        field = [0.0 1.0; 1.0 1.0; 0.0 1.0]
        d = SpinDatum(0.0, moments, field)
        @test d.directions[:, 2] ≈ [0.0, 0.0, 1.0]     # the documented ẑ placeholder
        @test d.torques[:, 2] ≈ zeros(3) atol = 1e-15  # m = 0 ⇒ τ = 0 exactly
        @test d.magmoms[2] == 0.0
    end

    @testset "constructor / dataset validation throws" begin
        @test_throws ArgumentError SpinDatum(0.0, zeros(2, 3), zeros(2, 3))   # not 3 × n
        @test_throws ArgumentError SpinDatum(0.0, zeros(3, 2), zeros(3, 3))   # shape mismatch
        @test_throws ArgumentError read_configs(_EmptySource())               # no method
        # dataset-level guards
        lat = Lattice(Matrix(3.0 * I(3)))
        cr = Crystal(lat, [0.2 -0.2; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
        basis = SCEBasis(cr, BasisSpec(; nbody = 2, cutoff = 1.5, lmax = [1],
                                       isotropy = true))
        @test_throws ArgumentError SCEDataset(basis, SpinDatum[])             # empty data
        # all-zero torque targets with use_torque = true must fail loudly
        rng = MersenneTwister(3)
        nofield = [SpinDatum(0.1, randn(rng, 3, 2), zeros(3, 2)) for _ = 1:3]
        @test_throws ArgumentError SCEDataset(basis, nofield)
        ds = SCEDataset(basis, nofield; use_torque = false)                   # energy-only OK
        @test !has_torque(ds)
    end

    @testset "zero-moment guard: referenced atoms must stay magnetic" begin
        lat = Lattice(Matrix(3.0 * I(3)))
        # Fe + B; B is removed from the basis via lmax = 0 (non-magnetic species),
        # so only the Fe single-ion SALCs reference an atom (no Fe–Fe pair exists:
        # a single Fe atom has no self-pair under MinimumImage).
        cr = Crystal(lat, [0.2 -0.2; 0.0 0.0; 0.0 0.0], [1, 2], ["Fe", "B"])
        basis = SCEBasis(cr, BasisSpec(cr; nbody = 2, cutoff = 1.5,
                                       lmax = ["Fe" => 2, "B" => 0]))
        # quenched B is fine — the basis never reads it
        m_okay = [2.0 0.0; 0.0 0.0; 0.0 0.0]
        ds = SCEDataset(basis, [SpinDatum(0.0, m_okay, zeros(3, 2))];
                        use_torque = false)
        @test length(ds) == 1
        # quenched Fe is an error naming the config, atom, and species
        m_bad = [0.0 0.0; 0.0 1.0; 0.0 0.0]
        err = try
            SCEDataset(basis, [SpinDatum(0.0, m_okay, zeros(3, 2)),
                               SpinDatum(0.0, m_bad, zeros(3, 2))];
                       use_torque = false)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("config 2", err.msg) && occursin("(Fe)", err.msg)
        # atom-count mismatch is caught at the same boundary
        @test_throws DimensionMismatch SCEDataset(
            basis, [SpinDatum(0.0, zeros(3, 3) .+ 1.0, zeros(3, 3))];
            use_torque = false)
    end

    @testset "zero_moment_atol reaches the guard through the source path" begin
        # The convenience constructor `SCEDataset(basis, src)` used to silently drop
        # `zero_moment_atol`, so an adapter that built its `SpinDatum`s with a custom
        # tolerance could not align the referenced-moment guard with it.
        # [Backported behavior from SLCE.jl 011e3c9.]
        lat = Lattice(Matrix(3.0 * I(3)))
        cr = Crystal(lat, [0.2 -0.2; 0.0 0.0; 0.0 0.0], [1, 2], ["Fe", "B"])
        basis = SCEBasis(cr, BasisSpec(cr; nbody = 2, cutoff = 1.5,
                                       lmax = ["Fe" => 2, "B" => 0]))
        # an Fe moment far below the default 1e-10 guard, built with a matching
        # loose build tolerance so the stored direction is the real one
        m_tiny = [1.0e-12 0.0; 0.0 0.0; 0.0 0.0]
        src = _MemSource([SpinDatum(0.0, m_tiny, zeros(3, 2);
                                    zero_moment_atol = 1.0e-14)])
        # default guard: the referenced Fe atom reads as quenched — refused
        @test_throws ArgumentError SCEDataset(basis, src; use_torque = false)
        # the forwarded loose guard accepts the same source
        ds = SCEDataset(basis, src; use_torque = false, zero_moment_atol = 1.0e-14)
        @test length(ds) == 1
    end

    @testset "source → dataset round trip carries directions / energies / torques" begin
        lat = Lattice(Matrix(3.0 * I(3)))
        cr = Crystal(lat, [0.2 -0.2; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
        basis = SCEBasis(cr, BasisSpec(; nbody = 2, cutoff = 1.5, lmax = [1],
                                       isotropy = true))
        rng = MersenneTwister(5)
        data = [SpinDatum(0.1 * k, randn(rng, 3, 2), 0.1 .* randn(rng, 3, 2)) for k = 1:4]
        src = _MemSource(data)
        @test read_configs(src) === data
        ds = SCEDataset(basis, src)
        @test has_torque(ds)
        @test ds.y_E == [0.1 * k for k = 1:4]
        @test ds.configs == [d.directions for d in data]
        # flattened config-major / atom-major / xyz torque layout
        @test ds.y_T == reduce(vcat, [vec(d.torques) for d in data])
    end

    @testset "adiabatic-moment trio: moments_bare / constraint_axes / constraint_mode" begin
        n = 3
        moments = [1.2 0.0 0.0; 0.0 1.2 0.0; 0.0 0.0 1.2]
        B = zeros(3, n); B[1, 1] = 0.05
        M = [0.1 -0.2 0.0; 0.0 0.3 -0.1; 1.1 -1.0 0.02]     # signed, one near-zero: legal
        ax = repeat([0.0, 0.0, 1.0], 1, n)
        mk(; kw...) = SpinDatum(-1.0, moments, B; kw...)

        # defaults: absent, and the datum is what it was before the extension
        d0 = mk()
        @test d0.moments_bare === nothing
        @test d0.constraint_axes === nothing
        @test d0.constraint_mode === nothing
        # the five-field direct form is the same datum
        d5 = SpinDatum(d0.energy, d0.directions, d0.magmoms, d0.field, d0.torques)
        @test d5.moments_bare === nothing && d5.constraint_mode === nothing
        # ... and it converts, as the default constructor it replaced did
        d5i = SpinDatum(0, d0.directions, [1, 1, 1], d0.field, d0.torques)
        @test d5i.energy === 0.0 && d5i.magmoms == [1.0, 1.0, 1.0]

        # moments_bare: finiteness is the ONLY value constraint (signed, zero-crossing
        # magnitudes are the point of the vector storage)
        d = mk(; moments_bare = M)
        @test d.moments_bare == M
        Mbad = copy(M); Mbad[2, 2] = Inf
        @test_throws ArgumentError mk(; moments_bare = Mbad)
        @test_throws ArgumentError mk(; moments_bare = zeros(3, n + 1))

        # constraint_mode: only the two physical classes exist
        for bad in (0, 2, 3, 5, -1)
            @test_throws ArgumentError mk(; constraint_mode = bad)
        end
        d4 = mk(; moments_bare = M, constraint_mode = 4)
        @test d4.constraint_mode == 4 && d4.constraint_axes === nothing

        # mode 1 (transverse-penalty type) requires the axes; axes without a declared
        # mode are refused (the axis rule is keyed by mode, never by field presence)
        @test_throws ArgumentError mk(; moments_bare = M, constraint_mode = 1)
        @test_throws ArgumentError mk(; constraint_axes = ax)
        d1 = mk(; moments_bare = M, constraint_axes = ax, constraint_mode = 1)
        @test d1.constraint_mode == 1 && d1.constraint_axes == ax
        # an all-zero axes matrix satisfies "axes present" in letter only: refused
        @test_throws ArgumentError mk(; moments_bare = M, constraint_axes = zeros(3, n),
                                      constraint_mode = 1)
        # a Bool is not a class
        @test_throws ArgumentError mk(; constraint_mode = true)
        # mode 1 without the bare moments is legal (axis without target)
        @test mk(; constraint_axes = ax, constraint_mode = 1).moments_bare === nothing

        # axes columns: unit vectors, or exactly zero (= "no axis for this atom");
        # near-zero noise and off-unit columns are refused, not normalized
        axz = copy(ax); axz[:, 2] .= 0.0
        dz = mk(; constraint_axes = axz, constraint_mode = 1)
        @test dz.constraint_axes[:, 2] == zeros(3)
        axbad = copy(ax); axbad[:, 2] .= 0.5
        @test_throws ArgumentError mk(; constraint_axes = axbad, constraint_mode = 1)
        axeps = copy(ax); axeps[:, 2] .= 1e-9
        @test_throws ArgumentError mk(; constraint_axes = axeps, constraint_mode = 1)
        axnan = copy(ax); axnan[1, 1] = NaN
        @test_throws ArgumentError mk(; constraint_axes = axnan, constraint_mode = 1)
        axoff = copy(ax); axoff[3, 1] = 1.0 + 2e-6       # outside the 1e-6 band
        @test_throws ArgumentError mk(; constraint_axes = axoff, constraint_mode = 1)
        # inside the norm band but with a component > 1: refused — the component
        # bound is the load-bearing half (the harmonic kernels' |z| ≤ 1 domain),
        # exactly as `_validate_config` / `Harmonics._validate_unit` refuse it
        axpole = copy(ax); axpole[3, 1] = 1.0 + 5e-7
        @test_throws ArgumentError mk(; constraint_axes = axpole, constraint_mode = 4)
        # inside the band AND inside the component bound: accepted
        axin = copy(ax); axin[:, 1] = [1.0, 1.0, 0.0] ./ sqrt(2) .* (1 + 5e-7)
        @test mk(; constraint_axes = axin, constraint_mode = 4).constraint_axes == axin

        # the two construction paths of one type carry the fields identically, and
        # the derived E/T fields are bit-identical with and without the trio
        dd = SpinDatum(d1.energy, d1.directions, d1.magmoms, d1.field, d1.torques,
                       M, ax, 1)
        for f in (:moments_bare, :constraint_axes, :constraint_mode)
            @test getfield(dd, f) == getfield(d1, f)
        end
        for f in (:energy, :directions, :magmoms, :field, :torques)
            @test getfield(d1, f) == getfield(d0, f)
        end

        # the E and T paths are inert to the new fields: identical design matrices
        # and targets through SCEDataset
        lat = Lattice(Matrix(3.0 * I(3)))
        cr = Crystal(lat, [0.2 -0.2; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
        basis = SCEBasis(cr, BasisSpec(; nbody = 2, cutoff = 1.5, lmax = [1],
                                       isotropy = true))
        rng = MersenneTwister(17)
        raw = [(0.1 * k, randn(rng, 3, 2), 0.1 .* randn(rng, 3, 2)) for k = 1:4]
        plain = [SpinDatum(e, m, b) for (e, m, b) in raw]
        withm = [SpinDatum(e, m, b; moments_bare = 1.5 .* m, constraint_mode = 4)
                 for (e, m, b) in raw]
        dsp = SCEDataset(basis, plain)
        dsm = SCEDataset(basis, withm)
        @test dsp.X_E == dsm.X_E && dsp.y_E == dsm.y_E
        @test dsp.X_T == dsm.X_T && dsp.y_T == dsm.y_T
    end
end
