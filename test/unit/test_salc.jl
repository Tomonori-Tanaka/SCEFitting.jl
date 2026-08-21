using Test
using SCEFitting
using SCEFitting: _assemble_spacegroup, evaluate_salc, _canonicalize_members
using StaticArrays
using LinearAlgebra
using Random

# Apply a space-group op to a spin config: (g·e)_a = R_g · e_{preimage(a)}.
function act_on_config(sg, g, e)
    nat = size(e, 2)
    R = Matrix(sg.ops[g].rotation_cart)
    out = similar(e)
    for a = 1:nat
        b = findfirst(==(a), @view sg.map_sym[:, g])   # preimage of a under g
        out[:, a] = R * e[:, b]
    end
    return out
end

@testset "SALC" begin
    lat = Lattice(Matrix(3.0 * I(3)))
    crystal = Crystal(lat, [0.2 -0.2; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
    # manual Z₂ space group: identity + inversion (swaps the centrosymmetric pair)
    rots = [SMatrix{3,3,Float64}(I), SMatrix{3,3,Float64}(-I)]
    trans = [SVector{3,Float64}(0, 0, 0), SVector{3,Float64}(0, 0, 0)]
    sg = _assemble_spacegroup(crystal, rots, trans, "manual", 0; tol = 1e-6)
    nl = build_neighbor_list(crystal, 1.5)
    clusters = build_clusters(crystal, nl, sg; nbody = 2)
    basis = build_salc_basis(crystal, sg, clusters; lmax_by_species = [2])

    @testset "basis is nonempty and keys are sorted/unique" begin
        @test length(basis) > 0
        @test issorted(basis.keys)
        @test allunique(basis.keys)
        @test basis.fingerprint == hash(basis.keys)
        @test any(s -> s.Lf > 0, basis.salcs)   # anisotropic channels present & tested
    end

    @testset "every SALC is space-group invariant: Φ(g·e) = Φ(e)" begin
        rng = MersenneTwister(7)
        for _ = 1:20
            e = Matrix(randcfg(rng, 2))
            for g = 1:length(sg.ops)
                ge = act_on_config(sg, g, e)
                for s in basis.salcs
                    @test isapprox(evaluate_salc(s, ge), evaluate_salc(s, e); atol = 1e-9, rtol = 1e-8)
                end
            end
        end
    end

    @testset "every SALC is time-reversal even: Φ(-e) = Φ(e)" begin
        rng = MersenneTwister(11)
        for _ = 1:20
            e = Matrix(randcfg(rng, 2))
            for s in basis.salcs
                @test isapprox(evaluate_salc(s, -e), evaluate_salc(s, e); atol = 1e-9, rtol = 1e-8)
            end
        end
    end

    @testset "no odd-Σl channels survive (time reversal)" begin
        for s in basis.salcs
            @test iseven(sum(spin_ls(s.key)))
        end
    end

    @testset "members are canonical and fold exactly" begin
        perms = Dict(1 => ([1], [1]), 2 => ([2, 1], [1, 2]))
        for s in basis.salcs
            check_canonical_members(s)
            @test same_members(_canonicalize_members(s.members), s.members)  # idempotent
            split_roundtrip_exact(s, perms)
        end
    end

    @testset "isotropy keeps only Lf = 0" begin
        iso = build_salc_basis(crystal, sg, clusters; lmax_by_species = [2], isotropy = true)
        @test all(s -> s.Lf == 0, iso.salcs)
        @test length(iso) > 0
    end

    # Function-space reduction: a supercell that folds an orbit's distinct cluster
    # instances onto one atom set aggregates their tensors, and SALCs whose aggregate
    # vanishes (or goes linearly dependent) must be dropped — loudly — or the design
    # matrix is silently rank deficient (bulk MnTe + SOC: 51 emitted, rank 37,
    # max|coef| ~1e7 from OLS with normal-looking R²; found 2026-08-12).
    #
    # Closed-form fixture: a CsCl-type cubic cell — A at (0,0,0), B at (½,½,½) —
    # under the full 48-op cubic point group. The A–B pair sits exactly at the WS-cell
    # corner (√3/2·L), so all 8 body-diagonal images are kept as one tie shell of ONE
    # orbit, all on the same atom pair (1,2). With lmax = 1 the pair channel is
    # bilinear, and the hand oracle is classical invariant theory, independent of the
    # implementation:
    #   • the stabilizer of one image is C₃ᵥ along its diagonal (order 6). Its A₁
    #     (invariant) content per Lf: Lf = 0 → 1 (e₁·e₂); Lf = 1 → 0 (of the three
    #     components of e₁×e₂, the axial one n̂·(e₁×e₂) is A₂ — odd under the
    #     vertical mirrors — and the two perpendicular ones span E, dying under
    #     C₃); Lf = 2 → 1 (the axial quadrupole e₁·(n̂n̂ᵀ−I/3)·e₂). So the builder
    #     emits 2 SALCs;
    #   • aggregated over the 8 diagonals n̂ₖ = (±1,±1,±1)/√3:
    #       Lf = 2 :  Σₖ n̂ₖn̂ₖᵀ = (8/3)·I  ⇒ its traceless part sums to zero → dropped
    #       Lf = 0 :  e₁·e₂                                               → survives
    #   so exactly ONE basis function must remain, and it must be a multiple of
    #   e₁·e₂ (the only cubic-invariant bilinear).
    @testset "tie-shell aggregation: redundant SALCs are dropped, loudly" begin
        L = 2.0
        cubic = Crystal(Lattice(Matrix(L * I(3))),
                        [0.0 0.5; 0.0 0.5; 0.0 0.5], [1, 2], ["A", "B"])
        # all 48 signed permutation matrices (the cubic point group, exact)
        crots = SMatrix{3,3,Float64}[]
        for p in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)),
            s1 in (1, -1), s2 in (1, -1), s3 in (1, -1)

            W = zeros(3, 3)
            W[1, p[1]] = s1
            W[2, p[2]] = s2
            W[3, p[3]] = s3
            push!(crots, SMatrix{3,3,Float64}(W))
        end
        ctrans = [SVector{3,Float64}(0, 0, 0) for _ in crots]
        csg = _assemble_spacegroup(cubic, crots, ctrans, "Pm-3m(manual)", 221;
                                   tol = 1e-6)
        cnl = build_neighbor_list(cubic, Inf, MinimumImage())
        ccl = build_clusters(cubic, cnl, csg; nbody = 2)
        # one pair orbit, and its members alias onto the single atom pair (1, 2)
        @test length(ccl.by_body[2]) == 1
        cb = @test_logs (:warn, r"SALC reduction: dropped 1 redundant.*Lf = 2, block 1: zero"s) match_mode = :any begin
            build_salc_basis(cubic, csg, ccl; lmax_by_species = [1, 1])
        end
        @test length(cb) == 1
        @test cb.salcs[1].Lf == 0
        # the survivor IS the Heisenberg invariant: Φ(e) / (e₁·e₂) is one constant
        rng = MersenneTwister(42)
        ratios = map(1:8) do _
            e = randn(rng, 3, 2)
            e ./= reshape(map(norm, eachcol(e)), 1, :)
            evaluate_salc(cb.salcs[1], e) / dot(e[:, 1], e[:, 2])
        end
        @test all(r -> isapprox(r, ratios[1]; rtol = 1e-12), ratios)
        @test abs(ratios[1]) > 1e-6            # a nonzero multiple, not 0/0 noise
    end

    @testset "build is deterministic / thread-safe" begin
        # `build_salc_basis` processes orbits in parallel (`Threads.@threads`), writing
        # disjoint per-orbit results then sorting by key, so the output is byte-for-byte
        # independent of thread count and of the parallel completion order. A rebuild
        # must reproduce keys *and* the folded tensors exactly; run the suite with
        # `julia -t N>1` to exercise the threaded path (a shared-state race fails here).
        Threads.nthreads() == 1 &&
            @warn "determinism test runs serial; launch `julia -t N>1` to exercise the threaded path"
        b2 = build_salc_basis(crystal, sg, clusters; lmax_by_species = [2])
        @test b2.keys == basis.keys
        @test b2.fingerprint == basis.fingerprint
        @test length(b2) == length(basis)
        for (s1, s2) in zip(basis.salcs, b2.salcs)
            @test s1.key == s2.key
            @test length(s1.members) == length(s2.members)
            for (m1, m2) in zip(s1.members, s2.members)
                for (t1, t2) in zip(m1.terms, m2.terms)
                    @test t1.slots == t2.slots
                    @test t1.folded == t2.folded     # exact equality, not isapprox
                end
            end
        end
    end
end
