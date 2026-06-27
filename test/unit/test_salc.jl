using Test
using MagestyRebuild
using MagestyRebuild: _assemble_spacegroup, evaluate
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

rand_config(rng, nat) = reduce(hcat, (v = randn(rng, 3); SVector{3,Float64}(v / norm(v)) for _ = 1:nat))

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
    end

    @testset "every SALC is space-group invariant: Φ(g·e) = Φ(e)" begin
        rng = MersenneTwister(7)
        for _ = 1:20
            e = Matrix(rand_config(rng, 2))
            for g = 1:length(sg.ops)
                ge = act_on_config(sg, g, e)
                for s in basis.salcs
                    @test isapprox(evaluate(s, ge), evaluate(s, e); atol = 1e-9, rtol = 1e-8)
                end
            end
        end
    end

    @testset "every SALC is time-reversal even: Φ(-e) = Φ(e)" begin
        rng = MersenneTwister(11)
        for _ = 1:20
            e = Matrix(rand_config(rng, 2))
            for s in basis.salcs
                @test isapprox(evaluate(s, -e), evaluate(s, e); atol = 1e-9, rtol = 1e-8)
            end
        end
    end

    @testset "no odd-Σl channels survive (time reversal)" begin
        for s in basis.salcs
            @test iseven(sum(s.ls))
        end
    end

    @testset "isotropy keeps only Lf = 0" begin
        iso = build_salc_basis(crystal, sg, clusters; lmax_by_species = [2], isotropy = true)
        @test all(s -> s.Lf == 0, iso.salcs)
        @test length(iso) > 0
    end
end
