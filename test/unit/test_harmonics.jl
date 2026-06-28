using Test
using SCEFitting: Harmonics
using StaticArrays
using LinearAlgebra
using Random

const Z = Harmonics.Zlm
const gZ = Harmonics.grad_Zlm

rand_unit(rng) = (v = SVector{3,Float64}(randn(rng), randn(rng), randn(rng)); v / norm(v))

@testset "harmonics" begin
    rng = MersenneTwister(20260627)

    # Convention-independent anchor #1: closed-form standard real solid harmonics
    # on the unit sphere. These pin the absolute normalization and signs (the
    # gauge-invariant / oracle tests cannot — a self-consistent wrong sign would
    # pass those).
    @testset "closed-form low-l values" begin
        c1 = sqrt(3 / (4π))
        c2 = sqrt(15 / (4π))
        c20 = sqrt(5 / (16π))
        c22 = sqrt(15 / (16π))
        for _ = 1:50
            u = rand_unit(rng)
            x, y, z = u
            @test Z(0, 0, u) ≈ 1 / (2 * sqrt(π))
            @test Z(1, -1, u) ≈ c1 * y
            @test Z(1, 0, u) ≈ c1 * z
            @test Z(1, 1, u) ≈ c1 * x
            @test Z(2, -2, u) ≈ c2 * x * y
            @test Z(2, -1, u) ≈ c2 * y * z
            @test Z(2, 0, u) ≈ c20 * (3z^2 - 1)
            @test Z(2, 1, u) ≈ c2 * x * z
            @test Z(2, 2, u) ≈ c22 * (x^2 - y^2)
        end
    end

    # Convention-independent anchor #2: the gradient is tangent (u·∇Z = 0) and
    # agrees with an on-sphere central difference. The FD perturbation is
    # renormalized onto the sphere (the harmonics are not degree-0 homogeneous),
    # which is exactly what the torque self-consistency relies on.
    @testset "gradient tangency + on-sphere central difference" begin
        ε = 1e-6
        for l = 0:4, m = -l:l
            for _ = 1:8
                u = rand_unit(rng)
                g = gZ(l, m, u)
                @test abs(dot(g, u)) < 1e-7
                a = SVector{3,Float64}(randn(rng), randn(rng), randn(rng))
                t1 = a - dot(a, u) * u
                t1 = t1 / norm(t1)
                t2 = cross(u, t1)
                for t in (t1, t2)
                    up = u + ε * t
                    up = up / norm(up)
                    um = u - ε * t
                    um = um / norm(um)
                    fd = (Z(l, m, up) - Z(l, m, um)) / (2ε)
                    @test isapprox(fd, dot(g, t); atol = 1e-5, rtol = 1e-5)
                end
            end
        end
    end

    @testset "validation" begin
        u = SVector{3,Float64}(0.0, 0.0, 1.0)
        @test_throws ArgumentError Z(-1, 0, u)
        @test_throws ArgumentError Z(1, 2, u)
        @test_throws ArgumentError Z(1, 0, SVector{3,Float64}(1.0, 1.0, 1.0))
        @test_throws ArgumentError gZ(2, -3, u)
    end

    @testset "lm_index contiguous and unique" begin
        idxs = [Harmonics.lm_index(l, m) for l = 0:4 for m = -l:l]
        @test sort(idxs) == collect(1:Harmonics.num_lm(4))
    end
end
