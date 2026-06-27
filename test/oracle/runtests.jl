# Oracle validation: MagestyRebuild's from-scratch numerics vs Magesty.jl.
#
# Conventions are pinned independently by the core suite (closed forms, on-sphere
# finite differences). Here we additionally confirm bit-for-bit agreement with
# Magesty, which guards the convention all the way to the design matrix and
# predictions downstream. Gauge-dependent quantities (raw SALC coefficients) are
# NOT compared directly — only convention-fixed kernels and gauge-invariant
# aggregates / predictions (added as later milestones land).

using Test
using Random
using StaticArrays
using LinearAlgebra
import MagestyRebuild
import Magesty
import WignerSymbols

const MRH = MagestyRebuild.Harmonics
const MRA = MagestyRebuild.AngularMomentum
const MTH = Magesty.TesseralHarmonics
const MTR = Magesty.RotationMatrix

rand_unit(rng) = (v = SVector{3,Float64}(randn(rng), randn(rng), randn(rng)); v / norm(v))

function rand_rotation(rng)
    a = SVector{3,Float64}(randn(rng), randn(rng), randn(rng))
    a = a / norm(a)
    θ = 2π * rand(rng)
    Kc = @SMatrix [0.0 -a[3] a[2]; a[3] 0.0 -a[1]; -a[2] a[1] 0.0]
    return SMatrix{3,3,Float64}(I) + sin(θ) * Kc + (1 - cos(θ)) * (Kc * Kc)
end

@testset "oracle vs Magesty" begin
    @testset "tesseral harmonics Zₗₘ and ∇Zₗₘ (M2)" begin
        rng = MersenneTwister(11)
        lmax = 4
        for _ = 1:200
            u = rand_unit(rng)
            for l = 0:lmax, m = -l:l
                @test MRH.Zlm(l, m, u) ≈ MTH.Zₗₘ(l, m, u) atol = 1e-13 rtol = 1e-12
                g_new = MRH.grad_Zlm(l, m, u)
                g_ref = SVector{3,Float64}(MTH.∂ᵢZlm(l, m, u))
                @test isapprox(g_new, g_ref; atol = 1e-12, rtol = 1e-11)
            end
        end
    end

    @testset "Clebsch-Gordan vs WignerSymbols (M3)" begin
        for j1 = 0:3, j2 = 0:3
            for J = abs(j1 - j2):(j1+j2), m1 = -j1:j1, m2 = -j2:j2
                M = m1 + m2
                abs(M) <= J || continue
                ref = Float64(WignerSymbols.clebschgordan(j1, m1, j2, m2, J, M))
                @test MRA.clebsch_gordan(j1, m1, j2, m2, J, M) ≈ ref atol = 1e-12
            end
        end
    end

    @testset "real Wigner-D vs Magesty Δl (M3)" begin
        rng = MersenneTwister(31)
        for _ = 1:15
            R = rand_rotation(rng)
            α, β, γ = MTR.rotmat2euler(Matrix(R))
            for l = 1:4
                mine = MRA.wignerD_real(l, R)
                ref = MTR.Δl(l, α, β, γ)
                # Both are the rotation's representation in the same tesseral basis;
                # they agree up to the active/passive (transpose) convention.
                @test isapprox(mine, ref; atol = 1e-9) ||
                      isapprox(mine, permutedims(ref); atol = 1e-9)
            end
        end
    end
end
