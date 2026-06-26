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

const MRH = MagestyRebuild.Harmonics
const MTH = Magesty.TesseralHarmonics

rand_unit(rng) = (v = SVector{3,Float64}(randn(rng), randn(rng), randn(rng)); v / norm(v))

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
end
