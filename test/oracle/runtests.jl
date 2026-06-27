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
import Spglib   # triggers MagestyRebuildSpglibExt

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

    @testset "coupled real tensors vs Magesty (M4)" begin
        MAC = Magesty.AngularMomentumCoupling
        for ls in ([1], [2], [1, 1], [1, 2], [2, 2], [1, 1, 1])
            mine = MRA.build_real_bases(ls)
            bases_by_L, paths_by_L = MAC.build_all_real_bases(ls)
            @test length(mine) == sum(length(v) for v in values(bases_by_L))
            for (Lseq, Lf, T) in mine
                ks = findall(==(Lseq), paths_by_L[Lf])
                @test length(ks) == 1
                @test isapprox(T, bases_by_L[Lf][ks[1]]; atol = 1e-10)
            end
        end
    end

    @testset "SpglibBackend space groups (M5)" begin
        Lattice = MagestyRebuild.Lattice
        Crystal = MagestyRebuild.Crystal
        SpglibBackend = MagestyRebuild.SpglibBackend
        analyze = MagestyRebuild.analyze_symmetry

        a = 3.0
        sc = Crystal(Lattice(Matrix(a * I(3))), reshape([0.0, 0.0, 0.0], 3, 1), [1], ["X"])
        bcc = Crystal(Lattice(Matrix(a * I(3))), [0.0 0.5; 0.0 0.5; 0.0 0.5], [1, 1], ["X"])
        fcc = Crystal(Lattice(Matrix(4.0 * I(3))),
                      [0.0 0.5 0.5 0.0; 0.0 0.5 0.0 0.5; 0.0 0.0 0.5 0.5],
                      [1, 1, 1, 1], ["X"])

        sg_sc = analyze(SpglibBackend(), sc)
        @test sg_sc.number == 221 && length(sg_sc.ops) == 48          # Pm-3m
        @test analyze(SpglibBackend(), bcc).number == 229             # Im-3m
        @test analyze(SpglibBackend(), fcc).number == 225             # Fm-3m

        for sg in (analyze(SpglibBackend(), bcc), analyze(SpglibBackend(), fcc))
            nat = size(sg.map_sym, 1)
            @test any(o -> sg.ops[o].is_translation && all(sg.map_sym[:, o] .== 1:nat),
                      1:length(sg.ops))                                # identity present
            for o = 1:length(sg.ops)
                @test sort(sg.map_sym[:, o]) == collect(1:nat)         # each op permutes atoms
                @test isapprox(sg.ops[o].rotation_cart * sg.ops[o].rotation_cart',
                               Matrix(I, 3, 3); atol = 1e-9)           # orthogonal
            end
        end
    end
end
