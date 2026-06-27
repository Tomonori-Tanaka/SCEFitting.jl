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
using OffsetArrays
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

    @testset "cluster orbits with Spglib (M6)" begin
        Lattice = MagestyRebuild.Lattice
        Crystal = MagestyRebuild.Crystal
        SpglibBackend = MagestyRebuild.SpglibBackend
        analyze = MagestyRebuild.analyze_symmetry
        a = 3.0

        # simple cubic (1 atom): nearest neighbors are the 6 edge atoms, all
        # symmetry-equivalent → one 1-body orbit and one 2-body orbit.
        sc = Crystal(Lattice(Matrix(a * I(3))), reshape([0.0, 0.0, 0.0], 3, 1), [1], ["X"])
        sg = analyze(SpglibBackend(), sc)
        cs = MagestyRebuild.build_clusters(sc, MagestyRebuild.build_neighbor_list(sc, a + 0.1),
                                           sg; nbody = 2)
        @test length(cs.by_body[1]) == 1
        @test length(cs.by_body[2]) == 1

        # bcc conventional cell (2 atoms): the 8 nearest neighbors are equivalent,
        # and corner/center sites are equivalent under body-centering.
        bcc = Crystal(Lattice(Matrix(a * I(3))), [0.0 0.5; 0.0 0.5; 0.0 0.5], [1, 1], ["X"])
        sgb = analyze(SpglibBackend(), bcc)
        csb = MagestyRebuild.build_clusters(bcc, MagestyRebuild.build_neighbor_list(bcc, 2.7),
                                            sgb; nbody = 2)
        @test length(csb.by_body[1]) == 1
        @test length(csb.by_body[2]) == 1
    end

    @testset "SALC basis on a real chain (M7)" begin
        Lattice = MagestyRebuild.Lattice
        Crystal = MagestyRebuild.Crystal
        SpglibBackend = MagestyRebuild.SpglibBackend
        analyze = MagestyRebuild.analyze_symmetry
        evaluate = MagestyRebuild.evaluate

        # 4-atom chain along z (so neighboring spins differ)
        lat = Lattice(Matrix(Diagonal([8.0, 8.0, 10.0])))
        frac = [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.25 0.5 0.75]
        chain = Crystal(lat, frac, [1, 1, 1, 1], ["Fe"])
        sg = analyze(SpglibBackend(), chain)
        nl = MagestyRebuild.build_neighbor_list(chain, 2.6)   # nearest neighbors only
        cs = MagestyRebuild.build_clusters(chain, nl, sg; nbody = 2)
        # full basis including anisotropic (Lf>0) channels
        basis = MagestyRebuild.build_salc_basis(chain, sg, cs; lmax_by_species = [2])
        @test any(s -> s.Lf > 0, basis.salcs)   # anisotropic channels present

        rng = MersenneTwister(5)
        randcfg() = reduce(hcat, [SVector{3,Float64}((v = randn(rng, 3); v / norm(v))) for _ = 1:4])

        # invariance under the real space group (all Lf, non-collinear spins)
        for _ = 1:5
            e = randcfg()
            for g = 1:length(sg.ops)
                R = Matrix(sg.ops[g].rotation_cart)
                ge = similar(e)
                for a = 1:4
                    b = findfirst(==(a), @view sg.map_sym[:, g])
                    ge[:, a] = R * e[:, b]
                end
                for s in basis.salcs
                    @test isapprox(evaluate(s, ge), evaluate(s, e); atol = 1e-8, rtol = 1e-7)
                end
            end
        end

        # the isotropic (1,1)/Lf=0 SALC is the Heisenberg invariant: Φ = √3 Σ_members e_i·e_j
        heis = only(filter(s -> s.ls == [1, 1] && s.Lf == 0, basis.salcs))
        for _ = 1:10
            e = randcfg()
            ref = sqrt(3.0) * sum(dot(e[:, m.atoms[1]], e[:, m.atoms[2]]) for m in heis.members)
            @test isapprox(evaluate(heis, e), ref; atol = 1e-9, rtol = 1e-8)
        end
    end

    @testset "Heisenberg chain end-to-end: recover J (M8/M9)" begin
        Lattice = MagestyRebuild.Lattice
        Crystal = MagestyRebuild.Crystal
        SpglibBackend = MagestyRebuild.SpglibBackend

        lat = Lattice(Matrix(Diagonal([8.0, 8.0, 10.0])))
        frac = [zeros(2, 4); [0.0 0.25 0.5 0.75]]
        chain = Crystal(lat, frac, [1, 1, 1, 1], ["Fe"])
        interaction = MagestyRebuild.Interaction(; nbody = 2, pair_cutoff = 2.6,
                                                 lmax = [1], isotropy = true)
        basis = MagestyRebuild.SCEBasis(chain, interaction; backend = SpglibBackend())
        @test length(basis.salcs) == 1          # only the nearest-neighbor Heisenberg SALC
        heis = basis.salcs.salcs[1]

        rng = MersenneTwister(3)
        cfg() = (M = Matrix{Float64}(undef, 3, 4);
                 for a = 1:4; v = randn(rng, 3); M[:, a] = v / norm(v); end; M)
        J_true = 0.0137
        configs = [cfg() for _ = 1:30]
        # physical Heisenberg energies E = J Σ_{undirected nn} e_i·e_j
        E = [J_true * 0.5 * sum(dot(c[:, m.atoms[1]], c[:, m.atoms[2]]) for m in heis.members)
             for c in configs]

        ds = MagestyRebuild.SCEDataset(basis, configs, E)
        f = MagestyRebuild.fit(MagestyRebuild.SCEFit, ds, MagestyRebuild.OLS())
        @test MagestyRebuild.r2_energy(f) ≈ 1.0 atol = 1e-10
        @test isapprox(MagestyRebuild.predict_energy(f, configs), E; atol = 1e-10)
        # recover J from the SALC coefficient: Φ = 2√3·Σ_{undirected} e_i·e_j
        J_recovered = 2 * sqrt(3.0) * MagestyRebuild.coef(f)[1]
        @test isapprox(J_recovered, J_true; rtol = 1e-8)

        # torque from the fitted model vs the Heisenberg closed form. With
        # E = J·0.5·Σ_members e_i·e_j, ∂E/∂e_a = J·0.5·Σ_members(δ_{a,i} e_j + δ_{a,j} e_i)
        # and τ_a = e_a × ∂E/∂e_a — convention-independent ground truth.
        function analytic_torque(c, J)
            nat = size(c, 2)
            G = zeros(3, nat)
            for mem in heis.members
                i, j = mem.atoms[1], mem.atoms[2]
                G[:, i] .+= (J * 0.5) .* c[:, j]
                G[:, j] .+= (J * 0.5) .* c[:, i]
            end
            return reduce(hcat, cross(SVector{3}(c[:, a]), SVector{3}(G[:, a])) for a = 1:nat)
        end
        for c in configs[1:5]
            @test isapprox(MagestyRebuild.predict_torque(f, c), analytic_torque(c, J_true);
                           atol = 1e-9)
        end
    end

    @testset "N-body SALC dimensions vs Magesty (kagome, 1/2/3-body)" begin
        # A kagome cell (P6/mmm): three equivalent sites whose triangle has C₃ᵥ site
        # symmetry permuting them — exercises coupling-path mixing and, for unequal-l
        # multisets like (1,1,2), l-ordering mixing. The number of independent
        # invariants per (body, ls, Lf) channel is convention-free, so the two
        # from-scratch constructions must agree exactly. This validates *both*
        # implementations (and would expose a bug in either) at N ≥ 3.
        a, c = 2.0, 8.0
        A = [a -a/2 0.0; 0.0 a*sqrt(3)/2 0.0; 0.0 0.0 c]
        frac = [0.0 0.5 0.5; 0.5 0.0 0.5; 0.0 0.0 0.0]
        kd = [1, 1, 1]
        cutoff, nbody = 1.2, 3

        # Magesty: count SALCs per channel by key-group (= design-matrix columns),
        # filtered to per-site l ≤ 2 to match the rebuild's subset.
        st = Magesty.Structures.Structure(A, SVector{3,Bool}(true, true, true), ["Fe"], kd,
                                          frac; verbosity = false)
        sym = Magesty.Symmetries.Symmetry(st, 1e-5; verbosity = false)
        cutarr = OffsetArray(fill(cutoff, nbody - 1, 1, 1), 2:nbody, 1:1, 1:1)
        clu = Magesty.Clusters.Cluster(st, sym, nbody, cutarr; verbosity = false)
        sb = Magesty.SALCBases.SALCBasis(st, sym, clu, [2], OffsetArray([4, 6], 2:nbody),
                                         nbody; isotropy = false, verbosity = false)
        magcount = Dict{Tuple{Int,Vector{Int},Int},Int}()
        for group in sb.salc_list
            cbc = group[1]
            all(l -> l <= 2, cbc.ls) || continue
            k = (length(cbc.atoms), sort(cbc.ls), cbc.Lf)
            magcount[k] = get(magcount, k, 0) + 1
        end

        # MagestyRebuild: count SALCs per channel.
        Lattice = MagestyRebuild.Lattice
        Crystal = MagestyRebuild.Crystal
        xtal = Crystal(Lattice(A), frac, kd, ["Fe"])
        sg = MagestyRebuild.analyze_symmetry(MagestyRebuild.SpglibBackend(), xtal)
        cs = MagestyRebuild.build_clusters(xtal, MagestyRebuild.build_neighbor_list(xtal, cutoff),
                                           sg; nbody = nbody)
        basis = MagestyRebuild.build_salc_basis(xtal, sg, cs; lmax_by_species = [2])
        mrcount = Dict{Tuple{Int,Vector{Int},Int},Int}()
        for s in basis.salcs
            k = (s.key.body, sort(s.key.ls), s.key.Lf)
            mrcount[k] = get(mrcount, k, 0) + 1
        end

        @test sym.spacegroup_number == 191
        @test mrcount == magcount                       # every (body, ls, Lf) dimension agrees
        @test sum(values(mrcount)) == 42
        @test mrcount[(3, [1, 1, 2], 2)] == 5           # a multi-ordering 3-body channel

        # Independent ground truth on *Magesty*: its own 3-body SALCs are invariant.
        magsc(e) = Magesty.SpinConfigs.SpinConfig(0.0, ones(size(e, 2)), e, zeros(3, size(e, 2)))
        function act_mag(e, g)
            nat = size(e, 2)
            R = Matrix(sym.symdata[g].rotation_cart)
            out = similar(e)
            for x = 1:nat
                b = findfirst(==(x), @view sym.map_sym[:, g])
                out[:, x] = R * e[:, b]
            end
            return out
        end
        rng = MersenneTwister(42)
        rc() = reduce(hcat, [SVector{3,Float64}((v = randn(rng, 3); v / norm(v))) for _ = 1:3])
        for _ = 1:4
            e = Matrix(rc())
            Xe = Magesty.Fitting.build_design_matrix_energy(sb.salc_list, [magsc(e)], sym;
                                                           verbosity = false)
            for g = 1:sym.nsym
                Xg = Magesty.Fitting.build_design_matrix_energy(sb.salc_list,
                                                               [magsc(act_mag(e, g))], sym;
                                                               verbosity = false)
                @test isapprox(Xg, Xe; atol = 1e-7)
            end
        end
    end

    # VASP I/O: the from-scratch POSCAR / OSZICAR readers agree bit-for-bit with
    # Magesty's proven parsers (synthetic files, since real VASP outputs are not
    # vendored). Default SAXIS = ẑ ⇒ no rotation, matching Magesty's raw parse.
    @testset "VASP POSCAR / OSZICAR parsing vs Magesty" begin
        dir = mktempdir()

        poscar = joinpath(dir, "POSCAR")
        write(poscar,
            "FePt\n1.0\n 3.1 0.2 0.0\n 0.0 3.0 0.0\n 0.0 0.0 4.0\nFe Pt\n1 1\nDirect\n 0.1 0.2 0.3\n 0.5 0.5 0.5\n")
        mr = MagestyRebuild.VASP.read_poscar(poscar)
        mg = Magesty.VaspIO.parse_poscar(poscar)
        for i = 1:3
            @test mr.lattice.vectors[:, i] ≈ mg.lattice_vectors[i]   # column i = lattice vector i
        end
        for a = 1:size(mr.frac_positions, 2)
            @test mr.frac_positions[:, a] ≈ mg.positions[a]
        end
        kd = reduce(vcat, [fill(i, n) for (i, n) in enumerate(mg.numbers)])
        @test mr.species == kd
        @test mr.species_labels == mg.element_symbols

        osz = joinpath(dir, "OSZICAR")
        write(osz,
            " ion   MW_int   M_int\n   1  1.0 0.2 0.0  1.1 0.2 0.0\n   2  0.0 0.0 2.0  0.0 0.0 2.1\n" *
            " lambda*MW_perp\n   1  0.0 0.02 0.0\n   2  0.03 0.0 0.0\n" *
            "   1 F= -.84314080E+02 E0= -.84200000E+02  d E = 0\n")
        mrd = MagestyRebuild.read_configs(MagestyRebuild.VASP.Oszicar(osz))[1]
        mgd = Magesty.VaspIO.parse_oszicar(osz; energy_kind = "f", mint = false)
        @test mrd.energy ≈ mgd.energy
        for i = 1:length(mrd.magmoms)
            @test mrd.magmoms[i] .* mrd.directions[:, i] ≈ mgd.magmom[i, :]   # moment vector
            @test mrd.field[:, i] ≈ mgd.magfield[i, :]
        end
    end
end
