# Mixed-channel (decor) SALC engine. This package's models are pure spin; the
# decor engine exists because the pointed site-moment channel marks one site
# with a displacement decor. Gates at engine level: bitwise agreement with the
# production pure-spin engine (the anti-drift gate, including the two shapes
# whose enumeration order the upstream review found broken), invariant counts
# against a Cartesian projector that shares no code with the SALC machinery,
# the u = 0 degeneracy, mixed space-group invariance, and the refusals that
# keep a decorated SALC out of the pure-spin kernels.

using Test
using SCEFitting
using SCEFitting: SiteDecor, SiteFactor, Slot, SPIN, DISP, spin_decors, spin_ls,
                  _orbit_salcs_decors, _orbit_salcs, _build_wig_cache,
                  evaluate_salc, build_clusters, build_neighbor_list,
                  _assemble_spacegroup, salcs, LSUM_UNCAPPED
using LinearAlgebra
using StaticArrays
using Random

# shared helpers (same_members, ...) — included once by runtests.jl; standalone
# runs pull them in here
isdefined(@__MODULE__, :same_members) || include("testutils.jl")

# Molecule-in-a-box triangles (local copies of the test_nbody.jl fixtures, so
# this file runs standalone): Cs isosceles (ls = [1,1,2] splits into two
# ordering orbits sharing one sorted label) and C3v equilateral (one orbit of
# three assignments — the gauge is column-order dependent). These are the two
# shapes that expose an enumeration-order drift between the engines.
function _mx_triangle_cs(; L = 8.0)
    c = SVector{3,Float64}(0.5, 0.5, 0.5)
    offs = [SVector{3,Float64}(0.0, 1.5, 0.0), SVector{3,Float64}(-1.0, 0.0, 0.0),
            SVector{3,Float64}(1.0, 0.0, 0.0)]
    frac = reduce(hcat, [c + o / L for o in offs])
    crystal = Crystal(Lattice(Matrix(L * I(3))), frac, [1, 1, 1], ["Fe"])
    σ = SMatrix{3,3,Float64}([-1.0 0 0; 0 1 0; 0 0 1])
    rots = [SMatrix{3,3,Float64}(I), σ]
    trans = [(SMatrix{3,3,Float64}(I) - R) * c for R in rots]
    return crystal, _assemble_spacegroup(crystal, rots, trans, "Cs", 0; tol = 1e-6)
end

# Two-species version of the same triangle: the apex is species 2, the
# mirror-equivalent base pair species 1, so a per-species `lmax` of [2, 1] makes
# the sites carry DIFFERENT caps — the heterogeneous-lmax cross-engine shape.
function _mx_triangle_cs2(; L = 8.0)
    c = SVector{3,Float64}(0.5, 0.5, 0.5)
    offs = [SVector{3,Float64}(0.0, 1.5, 0.0), SVector{3,Float64}(-1.0, 0.0, 0.0),
            SVector{3,Float64}(1.0, 0.0, 0.0)]
    frac = reduce(hcat, [c + o / L for o in offs])
    crystal = Crystal(Lattice(Matrix(L * I(3))), frac, [2, 1, 1], ["Fe", "Co"])
    σ = SMatrix{3,3,Float64}([-1.0 0 0; 0 1 0; 0 0 1])
    trans = [(SMatrix{3,3,Float64}(I) - R) * c
             for R in [SMatrix{3,3,Float64}(I), σ]]
    return crystal, _assemble_spacegroup(crystal, [SMatrix{3,3,Float64}(I), σ],
                                         trans, "Cs", 0; tol = 1e-6)
end

function _mx_triangle_c3v(; a = 6.0, c = 6.0, r = 1.2)
    lat = Lattice(SMatrix{3,3,Float64}([a -a/2 0; 0 a*√3/2 0; 0 0 c]))
    ctr = SVector{3,Float64}(0.5, 0.5, 0.5)             # the C3 axis passes here
    ang = deg2rad.([90.0, 210.0, 330.0])
    offs = [SVector{3,Float64}(r * cos(t), r * sin(t), 0.0) for t in ang]
    frac = reduce(hcat, [ctr + lat.reciprocal * o for o in offs])
    crystal = Crystal(lat, frac, [1, 1, 1], ["Fe"])
    # a hexagonal cell, so C3 and the mirror are EXACTLY integral in fractional
    # coordinates — `_assemble_spacegroup` refuses a non-integer fractional rotation
    C3 = SMatrix{3,3,Float64}([0 -1 0; 1 -1 0; 0 0 1])
    M1 = SMatrix{3,3,Float64}([-1 1 0; 0 1 0; 0 0 1])
    rots = [SMatrix{3,3,Float64}(I), C3, C3 * C3, M1, C3 * M1, C3 * C3 * M1]
    trans = [(SMatrix{3,3,Float64}(I) - R) * ctr for R in rots]
    return crystal, _assemble_spacegroup(crystal, rots, trans, "C3v", 0; tol = 1e-6)
end

# All 48 signed permutation matrices = O_h in Cartesian (= fractional for sc).
function _oh48()
    mats = SMatrix{3,3,Float64,9}[]
    for p in [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
        for sx in (-1, 1), sy in (-1, 1), sz in (-1, 1)
            M = zeros(3, 3)
            signs = (sx, sy, sz)
            for i = 1:3
                M[i, p[i]] = signs[i]
            end
            push!(mats, SMatrix{3,3,Float64,9}(M))
        end
    end
    return mats
end

# ── independent invariant counters ────────────────────────────────────────────
#
# Both are plain Cartesian projectors: they average the group action over the
# multilinear function space written in the components of ê and u, and read the
# invariant count off the rank. They share no code with the SALC engine — no
# Clebsch–Gordan coefficients, no Wigner-D, no spherical harmonics — so they are
# an oracle for its counts rather than a restatement of them. They are legitimate
# because the rank-1 factors span exactly the Cartesian components
# (`R₁₋₁, R₁₀, R₁₁ = y, z, x`, and `Z₁ₘ(ê)` is a fixed linear bijection of `ê`)
# and the rank-2 factors span the traceless symmetric quadratic forms; with
# Σl_spin even the axial `det(R)` factors square to +1, so spin and displacement
# axes rotate the same way.

# Invariants of the degree-(1,1,1,1) forms Φ(e₁,e₂,u₁,u₂) under ops that carry a
# site swap: `swap(R) === true` means the op exchanges the two sites.
function _count_bond_1111_invariants(rots, swap)
    n = 81
    P = zeros(n, n)
    lin = LinearIndices((3, 3, 3, 3))
    for R in rots
        # coefficients of Φ(g·e, g·u) in terms of Φ's own coefficients
        M = zeros(n, n)
        for a = 1:3, b = 1:3, c = 1:3, d = 1:3
            for p = 1:3, q = 1:3, r = 1:3, s = 1:3
                w = R[a, p] * R[b, q] * R[c, r] * R[d, s]
                w == 0.0 && continue
                # e_{σ(1)}[p] e_{σ(2)}[q] u_{σ(1)}[r] u_{σ(2)}[s]: a swapping op
                # exchanges the roles of the two sites in BOTH channels
                out = swap(R) ? lin[q, p, s, r] : lin[p, q, r, s]
                M[out, lin[a, b, c, d]] += w
            end
        end
        P .+= M
    end
    P ./= length(rots)
    return rank(P; atol = 1e-8)
end

# Orthonormal (Frobenius) basis of the traceless symmetric 3×3 matrices — the
# l = 2 factor of a single site, whose group action is the congruence R·Q·Rᵀ.
const _Q5 = let
    b = [[0.0 1 0; 1 0 0; 0 0 0], [0.0 0 0; 0 0 1; 0 1 0], [0.0 0 1; 0 0 0; 1 0 0],
         [1.0 0 0; 0 -1 0; 0 0 0], [1.0 0 0; 0 1 0; 0 0 -2]]
    [Q ./ sqrt(sum(abs2, Q)) for Q in b]
end

_q5_action(R) = [sum((R * _Q5[j] * R') .* _Q5[i]) for i = 1:5, j = 1:5]

# Invariants of Φ(ê, u) built from `n` independent rank-2 factors at one site.
function _count_single_site_l2_invariants(rots, n::Int)
    dim = 5^n
    P = zeros(dim, dim)
    for R in rots
        D = _q5_action(R)
        M = n == 1 ? D : kron(D, D)
        P .+= M
    end
    P ./= length(rots)
    return rank(P; atol = 1e-8)
end

@testset "mixed-channel SALC engine" begin
    rng = MersenneTwister(0x51ce)

    # -- fixture A: 1-atom sc crystal with the full hand-assembled O_h group --
    # (a 1-atom cell hosts no MinimumImage pair — self-images are dropped — so
    # this crystal serves the single-site gates only)
    lat = Lattice(Matrix(2.0 * I(3)))
    xtal = Crystal(lat, zeros(3, 1), [1], ["Fe"])
    rots = _oh48()
    trs = [SVector{3,Float64}(0, 0, 0) for _ in rots]
    sg = _assemble_spacegroup(xtal, rots, trs, "Pm-3m", 221; tol = 1e-5)
    wc = _build_wig_cache(sg, 4)
    cs1 = build_clusters(xtal, build_neighbor_list(xtal, 2.1), sg; nbody = 1)
    O1 = cs1.by_body[1][1]

    # -- fixture B: two Fe at cart ±(0.5, 0, 0) in a 3.0 cube, D4h ops about the
    # origin (= the bond midpoint): the x-reversing half SWAPS the two sites, so
    # the bond stabilizer is the full 16-op d4h.
    latB = Lattice(Matrix(3.0 * I(3)))
    xtalB = Crystal(latB, [1 / 6 -1 / 6; 0.0 0.0; 0.0 0.0], [1, 1], ["Fe"])
    rotsB = [R for R in _oh48() if abs(R[1, 1]) == 1.0]     # D4h about x
    @test length(rotsB) == 16
    trsB = [SVector{3,Float64}(0, 0, 0) for _ in rotsB]
    sgB = _assemble_spacegroup(xtalB, rotsB, trsB, "P4/mmm", 123; tol = 1e-5)
    wcB = _build_wig_cache(sgB, 4)
    csB = build_clusters(xtalB, build_neighbor_list(xtalB, 1.1), sgB; nbody = 2)
    O2 = csB.by_body[2][1]                      # the 1.0-long x bond orbit

    # Per pure-spin label: the decor engine, handed that one label, must
    # reproduce the production engine's SALCs bitwise — same keys (block indices
    # included) and same members (slots and folded tensors, `==` on Float64).
    function anti_drift(xt, sgx, N, O, lmax, isotropy, wcx)
        old = _orbit_salcs(xt, sgx, N, 1, O, lmax, LSUM_UNCAPPED, isotropy, wcx)
        labels = unique(spin_ls(s.key) for s in old)
        @test !isempty(labels)
        # The decor engine sees a sorted multiset, so which site may carry which
        # rank has to be handed in: this is the per-species cap of `_enumerate_ls`,
        # written out here instead of called. That is not an independent
        # derivation — it is the same rule transcribed; get it wrong and the gate
        # below goes red, which is the safe direction.
        admit = t -> all(t[i].spin_l <= lmax[O.species[i]] for i in eachindex(t))
        for lab in labels
            olab = [s for s in old if spin_ls(s.key) == lab]
            new = _orbit_salcs_decors(xt, sgx, N, 1, O, [spin_decors(lab)],
                                      isotropy, wcx; admit = admit)
            @test [s.key for s in olab] == [s.key for s in new]
            @test all(s.key.L_S == s.key.Lf for s in new)
            for (a, b) in zip(olab, new)
                @test same_members(a.members, b.members)
            end
        end
        return length(old)
    end

    @testset "engines agree on pure spin (anti-drift gate)" begin
        @test anti_drift(xtal, sg, 1, O1, [4], false, wc) > 0
        @test anti_drift(xtalB, sgB, 2, O2, [2], false, wcB) > 0
        # The two review-blocker shapes: on the Cs triangle ls = [1,1,2] SPLITS
        # into two ordering orbits sharing one sorted label (block indices must
        # match the pure-spin emission order), and on the C3v triangle one orbit
        # carries three assignments (gauge order — bitwise match required, not
        # just span equality).
        for (xt, st) in (_mx_triangle_cs(), _mx_triangle_c3v())
            wc3 = _build_wig_cache(st, 2)
            cs3 = build_clusters(xt, build_neighbor_list(xt, 2.2), st; nbody = 3)
            @test anti_drift(xt, st, 3, cs3.by_body[3][1], [2], false, wc3) > 0
        end
    end

    @testset "heterogeneous lmax and the isotropy screen" begin
        # Per-species caps [2, 1]: the apex reaches l = 2, the base pair l = 1,
        # so a label's sites carry different ranks and the l = 1 sites are the
        # permuted pair.
        xt2, st2 = _mx_triangle_cs2()
        wc2 = _build_wig_cache(st2, 2)
        cs2 = build_clusters(xt2, build_neighbor_list(xt2, 2.2), st2; nbody = 3)
        O3 = cs2.by_body[3][1]
        @test anti_drift(xt2, st2, 3, O3, [2, 1], false, wc2) > 0
        # ...and the hook is load-bearing, not decoration: WITHOUT it the label
        # [1,1,2] also builds the ordering orbit that puts l = 2 on the l ≤ 1
        # apex, which shifts every `block` index of the admissible orbit. The
        # count below is the one the apex cap forbids.
        lab112 = spin_decors([1, 1, 2])
        capped = _orbit_salcs_decors(xt2, st2, 3, 1, O3, [lab112], false, wc2;
            admit = t -> all(t[i].spin_l <= [2, 1][O3.species[i]]
                             for i in eachindex(t)))
        uncapped = _orbit_salcs_decors(xt2, st2, 3, 1, O3, [lab112], false, wc2)
        @test length(uncapped) > length(capped) > 0
        @test [s.key for s in uncapped] != [s.key for s in capped]
        # `isotropy = true` is the same screen in both engines: the pure-spin
        # one keeps Lf == 0, the decor one L_S == 0, and on a pure-spin label
        # L_S ≡ Lf. It must also be the bitwise Lf == 0 SUBSET of the
        # anisotropic build, not a separate construction.
        @test anti_drift(xt2, st2, 3, O3, [2, 1], true, wc2) > 0
        @test anti_drift(xtalB, sgB, 2, O2, [2], true, wcB) > 0
        full = _orbit_salcs_decors(xtalB, sgB, 2, 1, O2,
                                   [spin_decors([1, 1]), spin_decors([2, 2])],
                                   false, wcB)
        iso = _orbit_salcs_decors(xtalB, sgB, 2, 1, O2,
                                  [spin_decors([1, 1]), spin_decors([2, 2])],
                                  true, wcB)
        sub = [s for s in full if s.key.L_S == 0]
        @test 0 < length(iso) < length(full)
        @test [s.key for s in iso] == [s.key for s in sub]
        for (a, b) in zip(iso, sub)
            @test same_members(a.members, b.members)
        end
    end

    @testset "invariant counts vs an independent Cartesian projector" begin
        # Single site, spin l = 2 × disp (k = 0, l = 2): two rank-2 factors on
        # one site under O_h.
        lab_e = [SiteDecor(; spin = 2, disp = (0, 2))]
        se = _orbit_salcs_decors(xtal, sg, 1, 1, O1, [lab_e], false, wc)
        @test length(se) == _count_single_site_l2_invariants(rots, 2)
        @test length(se) == 2                       # E_g ⊕ T_2g, one square each
        @test all(s.key.L_S == 2 for s in se)
        # |u|² trace channel: one rank-2 factor times a scalar ⇒ no invariant.
        lab_t = [SiteDecor(; spin = 2, disp = (1, 0))]
        @test _count_single_site_l2_invariants(rots, 1) == 0
        @test isempty(_orbit_salcs_decors(xtal, sg, 1, 1, O1, [lab_t], false, wc))
        # `isotropy = true` (L_S = 0 only) empties the L_S = 2 sector.
        @test isempty(_orbit_salcs_decors(xtal, sg, 1, 1, O1, [lab_e], true, wc))

        # Both bond sites carry spin l = 1 AND disp (0, 1): the degree-(1,1,1,1)
        # form on the d4h bond.
        lab = [SiteDecor(; spin = 1, disp = (0, 1)),
               SiteDecor(; spin = 1, disp = (0, 1))]
        sall = _orbit_salcs_decors(xtalB, sgB, 2, 1, O2, [lab], false, wcB)
        @test length(sall) ==
              _count_bond_1111_invariants(rotsB, R -> R[1, 1] < 0)
        # The literal is Burnside by hand, and it is NOT redundant with the
        # projector above: `rotsB` feeds both the space group and the oracle, so
        # a broken op list would put them in agreement on the wrong answer. With
        # R = (s₁) ⊕ B over the 8 signed 2×2 permutations B: the 8 site-fixing
        # ops (s₁ = +1) contribute (tr R)⁴ = 81 + 7·1 = 88, and the 8 swapping
        # ops (s₁ = −1) contribute tr(R²)² = 6·9 + 2·1 = 56, so the invariant
        # count is (88 + 56)/16 = 9.
        @test length(sall) == 9
        # The chirality twist (ê₁×ê₂)·(u₁×u₂) lives in L_S = 1 and survives the
        # centrosymmetric bond; `isotropy = true` (L_S = 0) removes it.
        @test any(s.key.L_S == 1 for s in sall)
        s0 = _orbit_salcs_decors(xtalB, sgB, 2, 1, O2, [lab], true, wcB)
        @test all(s.key.L_S == 0 for s in s0)
        subset = [s for s in sall if s.key.L_S == 0]
        @test [s.key for s in s0] == [s.key for s in subset]
        for (a, b) in zip(s0, subset)
            @test same_members(a.members, b.members)
        end
        @test length(s0) < length(sall)
        # and the (L_S = 1, Lf = 0) SALC IS that twist: Φ ∝ (ê₁×ê₂)·(u₁×u₂),
        # a closed form this file writes down independently of the engine.
        tw = only(s for s in sall if s.key.L_S == 1 && s.key.Lf == 0)
        rng2 = MersenneTwister(0x7715)
        ratios = Float64[]
        for _ = 1:6
            e = reduce(hcat, [normalize(randn(rng2, 3)) for _ = 1:2])
            u = randn(rng2, 3, 2) * 0.4
            F = dot(cross(e[:, 1], e[:, 2]), cross(u[:, 1], u[:, 2]))
            abs(F) < 1e-3 && continue
            push!(ratios, evaluate_salc(tw, e, u) / F)
        end
        @test !isempty(ratios)
        @test all(r -> isapprox(r, ratios[1]; rtol = 1e-9), ratios)
    end

    @testset "u = 0 degeneracy, pure-spin consistency, and the refusals" begin
        lab = [SiteDecor(; spin = 1, disp = (0, 1)),
               SiteDecor(; spin = 1, disp = (0, 1))]
        sall = _orbit_salcs_decors(xtalB, sgB, 2, 1, O2, [lab], false, wcB)
        e = reduce(hcat, [normalize(randn(rng, 3)) for _ = 1:2])
        u0 = zeros(3, 2)
        for s in sall
            @test evaluate_salc(s, e, u0) == 0.0     # exact: homogeneous ≥ 1
        end
        # pure-spin SALCs evaluate identically through the joint form, ∀u
        pure = _orbit_salcs_decors(xtalB, sgB, 2, 1, O2, [spin_decors([1, 1])],
                                   false, wcB)
        u = randn(rng, 3, 2) * 0.3
        for s in pure
            @test evaluate_salc(s, e, u) === evaluate_salc(s, e)
        end
        # the spin-only forms refuse a decorated SALC instead of mis-scaling
        mixed = first(sall)
        @test_throws ArgumentError evaluate_salc(mixed, e)
        @test_throws ArgumentError SCEFitting.accumulate_grad!(zeros(3, 2), mixed,
                                                               e, 1.0)
        @test_throws ArgumentError evaluate_salc(mixed, e, zeros(3, 1))  # size
        # duplicate labels are rejected (collinear-column guard)
        lab2 = [SiteDecor(; spin = 1, disp = (0, 1)),
                SiteDecor(; spin = 1, disp = (0, 1))]
        @test_throws ArgumentError _orbit_salcs_decors(xtalB, sgB, 2, 1, O2,
                                                       [lab2, lab2], false, wcB)
    end

    @testset "mixed space-group invariance + time reversal" begin
        lab = [SiteDecor(; spin = 1, disp = (0, 1)),
               SiteDecor(; spin = 1, disp = (0, 1))]
        sall = _orbit_salcs_decors(xtalB, sgB, 2, 1, O2, [lab], false, wcB)
        lab1 = [SiteDecor(; spin = 2, disp = (0, 2))]
        s1 = _orbit_salcs_decors(xtal, sg, 1, 1, O1, [lab1], false, wc)
        for _ = 1:6
            # bond crystal: ops with x → −x swap the two atom columns
            e = reduce(hcat, [normalize(randn(rng, 3)) for _ = 1:2])
            u = randn(rng, 3, 2) * 0.4
            for s in sall
                v0 = evaluate_salc(s, e, u)
                for R in rotsB
                    p = R[1, 1] > 0 ? (1:2) : (2:-1:1)
                    eg = det(R) * R * e[:, p]      # axial spin action + site map
                    ug = R * u[:, p]               # polar displacement action
                    @test evaluate_salc(s, Matrix(eg), Matrix(ug)) ≈ v0 atol = 1e-10
                end
                @test evaluate_salc(s, -e, u) ≈ v0 atol = 1e-12
            end
            # single-site crystal: all 48 ops fix the lone atom
            e1 = reshape(normalize(randn(rng, 3)), 3, 1)
            u1 = randn(rng, 3, 1) * 0.4
            for s in s1
                v0 = evaluate_salc(s, e1, u1)
                for R in _oh48()
                    @test evaluate_salc(s, Matrix(det(R) * R * e1),
                                        Matrix(R * u1)) ≈ v0 atol = 1e-10
                end
                @test evaluate_salc(s, -e1, u1) ≈ v0 atol = 1e-12
            end
        end
    end

    @testset "label validation" begin
        @test_throws ArgumentError _orbit_salcs_decors(xtalB, sgB, 2, 1, O2,
            [[SiteDecor(; spin = 1), SiteDecor(; spin = 2)]], false, wcB)  # Σl odd
        @test_throws ArgumentError _orbit_salcs_decors(xtalB, sgB, 2, 1, O2,
            [[SiteDecor(; spin = 1)]], false, wcB)                         # wrong N
        @test_throws ArgumentError _orbit_salcs_decors(xtalB, sgB, 2, 1, O2,
            [[SiteDecor(; spin = 2), SiteDecor(; spin = 1, disp = (0, 1))]],
            false, wcB)                                                    # unsorted
    end

    @testset "group_costs refuses a decorated basis" begin
        # the MC entry-key model prices SPIN axes; a decorated basis has none of
        # its own contract, so the surface refuses rather than under-counting.
        lab = [SiteDecor(; spin = 1, disp = (0, 1)),
               SiteDecor(; spin = 1, disp = (0, 1))]
        sall = _orbit_salcs_decors(xtalB, sgB, 2, 1, O2, [lab], false, wcB)
        sb = SCEFitting.SALCBasis(sall, [s.key for s in sall])
        spec = SCEFitting.BasisSpec(xtalB; nbody = 2, lmax = 1, cutoff = 1.1)
        @test_throws ArgumentError SCEFitting.group_costs(
            SCEFitting.SCEBasis(xtalB, sgB, sb, spec))
    end

    @testset "joint evaluation values: Parseval closed form + hand contraction" begin
        # The gates above never pin a VALUE of the mixed kernel: the twist
        # `(ê₁×ê₂)·(u₁×u₂)` is symmetric under e ↔ u and under the site swap, the
        # u = 0 test is exact for any homogeneous kernel, and invariance holds
        # for any scale. So a kernel that dropped `|u|^{2k}` (the pointed mark is
        # k = 1), scaled by `(4π)^(n_slots/2)` instead of `(4π)^(n_spin/2)` (the
        # mark-only site has n_spin = 0), or read the displacement of the wrong
        # site would pass them all. This testset fixes the value two ways.
        #
        # (1) Parseval. Under the trivial group every coupling path is invariant,
        # so a label's SALCs form an orthonormal basis of the full product space
        # and, for one ordered instance, Σ_s Φ_s² = ‖⊗ factor tables‖² · (4π)^n_spin
        # — by the addition theorems Σ_m Z_{lm}(ê)² = (2l+1)/4π and
        # Σ_m R_{lm}(u)² = |u|^{2l} (Racah-normalised solid harmonics), that is
        #     Π_spin (2l+1) · Π_disp |u_a|^{2l + 4k},
        # with no folded tensor, Clebsch–Gordan coefficient, or Wigner matrix in
        # sight. One documented convention enters: the same physical cluster
        # arrives as N! ordered images which `_canonicalize_members` SUMS into one
        # member (`basis/salc.jl`), so the value carries N! and the sum (N!)².
        # Under C₁ the two-site orbit is exactly those images (asserted below).
        # Each distinct site assignment of a label is its own orthonormal block,
        # so the closed form sums over the label's arrangements.
        #
        # (2) A per-SALC contraction written from the public kernels `Zlm` / `Rlm`
        # and the slot list, which pins the channel dispatch, the `m ↔ axis`
        # offset, `|u|^{2k}` and the site lookup of `_fill_ztables_mixed!` term
        # by term (it shares `folded` with the engine, so it is a check on the
        # evaluation, not on the projection).
        I3 = SMatrix{3,3,Float64}(I)
        z3 = SVector{3,Float64}(0, 0, 0)
        sgP1 = _assemble_spacegroup(xtal, [I3], [z3], "P1", 1; tol = 1e-5)
        wcP1 = _build_wig_cache(sgP1, 2)
        O1P1 = build_clusters(xtal, build_neighbor_list(xtal, 2.1), sgP1;
                              nbody = 1).by_body[1][1]
        sgP1B = _assemble_spacegroup(xtalB, [I3], [z3], "P1", 1; tol = 1e-5)
        wcP1B = _build_wig_cache(sgP1B, 2)
        O2P1 = build_clusters(xtalB, build_neighbor_list(xtalB, 1.1), sgP1B;
                              nbody = 2).by_body[2][1]
        # the orbit is the one bond as its 2! ordered images, nothing else
        @test length(O2P1.members) == 2
        @test sort([m.atoms for m in O2P1.members]) == [[1, 2], [2, 1]]

        function closed_form(label, u)
            N = length(label)
            perms = N == 1 ? [[1]] : [[1, 2], [2, 1]]
            arrangements = unique([label[p] for p in perms])
            tot = 0.0
            for a in arrangements
                f = 1.0
                for s in eachindex(a)
                    SCEFitting.has_spin(a[s]) && (f *= 2 * a[s].spin_l + 1)
                    SCEFitting.has_disp(a[s]) &&
                        (f *= norm(u[:, s])^(2 * a[s].disp_l + 4 * a[s].disp_k))
                end
                tot += f
            end
            return factorial(N)^2 * tot
        end

        function hand_value(s, e, u)
            n_spin = count(SCEFitting.has_spin, s.decors)
            tot = 0.0
            for m in s.members, t in m.terms
                for idx in CartesianIndices(t.folded)
                    w = t.folded[idx]
                    w == 0.0 && continue
                    for (i, sl) in enumerate(t.slots)
                        a = m.atoms[sl.site]
                        l = sl.factor.l
                        μ = idx[i] - l - 1
                        if sl.factor.channel == SPIN
                            w *= SCEFitting.Harmonics.Zlm(l, μ, e[:, a])
                        else
                            w *= norm(u[:, a])^(2 * sl.factor.k) *
                                 SCEFitting.SolidHarmonics.Rlm(l, μ, u[:, a])
                        end
                    end
                    tot += w
                end
            end
            return (4π)^(n_spin / 2) * tot
        end

        # Labels chosen so every escape changes the closed form: a k = 1 mark on
        # a ranked site, the bare pointed mark (n_spin = 0, Σ = |u|⁴), a k = 0
        # rank-2 decor, and the two-site shapes — mark and rank on DIFFERENT
        # sites (n_spin = 1 of 2 slots), and both sites decorated.
        cases = [
            (xtal, sgP1, 1, O1P1, wcP1, [SiteDecor(; spin = 2, disp = (1, 1))]),
            (xtal, sgP1, 1, O1P1, wcP1, [SiteDecor(; disp = (1, 0))]),
            (xtal, sgP1, 1, O1P1, wcP1, [SiteDecor(; spin = 2, disp = (0, 2))]),
            (xtalB, sgP1B, 2, O2P1, wcP1B,
             sort([SiteDecor(; spin = 1), SiteDecor(; spin = 1, disp = (1, 1))])),
            (xtalB, sgP1B, 2, O2P1, wcP1B,
             sort([SiteDecor(; disp = (1, 0)), SiteDecor(; spin = 2)])),
            (xtalB, sgP1B, 2, O2P1, wcP1B,
             [SiteDecor(; spin = 1, disp = (0, 1)),
              SiteDecor(; spin = 1, disp = (0, 1))]),
        ]
        rngP = MersenneTwister(0x9a25)
        for (xt, sgx, N, O, wcx, label) in cases
            ss = _orbit_salcs_decors(xt, sgx, N, 1, O, [label], false, wcx)
            # completeness: Π(2l+1) over every factor, per arrangement
            nfull = sum(prod((SCEFitting.has_spin(d) ? 2d.spin_l + 1 : 1) *
                             (SCEFitting.has_disp(d) ? 2d.disp_l + 1 : 1) for d in a)
                        for a in unique([label[p] for p in
                                         (N == 1 ? [[1]] : [[1, 2], [2, 1]])]))
            @test length(ss) == nfull
            for _ = 1:3
                e = reduce(hcat, [normalize(randn(rngP, 3)) for _ = 1:N])
                # column norms kept away from 1 and from each other, so a dropped
                # `|u|^{2k}` or a swapped site cannot hide behind |u| ≈ 1
                u = reduce(hcat, [normalize(randn(rngP, 3)) * r
                                  for r in (N == 1 ? (0.6,) : (0.6, 1.7))])
                P = sum(evaluate_salc(s, e, u)^2 for s in ss)
                @test P ≈ closed_form(label, u) rtol = 1e-12
                for s in ss
                    v = evaluate_salc(s, e, u)
                    @test v ≈ hand_value(s, e, u) atol = 1e-12 * max(1.0, abs(v))
                end
            end
        end
        # the (4π) exponent counts SPIN slots, not slots: the spin-free mark
        # evaluates to |u|² R₀₀ = |u|² exactly (here 0.09 + 0.16 + 1.44), no 4π
        sm = only(_orbit_salcs_decors(xtal, sgP1, 1, 1, O1P1,
                                      [[SiteDecor(; disp = (1, 0))]], false, wcP1))
        um = reshape([0.3, -0.4, 1.2], 3, 1)
        @test evaluate_salc(sm, reshape([0.0, 0.0, 1.0], 3, 1), um) ≈ 1.69 rtol = 1e-13
    end
end
