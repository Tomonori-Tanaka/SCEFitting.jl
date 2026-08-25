# The pointed (site-marked) moment basis (src/basis/momentbasis.jl): MomentSpec
# validation, the FeGe B20 oracles ported from the design-record prototype (star
# closed form against an INDEPENDENT geometric enumeration, the 2√3 shell-sum
# normalization, G_i covariance including the evaluation axes, time reversal,
# marked-column substitution locality), and the pointed resolvability gate
# (signature rank ≡ independent random-design rank; null combinations annihilate
# the actual design; the repeated-image refusal). Ported from SLCE.jl with the
# `isotropy` spelling of the screen.

using Test
using SCEFitting
using SCEFitting: _assemble_spacegroup, n_ops, cartesian_positions, _design_moment,
                  salcs, n_salcs, has_disp, _orbit_salcs_decors, _build_wig_cache,
                  SiteDecor, UnclassifiableBasis, SpaceGroup, _orbits_from_members,
                  build_clusters, build_neighbor_list, candidate_clusters
using LinearAlgebra
using StaticArrays
using Random

# ── FeGe B20 primitive with the handwritten P2₁3 group (no Spglib in the core) ──
const _MB_WS = [[1 0 0; 0 1 0; 0 0 1], [-1 0 0; 0 -1 0; 0 0 1],
                [-1 0 0; 0 1 0; 0 0 -1], [1 0 0; 0 -1 0; 0 0 -1],
                [0 0 1; 1 0 0; 0 1 0], [0 0 1; -1 0 0; 0 -1 0],
                [0 0 -1; -1 0 0; 0 1 0], [0 0 -1; 1 0 0; 0 -1 0],
                [0 1 0; 0 0 1; 1 0 0], [0 -1 0; 0 0 1; -1 0 0],
                [0 1 0; 0 0 -1; -1 0 0], [0 -1 0; 0 0 -1; 1 0 0]]
const _MB_TS = [[0, 0, 0], [0.5, 0, 0.5], [0, 0.5, 0.5], [0.5, 0.5, 0],
                [0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                [0, 0, 0], [0, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5]]

_mb_wrap01(x) = (y = mod(x, 1.0); y >= 1.0 - 1e-12 ? 0.0 : y)

function _mb_orbit4a(x::Float64)
    out = Vector{Vector{Float64}}()
    for (W, t) in zip(_MB_WS, _MB_TS)
        p = _mb_wrap01.(W * [x, x, x] .+ t)
        any(q -> maximum(abs.(q .- p)) < 1e-9, out) || push!(out, p)
    end
    sort!(out)
    return hcat(out...)
end

function _mb_fege()
    pos = hcat(_mb_orbit4a(0.1352), _mb_orbit4a(0.8414))
    xt = Crystal(Lattice(Matrix(4.7 * I(3))), pos, [1, 1, 1, 1, 2, 2, 2, 2],
                 ["Fe", "Ge"])
    sg = _assemble_spacegroup(xt,
        [SMatrix{3,3,Float64}(Float64.(W)) for W in _MB_WS],
        [SVector{3,Float64}(Float64.(t)) for t in _MB_TS], "P2_13", 198; tol = 1e-5)
    return xt, sg
end

struct _MBFixedSG <: SCEFitting.AbstractSymmetryBackend
    sg::SpaceGroup
end
SCEFitting.analyze_symmetry(b::_MBFixedSG, ::Crystal; tol::Real = 1e-5) = b.sg

_mb_unit(rng, nat) = (m = randn(rng, 3, nat);
                      for a = 1:nat; m[:, a] ./= norm(m[:, a]); end; m)

@testset "pointed moment basis" begin
    xt, sg = _mb_fege()
    bk = _MBFixedSG(sg)
    nat = 8
    cart = Matrix(cartesian_positions(xt))
    A = Matrix(xt.lattice.vectors)

    @testset "MomentSpec validation" begin
        ok = MomentSpec(; lmax_env = [2, 0], sampled = [true, false],
                        cutoff_pair = 3.0)
        @test ok.nbody == 3 && ok.cutoff_star == [ok.cutoff_pair]
        @test ok.isotropy                           # L_S = 0 only, by default
        # per-star-order radii: one matrix per order, body N at index N - 2. A scalar
        # or a matrix broadcasts to every order; a vector must have exactly nbody - 2
        # entries; below body order 3 there is no star to cut, so a radius that would
        # never be read is refused rather than stored.
        pb = MomentSpec(; lmax_env = [2], sampled = [true], nbody = 4,
                        cutoff_pair = 3.0, cutoff_star = [3.0, [1.5;;]])
        @test length(pb.cutoff_star) == 2
        @test pb.cutoff_star[1] == fill(3.0, 1, 1) && pb.cutoff_star[2] == fill(1.5, 1, 1)
        @test SCEFitting._star_cutoff(pb, 4) == pb.cutoff_star[2]
        @test SCEFitting._star_cutoff_envelope(pb) == fill(3.0, 1, 1)
        bcast = MomentSpec(; lmax_env = [2], sampled = [true], nbody = 4,
                           cutoff_pair = 3.0, cutoff_star = 2.0)
        @test bcast.cutoff_star == [fill(2.0, 1, 1), fill(2.0, 1, 1)]
        @test_throws ArgumentError MomentSpec(; lmax_env = [2], sampled = [true],
                                              nbody = 4, cutoff_pair = 3.0,
                                              cutoff_star = [3.0])        # too few
        @test_throws ArgumentError MomentSpec(; lmax_env = [2], sampled = [true],
                                              nbody = 3, cutoff_pair = 3.0,
                                              cutoff_star = [3.0, 1.5])   # too many
        @test_throws ArgumentError MomentSpec(; lmax_env = [2], sampled = [true],
                                              nbody = 2, cutoff_pair = 3.0,
                                              cutoff_star = 3.0)          # no star order
        @test_throws ArgumentError MomentSpec(; lmax_env = [2], sampled = [true],
                                              nbody = 4, cutoff_pair = 3.0,
                                              cutoff_star = [3.0, -1.0])  # negative
        # the M3-1 assert: environment spin factors only on sampled species
        @test_throws ArgumentError MomentSpec(; lmax_env = [2, 1],
                                              sampled = [true, false],
                                              cutoff_pair = 3.0)
        @test_throws ArgumentError MomentSpec(; lmax_env = [2], sampled = [true],
                                              cutoff_pair = 3.0, lmax_mark = -1)
        @test_throws ArgumentError MomentSpec(; lmax_env = [2], sampled = [true],
                                              cutoff_pair = 3.0, nbody = 5)
        # 4 is inside the door: the enumeration is general in N, the cap is where the
        # oracles stop
        @test MomentSpec(; lmax_env = [2], sampled = [true], cutoff_pair = 3.0,
                         nbody = 4).nbody == 4
        @test_throws ArgumentError MomentSpec(; lmax_env = [2], sampled = [true],
                                              cutoff_pair = 3.0,
                                              marked = [false])
        @test_throws ArgumentError MomentSpec(; lmax_env = [2, 2],
                                              sampled = [true, true],
                                              cutoff_pair = [3.0 2.0; 2.5 3.0])
        # species-count mismatch dies at the basis door
        @test_throws ArgumentError MomentBasis(xt,
            MomentSpec(; lmax_env = [2], sampled = [true], cutoff_pair = 3.0);
            backend = bk)
        # a cutoff below every shell leaves exactly the 1-body content (the
        # per-orbit intercepts μ₀ and single-site ê invariants never need a bond)
        tiny = MomentBasis(xt, MomentSpec(; lmax_env = [1, 1],
                                          sampled = [true, true],
                                          cutoff_pair = 0.5); backend = bk)
        @test all(k -> k.body == 1, tiny.salc_basis.keys)
    end

    @testset "_orbits_from_members is build_clusters' core (factoring is inert)" begin
        nl = build_neighbor_list(xt, 3.3)
        cs = build_clusters(xt, nl, sg; nbody = 3)
        cand = candidate_clusters(xt, nl, 3)
        for body in (1, 2, 3)
            direct = _orbits_from_members(xt, sg, cand[body], body)
            @test length(direct) == length(cs.by_body[body])
            for (a, b) in zip(direct, cs.by_body[body])
                @test a.representative.atoms == b.representative.atoms
                @test a.representative.shifts == b.representative.shifts
                @test [m.atoms for m in a.members] == [m.atoms for m in b.members]
                @test a.species == b.species
            end
        end
    end

    spec = MomentSpec(; lmax_env = [2, 2], sampled = [true, true], lmax_mark = 2,
                      nbody = 3, cutoff_pair = 3.3, cutoff_star = 3.3)
    mb = MomentBasis(xt, spec; backend = bk)
    ks = mb.salc_basis.keys
    mark_l(k) = (i = findfirst(has_disp, k.decors); k.decors[i].spin_l)
    env_ls(k) = sort([d.spin_l for d in k.decors if !has_disp(d)])

    rng = MersenneTwister(0x20260819)
    e = _mb_unit(rng, nat)

    # independent geometric references (ported from the design-record prototype 01;
    # nothing here touches the SALC machinery)
    imgs = [A * Float64.([i, j, k]) for i = -1:1, j = -1:1, k = -1:1]
    function nn_neighbors(a; d0 = 2.881, tol = 0.05)
        nbrs = Tuple{Int,Vector{Float64}}[]
        for j = 1:4, R in imgs
            p = cart[:, j] + R
            abs(norm(p - cart[:, a]) - d0) < tol && push!(nbrs, (j, p))
        end
        return nbrs
    end
    star_ref(a) = begin
        nbrs = nn_neighbors(a)
        acc = 0.0
        for i = 1:length(nbrs), j = (i + 1):length(nbrs)
            (jj, pj) = nbrs[i]
            (kk, pk) = nbrs[j]
            abs(norm(pj - pk) - 2.881) < 0.05 || continue
            acc += dot(e[:, jj], e[:, kk])
        end
        acc
    end
    p1_ref(a) = sum(dot(e[:, a], e[:, j]) for (j, _) in nn_neighbors(a))

    @testset "FeGe oracles: closed forms vs independent geometry" begin
        @test mb.marked_atoms == collect(1:8)
        @test all(s.key.L_S == 0 for s in salcs(mb))   # the isotropy screen held
        X1 = _design_moment(mb, [e], [copy(e)])          # mode-4 identity axes
        # the fast (mark → term index) path is value-identical to the full one
        @test X1 == _design_moment(mb, [e], [copy(e)]; member_index = false)
        # star (0,1,1) on the closed nn Fe₃ triangle: SALC = 6 · Σ e_j·e_k over the
        # triangles at atom a — the prototype's headline geometric oracle
        jstars = [j for j in 1:length(ks) if ks[j].body == 3 && mark_l(ks[j]) == 0 &&
                  env_ls(ks[j]) == [1, 1] && all(==(1), mb.records[j].species) &&
                  all(abs(x - 2.881) < 0.05 for x in mb.records[j].edges)]
        @test length(jstars) == 1
        # The INDEPENDENT oracle is the geometry: the SALC column is one constant
        # times the triangle sum, the same constant on every atom and on a second
        # random configuration (a missing ordering or a doubled triangle breaks
        # this ratio, which is what caught the enumeration bug upstream).
        rs = [X1[a, jstars[1]] / star_ref(a) for a = 1:4]
        @test all(r -> isapprox(r, rs[1]; rtol = 1e-10), rs)
        e_b = _mb_unit(rng, nat)
        X1b = _design_moment(mb, [e_b], [copy(e_b)])
        star_ref_b(a) = begin
            nbrs = nn_neighbors(a)
            acc = 0.0
            for i = 1:length(nbrs), j = (i + 1):length(nbrs)
                (jj, pj) = nbrs[i]; (kk, pk) = nbrs[j]
                abs(norm(pj - pk) - 2.881) < 0.05 || continue
                acc += dot(e_b[:, jj], e_b[:, kk])
            end
            acc
        end
        @test all(isapprox(X1b[a, jstars[1]] / star_ref_b(a), rs[1]; rtol = 1e-10)
                  for a = 1:4)
        # PIN (change detector, not an oracle): the constant 6.0 was read off the
        # design-record prototype's run of the upstream machinery
        # (prototype/01_basis_oracles.jl, SLCE-pinned 2026-08-19). Recapture only
        # on an explicit decision to change the SALC normalization.
        @test rs[1] ≈ 6.0 rtol = 1e-10
        # nn Fe–Fe pointed (1,1): exactly 2 SALC blocks (the C₃ split of the 6-shell
        # into 3+3), whose SUM is the 2√3-normalized shell sum
        jnn = [j for j in 1:length(ks) if ks[j].body == 2 && mark_l(ks[j]) == 1 &&
               env_ls(ks[j]) == [1] && all(==(1), mb.records[j].species) &&
               abs(mb.records[j].edges[1] - 2.881) < 0.05]
        @test length(jnn) == 2
        # same discipline: the shell sum is the oracle (ratio constant across
        # atoms and configurations); 2√3 is the pinned constant (change detector,
        # captured from the same prototype run)
        rp = [sum(X1[a, j] for j in jnn) / p1_ref(a) for a = 1:4]
        @test all(r -> isapprox(r, rp[1]; rtol = 1e-10), rp)
        p1_ref_b(a) = sum(dot(e_b[:, a], e_b[:, j]) for (j, _) in nn_neighbors(a))
        @test all(isapprox(sum(X1b[a, j] for j in jnn) / p1_ref_b(a), rp[1];
                           rtol = 1e-10) for a = 1:4)
        @test rp[1] ≈ 2 * sqrt(3) rtol = 1e-10
        # the marked-multiplicity structure of the nn pointed pair: 12 members
        # (6 bonds × 2 orderings), mark histogram 3 per Fe atom of the split
        s = salcs(mb)[jnn[1]]
        @test length(s.members) == 12
        marks = Int[]
        for m in s.members, t in m.terms, sl in t.slots
            sl.factor.channel == SCEFitting.DISP && push!(marks, m.atoms[sl.site])
        end
        @test sort(unique(marks)) == [1, 2, 3, 4]
        @test all(count(==(a), marks) == length(marks) ÷ 4 for a = 1:4)
    end

    @testset "covariance / TR / substitution locality" begin
        axr = _mb_unit(rng, nat)
        Xr = _design_moment(mb, [e], [axr])
        # G_i covariance including the axes: Φ(g·a; g∘e, g∘ê) = Φ(a; e, ê)
        dev = 0.0
        for g = 1:n_ops(sg)
            R = sg.ops[g].rotation_cart
            eg = similar(e)
            axg = similar(axr)
            for a = 1:nat
                eg[:, sg.map_sym[a, g]] = R * e[:, a]
                axg[:, sg.map_sym[a, g]] = R * axr[:, a]
            end
            Xg = _design_moment(mb, [eg], [axg])
            for a = 1:nat, j = 1:size(Xr, 2)
                dev = max(dev, abs(Xg[sg.map_sym[a, g], j] - Xr[a, j]))
            end
        end
        @test dev < 1e-12
        # An ARBITRARY rotation, not just a space-group operation: an L_S = 0 column is
        # an isotropic invariant, so it survives any R ∈ SO(3) applied to the spins AND
        # the axes together. Rotating only the spins leaves the mark axis behind, and
        # this basis's mark carries rank up to 2, so the columns then move — that
        # asymmetry is the whole content of the marked-column substitution.
        let iso = findall(k -> k.L_S == 0, mb.salc_basis.keys)
            q = qr(randn(rng, 3, 3))
            R = Matrix(q.Q) * (det(Matrix(q.Q)) < 0 ? Diagonal([-1.0, 1, 1]) : I)
            @test det(R) ≈ 1.0 rtol = 1e-12
            @test _design_moment(mb, [R * e], [R * axr])[:, iso] ≈
                  Xr[:, iso] rtol = 1e-12
            @test !isapprox(_design_moment(mb, [R * e], [axr])[:, iso], Xr[:, iso];
                            rtol = 1e-6)
        end
        # time reversal is bitwise (every label has even total spin rank)
        @test _design_moment(mb, [-e], [-axr]) == Xr
        # marked-column substitution locality: changing the axis of atom b changes
        # rows of b and nothing else (the exactness of the substitution)
        ax2 = copy(axr)
        ax2[:, 3] = normalize(randn(rng, 3))
        X2 = _design_moment(mb, [e], [ax2])
        @test findall([any(X2[a, :] .!= Xr[a, :]) for a = 1:nat]) == [3]
        # mode-4 identity: substituting e's own columns is a no-op
        @test _design_moment(mb, [e], [copy(e)]) ==
              _design_moment(mb, [e], [Matrix(e)])
        # deterministic rebuild
        mb2 = MomentBasis(xt, spec; backend = bk)
        @test mb2.salc_basis.keys == mb.salc_basis.keys
        # isotropy = false keeps every L_S and contains the isotropic basis
        full = MomentBasis(xt, MomentSpec(; lmax_env = [2, 2], sampled = [true, true],
                                          lmax_mark = 2, nbody = 2, cutoff_pair = 3.3,
                                          isotropy = false); backend = bk)
        iso = MomentBasis(xt, MomentSpec(; lmax_env = [2, 2], sampled = [true, true],
                                         lmax_mark = 2, nbody = 2, cutoff_pair = 3.3);
                          backend = bk)
        @test any(k.L_S != 0 for k in full.salc_basis.keys)
        @test [k for k in full.salc_basis.keys if k.L_S == 0] == iso.salc_basis.keys
    end

    @testset "pointed resolvability gate" begin
        # pairs-only basis: no repeated-image members, so the gate classifies —
        # and its symbolic rank must equal the rank of an INDEPENDENT random
        # design (two implementations of one question)
        psp = MomentSpec(; lmax_env = [2, 2], sampled = [true, true],
                         lmax_mark = 2, nbody = 2, cutoff_pair = 3.3)
        pmb = MomentBasis(xt, psp; backend = bk)
        res = moment_resolvability(pmb)
        cfgs = [_mb_unit(rng, nat) for _ = 1:60]
        axes = [_mb_unit(rng, nat) for _ = 1:60]
        X = _design_moment(pmb, cfgs, axes)
        sv = svd(X).S
        @test count(>(1e-9 * sv[1]), sv) == res.rank
        # the two thresholds (1e-9 here, the gate's 1e-10) agree because the
        # spectrum has a clear gap at the rank — asserted, so the equality above
        # is not threshold luck
        @test res.rank == length(sv) || sv[res.rank] / sv[res.rank + 1] > 1e3
        # every symbolic null combination annihilates the actual design
        for comb in res.null_combinations
            v = zeros(n_salcs(pmb))
            for (j, w) in comb
                v[j] = w
            end
            @test norm(X * v) < 1e-10 * norm(X) * norm(v)
        end
        # the null report NAMES columns (nonempty combinations) when rank-deficient
        res.rank < length(res.kept) &&
            @test all(!isempty, res.null_combinations)
        # census: every pointed pair orbit on this cell carries ≥ 2 mark classes
        # (both ends markable) — the face-(b) hazard preregistration
        @test all(c -> c.n_mark_atoms >= 2, res.census)
        # default-rtol cache: the SAME object comes back (=== is the contract);
        # a non-default rtol always recomputes and is never cached; a rebuilt
        # basis recomputes from scratch and agrees (two computations, one answer)
        @test moment_resolvability(pmb) === res
        fresh = moment_resolvability(pmb; rtol = 1e-9)
        @test fresh !== res
        @test fresh.vanishing == res.vanishing     # rtol-independent classification
        pmb2 = MomentBasis(xt, psp; backend = bk)
        res2 = moment_resolvability(pmb2)
        @test res2 !== res && res2.rank == res.rank &&
              res2.vanishing == res.vanishing
        # The star basis on the PRIMITIVE cell used to fold two images of one
        # neighbour into a single environment sphere, and the gate refused the WHOLE
        # basis for it. `_pointed_star_candidates` now refuses that shape at the
        # enumeration, so the cell keeps exactly the stars it can resolve and the gate
        # classifies. Not a weakening: what left were Clebsch–Gordan-reducible products
        # on one sphere — never independent N-body functions — and their presence used
        # to cost every sound column alongside them.
        @test all(allunique(m.atoms) for s in salcs(mb) for m in s.members)
        @test moment_resolvability(mb) isa NamedTuple
        # a fourth spoke: same rule, same outcome at N = 4
        mb4f = MomentBasis(xt, MomentSpec(; lmax_env = [2, 2], sampled = [true, true],
                                          lmax_mark = 2, nbody = 4, cutoff_pair = 3.3,
                                          cutoff_star = 3.3, lsum = 4); backend = bk)
        @test all(allunique(m.atoms) for s in salcs(mb4f) for m in s.members)
        @test moment_resolvability(mb4f) isa NamedTuple
        @test occursin("UnclassifiableBasis", sprint(showerror,
                                                     UnclassifiableBasis("x")))
        # ... and the classifier's guard is still there, for a candidate source that
        # does not apply the enumeration's rule
        @test_throws UnclassifiableBasis moment_resolvability(
            fold_members_onto_one_atom(mb))
    end

    @testset "a star never repeats a reference-cell atom" begin
        # Two cells of ONE geometry — a chain of atoms 1.5 Å apart along x — differing
        # only in how much of it the reference cell holds. The counts are hand-derived
        # from that geometry, not read off the code.
        #
        # L = 3 Å, 2 atoms: atom 1's only neighbours within 2 Å are atom 2 at shift
        #   (0,0,0) and at (-1,0,0) — one atom reached two ways. The single 2-of-2
        #   choice therefore repeats it, and there is NO resolvable 3-body star: both
        #   environment factors would read the one spin `e₂`.
        # L = 6 Å, 4 atoms: the same two spokes now land on DIFFERENT atoms (2 and 4).
        #   Each of the four atoms anchors one triple, the four triples are distinct
        #   atom sets, and each is emitted once per site ordering: 4 × 3! = 24.
        #
        # (Before the enumeration refused repeats: 12 and 24. The small cell emitted
        # 2 translation classes × 3! of members that the dataset door then refused
        # wholesale, taking the sound 1- and 2-body columns with them.)
        chain(L, xs) = (n = length(xs); f = zeros(3, n); f[1, :] = xs ./ L;
                        Crystal(Lattice(Matrix(L * I(3))), f, fill(1, n), ["Fe"]))
        csp = MomentSpec(; lmax_env = [2], sampled = [true], lmax_mark = 2, nbody = 3,
                         cutoff_pair = 2.0, cutoff_star = 2.0, marked = [true])
        function stars(cr)
            nl = SCEFitting.build_neighbor_list(
                                     cr, SCEFitting._star_cutoff_envelope(csp),
                                     SCEFitting.MinimumImage();
                                     tol = SCEFitting._SAME_DIST_RTOL)
            return nl, SCEFitting._pointed_star_candidates(cr, nl, csp, 3)
        end
        small = chain(3.0, [0.0, 1.5])
        big = chain(6.0, [0.0, 1.5, 3.0, 4.5])
        nls, ss = stars(small)
        nlb, sb = stars(big)
        # the fixtures really are the two cases claimed: one atom reached twice, vs
        # two different atoms
        @test sort([p.j for p in nls.pairs if p.i == 1]) == [2, 2]
        @test sort([p.j for p in nlb.pairs if p.i == 1]) == [2, 4]
        @test isempty(ss)
        @test length(sb) == 24
        @test all(m -> allunique(m.atoms), sb)
        # end to end: the small cell keeps the sectors it can resolve and says so,
        # rather than building a 3-body sector its dataset door would refuse
        mbs = @test_logs((:warn, r"body order 3 contributes no SALC"),
                         match_mode = :any, MomentBasis(small, csp))
        @test sort(unique(k.body for k in mbs.salc_basis.keys)) == [1, 2]
        @test moment_resolvability(mbs) isa NamedTuple
        mbb = MomentBasis(big, csp)
        @test sort(unique(k.body for k in mbb.salc_basis.keys)) == [1, 2, 3]
        @test moment_resolvability(mbb) isa NamedTuple
    end

    @testset "saboteur closures: marked ≠ 1:n, species rules, labels, census, chain" begin
        # (1) marked_atoms ≠ 1:n — every other fixture has ai == a, so a slot/row
        # index used as an atom number was invisible. Ge only: marked 5:8.
        sge = MomentSpec(; lmax_env = [1, 1], sampled = [true, true], lmax_mark = 1,
                         nbody = 2, cutoff_pair = 3.0, marked = [false, true])
        mge = MomentBasis(xt, sge; backend = bk)
        @test mge.marked_atoms == [5, 6, 7, 8]
        e2 = _mb_unit(rng, nat)
        axg = _mb_unit(rng, nat)
        Xg = _design_moment(mge, [e2], [axg])
        @test size(Xg) == (4, n_salcs(mge))
        @test Xg == _design_moment(mge, [e2], [axg]; member_index = false)
        # locality in atom numbers: changing the axis of atom 7 changes row 3 only
        ax7 = copy(axg); ax7[:, 7] = normalize(randn(rng, 3))
        X7 = _design_moment(mge, [e2], [ax7])
        @test findall([any(X7[r, :] .!= Xg[r, :]) for r = 1:4]) == [3]
        # and changing the axis of an UNMARKED atom (Fe 2) changes nothing
        ax2 = copy(axg); ax2[:, 2] = normalize(randn(rng, 3))
        @test _design_moment(mge, [e2], [ax2]) == Xg
        # every mark in this basis sits on a Ge atom
        for sx in salcs(mge), m in sx.members, t in m.terms, sl in t.slots
            sl.factor.channel == SCEFitting.DISP && @test m.atoms[sl.site] in 5:8
        end
        # (2) the row order is configuration-major (the docstring's contract):
        # a two-config design is the single-config designs stacked
        e3 = _mb_unit(rng, nat)
        @test _design_moment(mb, [e2, e3], [copy(e2), copy(e3)]) ==
              vcat(_design_moment(mb, [e2], [copy(e2)]),
                   _design_moment(mb, [e3], [copy(e3)]))
        # (3) the M3-1 species rule on a REAL basis: lmax_env = [2, 0] — Ge is
        # never an environment spin, yet Ge is still marked (its induced moment is
        # what the channel predicts). Read off the slots, not the keys.
        s20 = MomentSpec(; lmax_env = [2, 0], sampled = [true, false], lmax_mark = 2,
                         nbody = 3, cutoff_pair = 3.3)
        m20 = MomentBasis(xt, s20; backend = bk)
        @test m20.marked_atoms == collect(1:8)
        n_env_ge = 0
        n_mark_ge = 0
        for sx in salcs(m20), m in sx.members, t in m.terms
            msite = findfirst(sl -> sl.factor.channel == SCEFitting.DISP, t.slots)
            marka = m.atoms[t.slots[msite].site]
            marka in 5:8 && (n_mark_ge += 1)
            for sl in t.slots
                sl.factor.channel == SCEFitting.SPIN || continue
                sl.site == t.slots[msite].site && continue      # the mark's own ê
                m.atoms[sl.site] in 5:8 && (n_env_ge += 1)
            end
        end
        @test n_env_ge == 0 && n_mark_ge > 0
        # Ge rows of the design are nonzero (the μ₀ intercept at least)
        X20 = _design_moment(m20, [e2], [copy(e2)])
        @test all(any(X20[r, :] .!= 0.0) for r = 5:8)
        # (4) labels against a hand enumeration of the four rules (one mark of rank
        # 0…lmax_mark, env ranks 1…lem, even total, ≤ lsum)
        # (encode a label as the sorted tuple of ranks, the mark's rank + 100)
        labs(spec, N) = Set(Tuple(sort([d.spin_l + (has_disp(d) ? 100 : 0) for d in l]))
                            for l in SCEFitting._moment_labels(spec, N))
        sp22 = MomentSpec(; lmax_env = [2, 2], sampled = [true, true], lmax_mark = 2,
                          cutoff_pair = 3.0)
        @test labs(sp22, 1) == Set([(100,), (102,)])
        @test labs(sp22, 2) == Set([(2, 100), (1, 101), (2, 102)])
        @test labs(sp22, 3) == Set([(1, 1, 100), (2, 2, 100), (1, 2, 101),
                                    (1, 1, 102), (2, 2, 102)])
        sp3 = MomentSpec(; lmax_env = [2, 2], sampled = [true, true], lmax_mark = 2,
                         cutoff_pair = 3.0, lsum = 3)
        @test labs(sp3, 1) == Set([(100,), (102,)])
        @test labs(sp3, 2) == Set([(2, 100), (1, 101)])
        @test labs(sp3, 3) == Set([(1, 1, 100)])
        @test_throws ArgumentError MomentSpec(; lmax_env = [2], sampled = [true],
                                              cutoff_pair = 3.0, lsum = -1)
        # (4b) `lsum` PER BODY ORDER. Oracle: the same hand enumeration, one cap at a
        # time. With `lmax_env = [2]`, `lmax_mark = 2` the labels of body order N are
        # {mark rank lm ∈ 0:2} × {non-decreasing env multiset e ∈ 1:2 of length N−1}
        # with `Σl = lm + Σe` even, so by hand
        #   N = 1: Σl ∈ {0, 2}                                              → 2 labels
        #   N = 2: Σl ∈ {2 (lm0e2), 2 (lm1e1), 4 (lm2e2)}                   → 3
        #   N = 3: Σl ∈ {2, 4, 4, 4, 6}   (env sums 2,3,4 × lm 0,1,2)       → 5, one at 6
        #   N = 4: Σl ∈ {4, 4, 6, 6, 6, 8} (env sums 3..6 × lm 0,1,2)       → 6
        #          (two at 4, three at 6, one at 8)
        # so a cap of 4 keeps 2/3/4/2 and a cap of 6 keeps 2/3/5/5. Note where the two
        # caps bite: 4 opens the 4-body sector at its floor while cutting the 3-body
        # sector's Σl = 6 content away — the starvation a single global cap cannot avoid.
        nlab(spec) = [length(SCEFitting._moment_labels(spec, N)) for N = 1:spec.nbody]
        sp4(ls) = MomentSpec(; lmax_env = [2], sampled = [true], lmax_mark = 2,
                             nbody = 4, cutoff_pair = 3.0, cutoff_star = 3.0, lsum = ls)
        @test nlab(sp4(nothing)) == [2, 3, 5, 6]
        @test nlab(sp4(4)) == [2, 3, 4, 2]
        @test nlab(sp4(6)) == [2, 3, 5, 5]
        # the mixed cap the global one cannot express: 3-body at 6, 4-body at 4
        @test nlab(sp4([3 => 6, 4 => 4])) == [2, 3, 5, 2]
        @test nlab(sp4(Dict(3 => 6, 4 => 4))) == [2, 3, 5, 2]
        # COMPOSITION (oracle: the per-order spec must carry, order by order, exactly
        # what the matching single-value spec carries — checked as label SETS, in both
        # directions, and the two orders must actually differ so this is not vacuous)
        let mixed = sp4([3 => 6, 4 => 4])
            @test labs(mixed, 3) == labs(sp4(6), 3)
            @test labs(mixed, 4) == labs(sp4(4), 4)
            @test labs(sp4(6), 4) != labs(sp4(4), 4)      # the orders really differ
            @test !isempty(labs(mixed, 4))
        end
        # a scalar broadcasts to every order — the pre-per-order behavior, unchanged
        @test sp4(4).lsum == fill(4, 4)
        @test sp4(nothing).lsum == fill(SCEFitting.LSUM_UNCAPPED, 4)
        @test labs(sp4(4), 3) == labs(sp4([1 => 4, 2 => 4, 3 => 4, 4 => 4]), 3)
        # validation is the energy side's resolver, so the refusals are shared
        @test_throws ArgumentError sp4([5 => 4])          # body order out of range
        @test_throws ArgumentError sp4([3 => 4, 3 => 4])  # duplicate body order
        @test_throws ArgumentError sp4([3 => -1])         # negative cap
        @test_throws ArgumentError sp4([0, 4, 4, 4])   # positional vector: refused
        # an order the spec never asked for has no cap to read
        @test_throws ArgumentError SCEFitting._moment_labels(sp4(4), 5)
        mls = MomentBasis(xt, MomentSpec(; lmax_env = [2, 2], sampled = [true, true],
                                         lmax_mark = 2, nbody = 2, cutoff_pair = 3.3,
                                         lsum = 2); backend = bk)
        @test all(sum(d.spin_l for d in k.decors) <= 2 for k in mls.salc_basis.keys)
        ranks(k) = sort([d.spin_l for d in k.decors])
        @test any(ranks(k) == [1, 1] for k in mls.salc_basis.keys)
        @test !any(ranks(k) == [2, 2] for k in mls.salc_basis.keys)
        # (5) isotropy = false on the 3-body STAR basis removes nothing that
        # isotropy = true keeps and adds L_S ≠ 0 star blocks
        full3 = MomentBasis(xt, MomentSpec(; lmax_env = [2, 2], sampled = [true, true],
                                           lmax_mark = 2, nbody = 3, cutoff_pair = 3.3,
                                           cutoff_star = 3.3, isotropy = false);
                            backend = bk)
        @test any(k.L_S != 0 && k.body == 3 for k in full3.salc_basis.keys)
        @test [k for k in full3.salc_basis.keys if k.L_S == 0] == mb.salc_basis.keys
        # (6) census with exact Wyckoff-derived counts: the nn Fe–Fe pair orbit
        # marks 4 atoms (4a, both ends), an Fe–Ge pair orbit marks 8, and with Ge
        # unmarked the same Fe–Ge orbit marks 4
        pm = MomentBasis(xt, MomentSpec(; lmax_env = [1, 1], sampled = [true, true],
                                        lmax_mark = 1, nbody = 2, cutoff_pair = 3.0);
                         backend = bk)
        cen = moment_resolvability(pm).census
        rec_of(o) = pm.records[findfirst(k -> (k.body, k.orbit_id) ==
                                              (o.body, o.orbit_id), pm.salc_basis.keys)]
        fefe = [c for c in cen if c.body == 2 && rec_of(c).species == (1, 1)]
        fege = [c for c in cen if c.body == 2 && rec_of(c).species == (1, 2)]
        gege = [c for c in cen if c.body == 2 && rec_of(c).species == (2, 2)]
        @test length(fefe) == 1 && all(c.n_mark_atoms == 4 for c in fefe)
        @test length(fege) == 2 && all(c.n_mark_atoms == 8 for c in fege)
        @test length(gege) == 1 && all(c.n_mark_atoms == 4 for c in gege)
        @test all(c.n_mark_atoms == 4 for c in cen if c.body == 1)
        pmf = MomentBasis(xt, MomentSpec(; lmax_env = [1, 1], sampled = [true, true],
                                         lmax_mark = 1, nbody = 2, cutoff_pair = 3.0,
                                         marked = [true, false]); backend = bk)
        cenf = moment_resolvability(pmf).census
        @test all(c.n_mark_atoms == 4 for c in cenf if c.body == 2)
        # (7) rank–nullity on the reported null space, and rtol monotonicity
        res = moment_resolvability(pm)
        @test length(res.null_combinations) == length(res.kept) - res.rank
        @test moment_resolvability(pm; rtol = 1.0).rank <= 1
        @test moment_resolvability(pm; rtol = 1e-12).rank == res.rank
        # (8) vanishing columns against the random design (none on FeGe — pinned
        # against data, not against the same code)
        cf = [_mb_unit(rng, nat) for _ = 1:40]
        Xr = _design_moment(pm, cf, cf)
        cn = [norm(@view Xr[:, j]) for j = 1:size(Xr, 2)]
        @test res.vanishing == findall(<=(1e-9 * maximum(cn)), cn)
        # (9) the mark–environment edge rule after RE-ANCHORING: a P1 chain
        # i(0) j(+d) k(−d) on a short axis L with 2d ≤ cutoff_star < L and
        # L − 2d < d. The triangle {i; j, k} is enumerated from i; the j–k edge of
        # that embedding (2d) is inside the radius but is NOT the minimum image
        # (L − 2d is), so only i may carry the mark: the (0,1,1) column is exactly
        # zero on the j and k rows. A reference with no SALC code in it.
        d, L = 1.0, 2.5
        lat = Lattice(Matrix(Diagonal([L, 10.0, 10.0])))
        xch = Crystal(lat, [0.0 d/L 1-d/L; 0.0 0.0 0.0; 0.0 0.0 0.0], [1, 1, 1], ["Fe"])
        sgch = _assemble_spacegroup(xch, [SMatrix{3,3,Float64}(I)],
                                    [SVector{3,Float64}(0, 0, 0)], "P1", 1; tol = 1e-6)
        mch = MomentBasis(xch, MomentSpec(; lmax_env = [1], sampled = [true],
                                          lmax_mark = 0, nbody = 3, cutoff_pair = 1.1,
                                          cutoff_star = 2.1); backend = _MBFixedSG(sgch))
        kch = mch.salc_basis.keys
        jt = [j for j in eachindex(kch) if kch[j].body == 3 &&
              mch.records[j].edges == (1.0, 1.0, 2.0)]
        @test !isempty(jt)
        ech = _mb_unit(rng, 3)
        Xch = _design_moment(mch, [ech], [copy(ech)])
        for j in jt
            @test Xch[2, j] == 0.0 && Xch[3, j] == 0.0 && Xch[1, j] != 0.0
        end
        cch = moment_resolvability(mch).census
        @test all(c.n_mark_atoms == 1 for c in cch
                  if c.body == 3 && any(j -> (kch[j].body, kch[j].orbit_id) ==
                                             (c.body, c.orbit_id), jt))
    end
    @testset "N = 4 pointed stars" begin
        # A P1 cell with a trivial site stabilizer, so the Reynolds projector acts on a
        # ONE-dimensional space per block and the absolute constant closes by hand. The
        # crystal and the whole `MomentSpec` are load-bearing: they are what makes
        # `D = 1`, and the constant below is wrong for any fixture whose stabilizer is
        # not trivial (the FeGe star at line 200 carries an extra 1/√3 for exactly that
        # reason).
        L4 = 12.0
        cr4 = Crystal(Lattice(Matrix(L4 * I(3))),
                      [0.0 0.20 0.0 0.0; 0.0 0.0 0.24 0.0; 0.0 0.0 0.0 0.28],
                      [1, 2, 2, 2], ["Fe", "X"])
        sg4 = _assemble_spacegroup(cr4, [SMatrix{3,3,Float64}(Matrix(1.0I(3)))],
                                   [SVector{3,Float64}(0, 0, 0)], "P1", 1; tol = 1e-5)
        spec4 = MomentSpec(; lmax_env = [0, 2], sampled = [true, true], lmax_mark = 0,
                           marked = [true, false], nbody = 4, cutoff_pair = 4.0,
                           cutoff_star = 4.0, lsum = 4, isotropy = true)
        mb4 = MomentBasis(cr4, spec4; backend = _MBFixedSG(sg4))
        keys4 = mb4.salc_basis.keys
        b4 = findall(k -> k.body == 4, keys4)

        # -- the labels. `Σl = 2⌈(N−1)/2⌉` puts the 4-body sector at Σl = 4, and with
        #    `lmax_mark = 0` / `lmax_env = 2` / `lsum = 4` exactly one label survives:
        #    mark l = 0 with environment (1, 1, 2). Its naive partner, mark l = 0 with
        #    (1, 1, 1), has Σl = 3 and dies on the time-reversal screen — the unique
        #    L_S = 0 invariant of three vectors is the pseudoscalar triple product.
        @test !isempty(b4)
        @test all(j -> sort([d.spin_l for d in keys4[j].decors]) == [0, 1, 1, 2], b4)
        @test all(j -> sum(d.spin_l for d in keys4[j].decors) == 4, b4)
        # one star, three blocks: which environment site carries the l = 2 factor
        @test length(b4) == 3
        @test length(unique(keys4[j].orbit_id for j in b4)) == 1
        @test sort([keys4[j].block for j in b4]) == [1, 2, 3]
        # the N! ordering expansion folds to ONE member carrying ONE term
        for j in b4
            @test length(mb4.salc_basis.salcs[j].members) == 1
            @test length(mb4.salc_basis.salcs[j].members[1].terms) == 1
        end

        # -- the absolute normalization, derived rather than captured.
        #    Geometry: the unique L_S = 0 invariant of ranks (1, 1, 2) is
        #      êⱼ·Q(ê_l)·ê_k = (êⱼ·ê_l)(ê_k·ê_l) − (1/3)(êⱼ·ê_k),   Q(ê) = êêᵀ − I/3
        #    (the same invariant the theory page publishes for the star (2, 1, 1)).
        #    Constant: N! from the ordering convention, times 1/√D = 1 here, times
        #      κ = (4π)^{n_spin/2} · (1/√5) · (3/4π) · √(15/8π) = 3√(3/2)
        #    with n_spin = 3 (the rank-0 mark contributes no spin factor, so the scale
        #    is NOT (4π)^{N/2}), 1/√5 the unit-Frobenius normalization of the (1,1,2)→0
        #    tensor (‖T‖² = Σ_r ‖Q_r‖_F² = 5), and the two tesseral constants of
        #    Z_{1m} and Z_{2m}. A uniform loss of orderings — emitting 12 of the 24, say
        #    — leaves every RATIO unchanged and moves this constant, which is why the
        #    gate is absolute and not a ratio.
        C4 = 24 * (4π)^(3 / 2) * (1 / sqrt(5)) * (3 / (4π)) * sqrt(15 / (8π))
        @test C4 ≈ 24 * 3 * sqrt(3 / 2) rtol = 1e-14
        inv112(e, j, k, l) = dot(e[:, j], e[:, l]) * dot(e[:, k], e[:, l]) -
                             dot(e[:, j], e[:, k]) / 3
        rng4 = MersenneTwister(20260824)
        for _ = 1:3
            e4 = _mb_unit(rng4, 4)
            X4 = _design_moment(mb4, [e4], [e4])
            ref = [C4 * inv112(e4, setdiff([2, 3, 4], [l])..., l) for l in (2, 3, 4)]
            got = X4[1, b4]
            # Both the block index and the column SIGN are gauge (`_sign_canon!` is
            # allowed to flip a column; test_normalization.jl puts sign out of scope
            # for the absolute oracles for the same reason), so the gauge-free
            # statement is the multiset of MAGNITUDES. What it pins is the constant:
            # a uniform loss of orderings would move every magnitude.
            @test sort(abs.(got)) ≈ sort(abs.(ref)) rtol = 1e-12
            # ...and the three blocks are the three assignments, not three copies of
            # one: their magnitudes are distinct on a generic configuration
            @test length(unique(round.(abs.(got); digits = 8))) == 3
        end

        # -- covariance under an ARBITRARY rotation, not just a space-group operation.
        #    The mark axis has to turn with the spins; rotating only the spins leaves it
        #    behind, and an L_S = 0 column would then move.
        for _ = 1:2
            e4 = _mb_unit(rng4, 4)
            X4 = _design_moment(mb4, [e4], [e4])
            q = qr(randn(rng4, 3, 3))
            R = Matrix(q.Q) * (det(Matrix(q.Q)) < 0 ? Diagonal([-1.0, 1, 1]) : I)
            @test _design_moment(mb4, [R * e4], [R * e4]) ≈ X4 rtol = 1e-12
            # this fixture's mark has rank 0, so its ê factor is the constant |u|²R₀₀
            # and the evaluation axis is never read — the axes argument is inert here,
            # which is why the "rotate the spins but not the axes" control lives on a
            # rank-1 mark instead (the covariance testset above)
            @test _design_moment(mb4, [e4], [_mb_unit(rng4, 4)]) == X4
            # time reversal is bitwise: every label has even total spin rank
            @test _design_moment(mb4, [-e4], [-e4]) == X4
        end

        # -- a requested body order that cannot be reached is LOUD, and `show`
        #    reports what was built rather than what was asked for. The sector's
        #    `Σl` floor is 4, so an `lsum` below it drops the whole 4-body sector
        #    while every cutoff stays generous — the silent-truncation shape.
        spec_lo = MomentSpec(; lmax_env = [0, 2], sampled = [true, true],
                             lmax_mark = 0, marked = [true, false], nbody = 4,
                             cutoff_pair = 4.0, cutoff_star = 4.0, lsum = 2,
                             isotropy = true)
        mb_lo = @test_logs (:warn, r"body order 4 contributes no SALC") match_mode =
            :any MomentBasis(cr4, spec_lo; backend = _MBFixedSG(sg4))
        @test !any(k -> k.body == 4, mb_lo.salc_basis.keys)
        @test occursin("of 4 requested", sprint(show, mb_lo))
        @test !occursin("requested", sprint(show, mb4))    # nothing to disclose

        # -- opening the door adds columns rather than replacing them
        spec3 = MomentSpec(; lmax_env = [0, 2], sampled = [true, true], lmax_mark = 0,
                           marked = [true, false], nbody = 3, cutoff_pair = 4.0,
                           cutoff_star = 4.0, lsum = 4, isotropy = true)
        mb3 = MomentBasis(cr4, spec3; backend = _MBFixedSG(sg4))
        @test n_salcs(mb3) == n_salcs(mb4) - length(b4)
        @test [k for k in keys4 if k.body <= 3] == mb3.salc_basis.keys
    end

    @testset "resolvability gate streams per marked atom (27-atom supercell)" begin
        # The gate used to assemble the signature expansion as ONE dense block
        # (rows = every distinct signature key on the cell) and was killed by
        # memory on a 1296-atom torus. It now folds one block per marked atom
        # into a triangular R by stacked QR; the ORACLE is unchanged — the
        # symbolic rank must equal the rank of an independent random design
        # and every reported null combination must annihilate it — but here
        # the blocks are many (27) and the signature rows far outnumber the
        # columns, which the 8-atom FeGe fixture above does not exercise.
        L = 3
        fr = hcat([[i, j, k] ./ L for i = 0:L-1, j = 0:L-1, k = 0:L-1]...)
        csc = Crystal(Lattice(Matrix(Float64(L) * I(3))), fr, fill(1, L^3), ["Fe"])
        perms = ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
        Wc = SMatrix{3,3,Float64}[]
        for p in perms, sx in (1, -1), sy in (1, -1), sz in (1, -1)
            W = zeros(3, 3)
            W[1, p[1]] = sx; W[2, p[2]] = sy; W[3, p[3]] = sz
            push!(Wc, SMatrix{3,3,Float64}(W))
        end
        tc = [SVector{3,Float64}(i / L, j / L, k / L) for i = 0:L-1, j = 0:L-1, k = 0:L-1]
        Wall = [W for W in Wc for _ in tc]
        tall = [t for _ in Wc for t in vec(tc)]
        sgc = _assemble_spacegroup(csc, Wall, tall, "Pm-3m (3x3x3)", 221; tol = 1e-5)
        @test n_ops(sgc) == 48 * L^3
        spc = MomentSpec(; lmax_env = [2], sampled = [true], lmax_mark = 2,
                         nbody = 2, cutoff_pair = 1.8, isotropy = true)   # 3 shells
        mbc = MomentBasis(csc, spc; backend = _MBFixedSG(sgc))
        resc = moment_resolvability(mbc)
        cfgc = [_mb_unit(rng, L^3) for _ = 1:40]
        Xc = _design_moment(mbc, cfgc, cfgc)
        svc = svd(Xc).S
        @test count(>(1e-9 * svc[1]), svc) == resc.rank
        @test resc.rank == length(svc) || svc[resc.rank] / svc[resc.rank + 1] > 1e3
        @test length(resc.null_combinations) == length(resc.kept) - resc.rank
        for comb in resc.null_combinations
            v = zeros(n_salcs(mbc))
            for (j, w) in comb
                v[j] = w
            end
            @test norm(Xc * v) < 1e-10 * norm(Xc) * norm(v)
        end
        @test all(c -> c.n_mark_atoms >= 2, resc.census)

        # the same gate at body order 4. `cutoff_star = 1.1` admits the six nearest
        # neighbours only, so the star count is C(6,3) = 20 per marked atom and the
        # member count is 27 · 20 · 4! — the N! ordering expansion, checked against the
        # closed form rather than assumed.
        sp4 = MomentSpec(; lmax_env = [2], sampled = [true], lmax_mark = 2, nbody = 4,
                         cutoff_pair = 1.8, cutoff_star = 1.1, lsum = 4,
                         isotropy = true)
        mb4c = MomentBasis(csc, sp4; backend = _MBFixedSG(sgc))
        nl4 = SCEFitting.build_neighbor_list(csc, SCEFitting._star_cutoff(sp4, 4),
                                             SCEFitting.MinimumImage(); tol = 1e-8)
        @test length(SCEFitting._pointed_star_candidates(csc, nl4, sp4, 4)) ==
              L^3 * binomial(6, 3) * factorial(4)
        @test any(k -> k.body == 4, mb4c.salc_basis.keys)
        res4 = moment_resolvability(mb4c)
        X4c = _design_moment(mb4c, cfgc, cfgc)
        sv4 = svd(X4c).S
        @test count(>(1e-9 * sv4[1]), sv4) == res4.rank
        @test res4.rank == length(sv4) || sv4[res4.rank] / sv4[res4.rank + 1] > 1e3
        @test isempty(res4.vanishing)        # this cell resolves the 4-body sector
        @test res4.rank == length(res4.kept)  # ...and resolves it FULLY
        @test isempty(res4.null_combinations)
        # the mark→term index path is value-identical to the full per-SALC evaluation
        @test _design_moment(mb4c, cfgc[1:2], cfgc[1:2]) ==
              _design_moment(mb4c, cfgc[1:2], cfgc[1:2]; member_index = false)

        # Per-star-order radii COMPOSE: a spec with `cutoff_star = [r3, r4]` must carry
        # exactly the 3-body content of a single-radius `r3` spec and exactly the
        # 4-body content of a single-radius `r4` spec, with no cross-talk. The oracle
        # is that composition property, stated from what a per-order cut MEANS — not a
        # captured column count. (`r3 = 1.5` reaches the 6 + 12 = 18 first two shells,
        # `r4 = 1.1` only the 6 nearest neighbours, so the two orders really do see
        # different neighbourhoods and a shared radius could not pass this.)
        r3, r4 = 1.5, 1.1
        _spsc(nb, cut) = MomentSpec(; lmax_env = [2], sampled = [true], lmax_mark = 2,
                                    nbody = nb, cutoff_pair = 1.8, cutoff_star = cut,
                                    lsum = 4, isotropy = true)
        mb_mix = MomentBasis(csc, _spsc(4, [r3, r4]); backend = _MBFixedSG(sgc))
        mb_r3 = MomentBasis(csc, _spsc(3, r3); backend = _MBFixedSG(sgc))
        mb_r4 = MomentBasis(csc, _spsc(4, r4); backend = _MBFixedSG(sgc))
        _bodykeys(mb, b) = [k for k in mb.salc_basis.keys if k.body == b]
        for b = 1:3
            @test _bodykeys(mb_mix, b) == _bodykeys(mb_r3, b)
        end
        @test _bodykeys(mb_mix, 4) == _bodykeys(mb_r4, 4)
        # ...and the cut is real at both ends: the wide order keeps more than the
        # narrow one would, the narrow order fewer than the wide one would.
        @test length(_bodykeys(mb_mix, 3)) > length(_bodykeys(mb_r4, 3))
        @test length(_bodykeys(mb_mix, 4)) <
              length(_bodykeys(MomentBasis(csc, _spsc(4, r3);
                                           backend = _MBFixedSG(sgc)), 4))
    end
end
