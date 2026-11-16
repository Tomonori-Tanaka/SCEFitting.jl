# The gauge inside a block (2026-11-16): the columns of a block are chosen along axes fixed
# to the crystal (ported from SLCE.jl a248d4b, itself from SphericalTensorFC c7b9957),
# and the per-orbit reduction decides what a reference cell drops block by block,
# against the block's gross mass, instead of column by column against each column's own
# norm. The pointed basis re-gauges its vanishing subspace into columns of its own and
# freezes them, as upstream does.
#
# Every oracle here is outside the construction: rotated columns are compared through the
# EVALUATOR on rotated configurations; the span a cell can see is the rank of a design
# SAMPLED through `evaluate_salc` on random spins, never the symbolic aggregation the
# reduction uses; fitted models are compared through readouts that do not depend on the
# gauge (`bilinear_terms` and the energy of every cluster instance), mapped between
# descriptions by Cartesian position. The groups are found by brute force (no Spglib in
# the core suite).
#
# Mutations run by hand on 2026-11-16 (each restored afterwards; the counts are in the
# CHANGELOG entry):
#   - the column-wise reduction restored (each column judged against its own norm, in
#     the crystal frame): the margin gate and the dropped-by-name gate fail;
#   - the block's gross mass replaced by the columns' own aggregated norms: the
#     cancellation of a whole block is no longer seen (CsCl tie shell, test_salc.jl,
#     and the span gates here);
#   - the crystal frame replaced by the Cartesian axes (`_build_frame_cache` fed I): the
#     column-rotation, description and orientation gates fail (also asserted below as
#     an in-suite mutation);
#   - the forbidden band moved (cut or top): the margin gate fails;
#   - the block weight of `penalty_metric` reverted to per column, and the reference
#     ensemble drawn in the lab frame: the metric gates fail.

using SCEFitting: _assemble_spacegroup, regauged_columns, salcs, evaluate_salc,
    symmetrized_lattice, SALC, build_clusters, build_neighbor_list
using Base.CoreLogging: with_logger, NullLogger, SimpleLogger
using LinearAlgebra
using Random
using StaticArrays

_gauge_rot(axis, θ) = (k = normalize(Float64.(axis));
                       K = [0 -k[3] k[2]; k[3] 0 -k[1]; -k[2] k[1] 0];
                       Matrix(I + sin(θ) * K + (1 - cos(θ)) * K * K))

# A supercell `M` (columns: its vectors in units of the primitive ones) of the primitive
# lattice `A0` with fractional basis `F`, turned by `R`. The atom order depends on `M`
# only, so the same `M` turned by any `R` lists the same atoms in the same order.
function _gauge_supercell(A0, F, M, R)
    Minv = inv(Float64.(M))
    pos = Vector{Vector{Float64}}()
    for i = -4:4, j = -4:4, k = -4:4, s = 1:size(F, 2)
        f = Minv * (F[:, s] .+ [i, j, k])
        f = f .- floor.(f .+ 1e-9)
        any(p -> norm(p .- f) < 1e-6, pos) || push!(pos, f)
    end
    frac = reduce(hcat, pos)
    return Crystal(Lattice(R * A0 * M), frac, fill(1, size(frac, 2)), ["X"]), frac
end

# The space group by brute force: the integer matrices with entries in -1:1 that keep the
# PRIMITIVE metric, conjugated into the cell, each with every translation that maps the
# atom set onto itself (that covers the nonsymmorphic hcp operations).
struct _GaugeOps <: SCEFitting.AbstractSymmetryBackend
    rots::Vector{SMatrix{3,3,Float64,9}}
    trans::Vector{SVector{3,Float64}}
end
SCEFitting.analyze_symmetry(b::_GaugeOps, cr::Crystal; tol::Real = 1e-5) =
    _assemble_spacegroup(cr, b.rots, b.trans, "brute force", 0; tol = tol)

function _gauge_ops(A0, M, frac)
    G0 = A0' * A0
    Mf = Float64.(M)
    on_atom(x) = any(p -> maximum(abs.(mod.(x .- p .+ 0.5, 1.0) .- 0.5)) < 1e-8,
                     eachcol(frac))
    rots, trans = SMatrix{3,3,Float64,9}[], SVector{3,Float64}[]
    for t in Iterators.product(ntuple(_ -> -1:1, 9)...)
        W0 = reshape(collect(Float64, t), 3, 3)
        abs(abs(det(W0)) - 1) < 1e-9 || continue
        maximum(abs, W0' * G0 * W0 - G0) < 1e-9 * maximum(abs, G0) || continue
        W = round.(Mf \ (W0 * Mf))
        for b in axes(frac, 2)
            tr = mod.(frac[:, b] .- W * frac[:, 1], 1.0)
            all(q -> on_atom(W * q .+ tr), eachcol(frac)) || continue
            push!(rots, SMatrix{3,3,Float64}(W))
            push!(trans, SVector{3,Float64}(tr))
        end
    end
    return _GaugeOps(rots, trans)
end

const _GAUGE_HCP_A0 = [2.5 -1.25 0; 0 2.5 * sqrt(3) / 2 0; 0 0 4.08]
const _GAUGE_HCP_F = [1/3 2/3; 2/3 1/3; 1/4 3/4]
function _gauge_rh_a0(α = deg2rad(80.0), a = 3.0)
    c = cos(α); s = sin(α); cy = (c - c^2) / s
    return a * [1 c c; 0 s cy; 0 0 sqrt(1 - c^2 - cy^2)]
end
const _GAUGE_FCC_A0 = 1.8 * [0 1 1; 1 0 1; 1 1 0]
const _GAUGE_ORIGIN = reshape([0.0, 0, 0], 3, 1)
const _GAUGE_M222 = [2 0 0; 0 2 0; 0 0 2]
const _GAUGE_M_TILTED = [2 0 2; 0 2 0; 0 0 2]          # lattice vectors (a1, a2, a1 + a3)

_gquiet(f) = with_logger(f, NullLogger())

# (crystal, basis) of one description of a crystal.
function _gauge_case(A0, F, M, R, spec)
    cr, frac = _gauge_supercell(A0, F, M, R)
    b = _gquiet(() -> SCEBasis(cr, spec; backend = _gauge_ops(A0, M, frac), tol = 1e-6))
    return cr, b
end

# Pure-spin pair bases on (2,2,2) cells with Wigner–Seitz boundary ties. Measured
# (2026-11-16): hcp drops 12 directions, each a single crystal-frame column; fcc
# primitive drops 3, also single columns; the rhombohedral cell drops 7 directions out
# of three blocks whose survivors are not single columns (the kept columns recombine
# each block).
const _GAUGE_CASES = (
    ("hcp (2,2,2)", _GAUGE_HCP_A0, _GAUGE_HCP_F,
     BasisSpec(; nbody = 2, cutoff = 4.2, lmax = [2], isotropy = false)),
    ("fcc primitive (2,2,2)", _GAUGE_FCC_A0, _GAUGE_ORIGIN,
     BasisSpec(; nbody = 2, cutoff = Inf, lmax = [2], isotropy = false)),
    ("rhombohedral (2,2,2)", _gauge_rh_a0(), _GAUGE_ORIGIN,
     BasisSpec(; nbody = 2, cutoff = Inf, lmax = [2], isotropy = false)))

# Cartesian site of atom `a` in cell `R` of `cr`, and the inverse map.
_gauge_cart(cr, a, R) = cr.lattice.vectors * (cr.frac_positions[:, a] .+ R)
function _gauge_locate(cr, x)
    A = cr.lattice.vectors
    f = A \ x
    for a in axes(cr.frac_positions, 2)
        d = f .- cr.frac_positions[:, a]
        R = round.(Int, d)
        norm(A * (d .- R)) < 1e-6 && return a, R
    end
    error("no atom at $x")
end
# perm[i]: the atom of `cb` sitting where atom i of `ca` sits after turning by R.
_gauge_atom_map(ca, cb, R) = [_gauge_locate(cb, R * _gauge_cart(ca, a, [0, 0, 0]))[1]
                              for a in axes(ca.frac_positions, 2)]

# The unreduced columns of every orbit, exactly as `build_salc_basis` hands them to the
# reduction: the construction's input, NOT an oracle — the oracles below read them only
# through the evaluator.
function _gauge_raw_orbits(b)
    cr, sg, spec = b.crystal, b.spacegroup, b.spec
    nl = build_neighbor_list(cr, SCEFitting._superset_cutoff(spec), MinimumImage())
    cs = build_clusters(cr, nl, sg; nbody = spec.nbody, selection = MinimumImage(),
                        cutoff = spec.cutoff)
    maxl = maximum(spec.lmax)
    wc = SCEFitting._build_wig_cache(sg, maxl)
    fc = SCEFitting._build_frame_cache(cr, sg, maxl)
    out = Vector{SALC}[]
    for N in sort(collect(keys(cs.by_body))), (oid, O) in enumerate(cs.by_body[N])
        push!(out, SCEFitting._orbit_salcs(cr, sg, N, oid, O, spec.lmax, spec.lsum[N],
                                           spec.isotropy, wc, fc))
    end
    return out
end

_gauge_values(cols, es) = [evaluate_salc(s, e) for e in es, s in cols]

_gauge_blocks(sal) = (g = Dict{Any,Vector{Int}}();
                      for (j, s) in enumerate(sal)
                          push!(get!(g, SCEFitting._block_key(s), Int[]), j)
                      end; collect(values(g)))

@testset "the columns of a block turn with the crystal" begin
    # Turning the crystal by any orthogonal matrix — proper or not — must turn every
    # column: Φⱼ of the turned basis at the turned configuration equals Φⱼ of the
    # original at the original one. Measured through the evaluator, so the oracle shares
    # nothing with the gauge. A spin enters through even total rank, so turning it as a
    # polar vector is the same as turning it as an axial one.
    rng = MersenneTwister(301)
    turns = (_gauge_rot([1, 2, 3], 0.3), _gauge_rot([0.3, -0.7, 0.5], 2.1),
             -_gauge_rot([1, -1, 2], 0.8))                  # the last one is improper
    for (name, A0, F, spec) in _GAUGE_CASES
        cr0, b0 = _gauge_case(A0, F, _GAUGE_M222, I(3), spec)
        nat = n_atoms(cr0)
        es = [randcfg(rng, nat) for _ = 1:6]
        axes0 = symmetrized_lattice(b0).crystal_axes
        @test axes0' * axes0 ≈ I atol = 1e-14
        v0 = _gauge_values(salcs(b0), es)
        for R in turns
            cr, b = _gauge_case(A0, F, _GAUGE_M222, R, spec)
            # the axes turn with the crystal (an improper R gives improper axes)
            @test maximum(abs, symmetrized_lattice(b).crystal_axes - R * axes0) <= 1e-13
            @test b.salc_basis.keys == b0.salc_basis.keys
            @test regauged_columns(b) == regauged_columns(b0)
            v = _gauge_values(salcs(b), [R * e for e in es])
            @test maximum(abs, v - v0) <= 1e-10 * maximum(abs, v0)
        end
    end
    # In-suite mutation: the Cartesian axes in place of the crystal frame. Built through
    # the orbit engine with the frame table fed the identity, the turned crystal's
    # columns are NOT the turned columns — the gate above has something to catch. A
    # block of several columns is needed for the gauge to matter: the first one of the
    # fcc primitive cell (its lattice is not upper triangular, so even the aligned cell
    # has a frame other than I).
    name, A0, F, spec = _GAUGE_CASES[2]
    cr0, b0 = _gauge_case(A0, F, _GAUGE_M222, I(3), spec)
    R = turns[1]
    cr, b = _gauge_case(A0, F, _GAUGE_M222, R, spec)
    nl(c) = build_neighbor_list(c, SCEFitting._superset_cutoff(spec), MinimumImage())
    raw(c, sg, frame) = (cs = build_clusters(c, nl(c), sg; nbody = 2,
                                             cutoff = spec.cutoff);
                         [SCEFitting._orbit_salcs(c, sg, 2, oid, O, [2],
                                                  SCEFitting.LSUM_UNCAPPED, false,
                                                  SCEFitting._build_wig_cache(sg, 2), frame)
                          for (oid, O) in enumerate(cs.by_body[2])])
    aligned = raw(cr0, b0.spacegroup, SCEFitting._build_frame_cache(cr0, b0.spacegroup, 2))
    proper = raw(cr, b.spacegroup, SCEFitting._build_frame_cache(cr, b.spacegroup, 2))
    mutant = raw(cr, b.spacegroup, SCEFitting._build_frame_cache(Matrix(1.0I, 3, 3), 2))
    @test any(cols -> length(cols) >= 2, _gauge_blocks(reduce(vcat, aligned)))
    nat = n_atoms(cr0)
    es = [randcfg(rng, nat) for _ = 1:4]
    gap(set) = maximum(maximum(abs, _gauge_values(s, [R * e for e in es]) -
                                    _gauge_values(a, es))
                       for (s, a) in zip(set, aligned))
    @test gap(proper) <= 1e-10
    @test gap(mutant) > 1e-3
end

@testset "isotropic columns do not see the frame" begin
    # The frame turns every site axis with the same matrices. An `Lf = 0` coupled tensor
    # (pure spin, `isotropy = true`) is a scalar, so turning all its axes together leaves
    # it unchanged: the frame cannot change an isotropic column. Checked on the coupled
    # tensors themselves for a generic frame, and on the built isotropic columns: built
    # along the Cartesian axes or along a generic frame, they are the same.
    Q = _gauge_rot([0.2, 0.9, -0.4], 1.3)
    frame = SCEFitting._build_frame_cache(Q, 4)
    checked = 0
    for ls in ([1, 1], [2, 2], [1, 1, 2], [2, 2, 2], [1, 3, 2])
        for cb in SCEFitting.coupled_bases(ls; isotropy = true)
            @test cb.Lf == 0
            v = dropdims(cb.tensor; dims = ndims(cb.tensor))
            w = v
            for i in eachindex(ls)
                w = SCEFitting.AngularMomentum.nmode_mul(w,
                        SCEFitting._wig(frame, ls[i], 1), i)
            end
            @test maximum(abs, w .- v) <= 1e-13
            checked += 1
        end
    end
    @test checked > 0
    # the decor engine's spin slots on an `L_S = 0` path (the pointed basis' screen)
    for label in ([SiteDecor(; spin = 1), SiteDecor(; spin = 1, disp = (1, 0))],
                  [SiteDecor(; spin = 2), SiteDecor(; spin = 2, disp = (1, 0))])
        slots = SCEFitting._assignment_slots(label)
        for (L_S, Lf, T) in SCEFitting._decor_coupled_bases(slots, true)
            @test L_S == 0
            for Mf = 1:(2Lf + 1)
                v = SCEFitting._mfslice(T, Mf)
                w = v
                for j in eachindex(slots)
                    slots[j].factor.channel == SCEFitting.SPIN || continue
                    w = SCEFitting.AngularMomentum.nmode_mul(w,
                            SCEFitting._wig(frame, slots[j].factor.l, 1), j)
                end
                @test maximum(abs, w .- v) <= 1e-13
            end
        end
    end
    cr, b = _gauge_case(_GAUGE_FCC_A0, _GAUGE_ORIGIN, _GAUGE_M222,
                        _gauge_rot([1, 2, 3], 0.3),
                        BasisSpec(; nbody = 2, cutoff = 2.6, lmax = [2], isotropy = true))
    sg = b.spacegroup
    cs = build_clusters(cr, build_neighbor_list(cr, 2.6), sg; nbody = 2)
    wc = SCEFitting._build_wig_cache(sg, 2)
    cartesian = SCEFitting._build_frame_cache(Matrix(1.0I, 3, 3), 2)
    generic = SCEFitting._build_frame_cache(Q, 2)
    for O in cs.by_body[2]
        a = SCEFitting._orbit_salcs(cr, sg, 2, 1, O, [2], SCEFitting.LSUM_UNCAPPED, true,
                                    wc, cartesian)
        c = SCEFitting._orbit_salcs(cr, sg, 2, 1, O, [2], SCEFitting.LSUM_UNCAPPED, true,
                                    wc, generic)
        @test length(a) == length(c) > 0
        @test all(maximum(abs, ta.folded .- tc.folded) <= 1e-14
                  for (sa, sc) in zip(a, c) for (ma, mc) in zip(sa.members, sc.members)
                  for (ta, tc) in zip(ma.terms, mc.terms))
    end
end

@testset "the kept span is the span the cell sees" begin
    # Per orbit, the unreduced columns and the kept ones are sampled through the
    # evaluator on random spins. The kept columns must span exactly what the unreduced
    # ones span (nothing the cell sees is lost), and must be independent on the cell (no
    # kept direction vanishes): their count is the rank of the unreduced design, and
    # their own design has full column rank. Measured (2026-11-16) on the three cells,
    # along the axes and turned: span residual ≤ 1.4e-15, smallest relative singular
    # value of a kept orbit design ≥ 0.57.
    rng = MersenneTwister(17)
    worst_residual = 0.0
    least_singular = Inf
    for (name, A0, F, spec) in _GAUGE_CASES, R in (I(3), _gauge_rot([1, 2, 3], 0.3))
        cr, b = _gauge_case(A0, F, _GAUGE_M222, R, spec)
        es = [randcfg(rng, n_atoms(cr)) for _ = 1:200]
        kept_by_orbit = Dict{Tuple{Int,Int},Vector{SALC}}()
        for s in salcs(b)
            push!(get!(kept_by_orbit, (s.key.body, s.key.orbit_id), SALC[]), s)
        end
        dropped = 0
        raws = _gauge_raw_orbits(b)
        for raw in raws
            isempty(raw) && continue
            orbit = (raw[1].key.body, raw[1].key.orbit_id)
            kept = get(kept_by_orbit, orbit, SALC[])
            Xr = _gauge_values(raw, es)
            sr = svdvals(Xr)
            rank_raw = count(>(1e-8 * sr[1]), sr)
            @test length(kept) == rank_raw
            dropped += length(raw) - length(kept)
            isempty(kept) && continue
            Xk = _gauge_values(kept, es)
            sk = svdvals(Xk)
            least_singular = min(least_singular, sk[end] / sk[1])
            # every unreduced column lies in the kept span
            Qk = Matrix(qr(Xk).Q)
            worst_residual = max(worst_residual,
                                 norm(Xr - Qk * (Qk' * Xr)) / norm(Xr))
        end
        @test dropped == sum(length, raws) - n_salcs(b)
        @test dropped > 0                       # the cell has something to drop
    end
    @test worst_residual <= 1e-10
    @test least_singular > 1e-4
end

@testset "the rank decision of a block sits far from its cut" begin
    # The reduction decides the rank of each block at `_CANCELLATION_RTOL` (1e-10) of the
    # block's gross mass and refuses a relative singular value between that and
    # `_REGAUGE_BAND_TOP` (1e-6); the pointed re-gauge reads the same band. Measured
    # 2026-11-16 on these cells, along the axes and turned: energy blocks cancel to at
    # most 1.6e-15 and carry content at least 1.0; pointed blocks (hcp, below) 5.3e-17
    # and 4.5e-2. The gate keeps three decades of headroom on each side of the band, so
    # moving either edge toward the measured values fails it.
    worst_zero = 0.0
    least_content = Inf
    record(r) = r <= SCEFitting._CANCELLATION_RTOL ? (worst_zero = max(worst_zero, r)) :
                                                     (least_content = min(least_content, r))
    for (name, A0, F, spec) in _GAUGE_CASES, R in (I(3), _gauge_rot([1, 2, 3], 0.3))
        cr, b = _gauge_case(A0, F, _GAUGE_M222, R, spec)
        for raw in _gauge_raw_orbits(b)
            foreach(record, last(SCEFitting._reduce_orbit_salcs(raw)))
        end
    end
    cr, frac = _gauge_supercell(_GAUGE_HCP_A0, _GAUGE_HCP_F, _GAUGE_M222, I(3))
    mspec = MomentSpec(; lmax_env = [1], sampled = [true], lmax_mark = 1, nbody = 2,
                       cutoff_pair = 4.2, isotropy = false)
    for R in (I(3), _gauge_rot([1, 2, 3], 0.3))
        c = Crystal(Lattice(R * cr.lattice.vectors), frac, fill(1, n_atoms(cr)), ["X"])
        mb = _gquiet(() -> MomentBasis(c, mspec; backend = _gauge_ops(_GAUGE_HCP_A0,
                                                                     _GAUGE_M222, frac)))
        groups = Dict{Any,Vector{Int}}()
        for (j, s) in enumerate(salcs(mb))
            push!(get!(groups, (SCEFitting._block_key(s), SCEFitting._mark_class(s)),
                       Int[]), j)
        end
        for cols in values(groups)
            S, gross = SCEFitting._moment_signature_of(salcs(mb)[cols])
            foreach(record, SCEFitting._block_split(S, gross).ratios)
        end
    end
    @test worst_zero < 1e-3 * SCEFitting._CANCELLATION_RTOL
    @test least_content > 1e3 * SCEFitting._REGAUGE_BAND_TOP
    # the gate itself: a value inside the band is refused, the edges are not
    @test_throws ErrorException SCEFitting._refuse_undecidable([1.0, 1e-8], :block)
    @test SCEFitting._refuse_undecidable([1.0, 1e-12, 0.0], :block) === nothing
    @test SCEFitting._refuse_undecidable([1.0, 1e-5], :block) === nothing
end

@testset "a block that loses nothing keeps its columns, bit for bit" begin
    # The reduction touches a block only when part of it goes. On the tied cells every
    # block that keeps all its columns is the projection's own output, bit for bit, and
    # a block whose survivors are single columns keeps exactly those columns under their
    # own keys. On a cell without a tie the reduction is the identity.
    for (name, A0, F, spec) in _GAUGE_CASES
        cr, b = _gauge_case(A0, F, _GAUGE_M222, I(3), spec)
        raw = Dict(s.key => s for s in reduce(vcat, _gauge_raw_orbits(b)))
        touched = Set(regauged_columns(b))
        recombined = 0
        for (j, s) in enumerate(salcs(b))
            r = raw[s.key]
            same = length(s.members) == length(r.members) &&
                   all(m.atoms == q.atoms && m.shifts == q.shifts &&
                       all(t.slots == u.slots && t.folded == u.folded
                           for (t, u) in zip(m.terms, q.terms))
                       for (m, q) in zip(s.members, r.members))
            j in touched || @test same
            recombined += !same
        end
        # which cells recombine (see `_GAUGE_CASES`)
        @test (recombined > 0) == (name == "rhombohedral (2,2,2)")
    end
    cr, b = _gauge_case(_GAUGE_FCC_A0, _GAUGE_ORIGIN, [3 0 0; 0 3 0; 0 0 3], I(3),
                        BasisSpec(; nbody = 2, cutoff = 2.6, lmax = [2], isotropy = false))
    @test isempty(regauged_columns(b))
    for raw in _gauge_raw_orbits(b)
        kept, drops, touched, _ = SCEFitting._reduce_orbit_salcs(raw)
        @test isempty(drops) && isempty(touched)
        @test length(kept) == length(raw) && all(kept[i] === raw[i] for i in eachindex(raw))
    end
end

@testset "dropped directions are named, and the table marks their blocks" begin
    # A single dropped column is named by its key; a block whose survivors recombine it
    # is named as a block with a count, since no single key names what went.
    for (name, A0, F, spec) in _GAUGE_CASES
        cr, frac = _gauge_supercell(A0, F, _GAUGE_M222, I(3))
        io = IOBuffer()
        b = with_logger(SimpleLogger(io, Base.CoreLogging.Warn)) do
            SCEBasis(cr, spec; backend = _gauge_ops(A0, _GAUGE_M222, frac), tol = 1e-6)
        end
        log = String(take!(io))
        @test occursin("SALC reduction: dropped", log)
        if name == "rhombohedral (2,2,2)"
            @test occursin(r"blocks 1, 2, 3: 2 combination\(s\) dropped \(2 zero, 0 " *
                           r"dependent\); the 1 kept column\(s\) recombine the block", log)
        else
            @test occursin(r"block \d+: zero", log)
            @test !occursin("recombine the block", log)
        end
        # coeftable flags exactly the re-gauged columns, and they are whole blocks
        flags = coeftable(SCEPredictor(b, 0.0, zeros(n_salcs(b)))).regauged
        @test findall(flags) == regauged_columns(b)
        for cols in _gauge_blocks(salcs(b))
            @test all(flags[cols]) || !any(flags[cols])
        end
        @test !isempty(regauged_columns(b))
    end
end

# The energy each cluster instance contributes: Σⱼ βⱼ·(member of column j) at `e`. This
# is what a supercell sees (a tile carries the members, not the orbit sums), it does not
# depend on the gauge, and it turns with the crystal. Keyed by the instance in the
# indexing of `home` (turned back by `R`), so descriptions can be compared.
function _gauge_member_energies(model, e, cr, home, R)
    out = Dict{Any,Float64}()
    for (j, s) in enumerate(salcs(model.basis))
        model.jphi[j] == 0.0 && continue
        for m in s.members
            one = SALC(s.key, s.body, s.decors, s.L_S, s.Lf, [m])
            loc = sort([_gauge_locate(home, R' * _gauge_cart(cr, m.atoms[i], m.shifts[i]))
                        for i in eachindex(m.atoms)]; by = l -> (l[1], Tuple(l[2])))
            key = ([l[1] for l in loc], [l[2] - loc[1][2] for l in loc])
            out[key] = get(out, key, 0.0) + model.jphi[j] * evaluate_salc(one, e)
        end
    end
    return out
end

# The bond and on-site matrices of `bilinear_terms`, turned back by `R` and keyed in the
# indexing of `home`, both orientations of each bond.
function _gauge_bond_table(model, cr, home, R)
    bt = _gquiet(() -> bilinear_terms(model))
    out = Dict{Any,Matrix{Float64}}()
    for ((a, b, S), M) in bt.pairs
        xa = R' * _gauge_cart(cr, a, [0, 0, 0])
        xb = R' * _gauge_cart(cr, b, S)
        a0 = first(_gauge_locate(home, xa))
        b0 = first(_gauge_locate(home, xb))
        d = round.(xb - xa; digits = 5) .+ 0.0                  # no -0.0 in a key
        Mh = R' * Matrix(M) * R
        out[(a0, b0, d)] = get(out, (a0, b0, d), zeros(3, 3)) + Mh
        out[(b0, a0, -d .+ 0.0)] = get(out, (b0, a0, -d .+ 0.0), zeros(3, 3)) + Mh'
    end
    for (a, A) in bt.onsites
        out[(first(_gauge_locate(home, R' * _gauge_cart(cr, a, [0, 0, 0]))), :onsite)] =
            R' * Matrix(A) * R
    end
    return out
end

_gauge_gap(ta, tb) = maximum(norm(get(ta, k, zeros(3, 3)) - get(tb, k, zeros(3, 3)))
                             for k in union(keys(ta), keys(tb)); init = 0.0)

@testset "one crystal, four descriptions, one fitted model" begin
    # The same crystal along the axes, turned, and with lattice vectors (a1, a2, a1+a3)
    # must give the same SALC count and re-gauged count, the same predictions, and the
    # same per-bond matrices (turned back) and per-instance energies under OLS,
    # AdaptiveRidge, metric Ridge and after a refit. Measured before the change
    # (50d2c82, Spglib): the turned hcp cell's OLS bond matrices were off by 14.2
    # relative, the turned rhombohedral cell's by 1.6e-2, the turned fcc cell's
    # AdaptiveRidge and metric-Ridge ones by 4e-5 and 9e-6. After: ≤ 2e-14 everywhere.
    #
    # Those three descriptions share their crystal frame (a1 and the a1a2 plane), so the
    # columns are the same functions in all three. A fourth, the cyclic relabelling
    # (a2, a3, a1), has another frame and other columns inside its blocks; the fitted
    # function and the part of each block the cell keeps (the orthogonal complement of
    # what vanishes, which no choice of columns moves) are the same, so an unpenalized
    # fit — OLS, and the refit of a support every description keeps — must agree there
    # too (measured ≤ 3.8e-15). This is the gate that sees the block decision itself: a
    # column-by-column reduction keeps one particular column of a block, which differs
    # between frames. AdaptiveRidge (per-column weights) and the metric (a Monte-Carlo
    # sample drawn in another frame) are not invariant to that choice and are not
    # compared there (measured 4e-4 and 6e-5).
    turn = _gauge_rot([1, 2, 3], 0.3)
    for (name, A0, F, spec) in _GAUGE_CASES
        cr0, b0 = _gauge_case(A0, F, _GAUGE_M222, I(3), spec)
        nat = n_atoms(cr0)
        rng = MersenneTwister(23)
        truth = SCEPredictor(b0, 0.0, 1e-2 .* randn(rng, n_salcs(b0)))
        cfgs = [randcfg(rng, nat) for _ = 1:(3 * n_salcs(b0) + 20)]
        E = predict_energy(truth, cfgs)
        T = predict_torque(truth, cfgs)
        probe = randcfg(rng, nat)
        probes = [randcfg(rng, nat) for _ = 1:5]
        Eprobe = predict_energy(truth, probes)
        tables = Dict{Symbol,Vector{Any}}()
        predictions = Dict{Symbol,Vector{Vector{Float64}}}()
        counts = Any[]
        cyclic = [0 0 2; 2 0 0; 0 2 0]                  # lattice vectors (a2, a3, a1)
        for (M, R) in ((_GAUGE_M222, Matrix(1.0I, 3, 3)), (_GAUGE_M222, turn),
                       (_GAUGE_M_TILTED, Matrix(1.0I, 3, 3)), (cyclic, Matrix(1.0I, 3, 3)))
            other_frame = M == cyclic
            cr, b = _gauge_case(A0, F, M, R, spec)
            perm = _gauge_atom_map(cr0, cr, R)
            place(e) = (out = similar(e); out[:, perm] = R * e; out)
            ds = _gquiet(() -> SCEDataset(b, [place(e) for e in cfgs], E,
                                          [place(t) for t in T]))
            push!(counts, (n_salcs(b), length(regauged_columns(b))))
            ests = (:OLS => OLS(),
                    :AdaptiveRidge => AdaptiveRidge(lambda = 1e-6, tol = 1e-13,
                                                    max_iter = 500),
                    :metric_Ridge => Ridge(b; lambda = 1e-3, torque_weight = 0.5))
            for (tag, est) in ests
                f = _gquiet(() -> fit(SCEFit, ds, est; torque_weight = 0.5))
                fits = [tag => f]
                tag === :AdaptiveRidge && push!(fits, :refit => _gquiet(() -> refit(f)))
                for (t2, g) in fits
                    other_frame && !(t2 in (:OLS, :refit)) && continue
                    m = SCEPredictor(g)
                    pred = predict_energy(m, [place(e) for e in probes])
                    t2 === :OLS &&
                        @test maximum(abs, pred - Eprobe) <= 1e-10 * maximum(abs, Eprobe)
                    push!(get!(predictions, t2, Vector{Float64}[]), pred)
                    push!(get!(tables, t2, Any[]),
                          (_gauge_bond_table(m, cr, cr0, R),
                           _gauge_member_energies(m, place(probe), cr, cr0, R)))
                end
            end
        end
        @test allequal(counts)
        for (tag, list) in predictions
            @test all(maximum(abs, p - list[1]) <= 1e-10 * maximum(abs, list[1])
                      for p in list[2:end])
        end
        for (tag, list) in tables
            (ba, ea) = list[1]
            bscale = maximum(norm, values(ba))
            escale = maximum(abs, values(ea))
            for (bb, eb) in list[2:end]
                gb = _gauge_gap(ba, bb)
                gb <= 1e-10 * bscale || @info "bond matrices differ" name tag gb bscale
                @test gb <= 1e-10 * bscale
                @test keys(ea) == keys(eb)
                ge = maximum(abs(ea[k] - eb[k]) for k in keys(ea))
                @test ge <= 1e-10 * escale
            end
        end
    end
end

@testset "the penalty metric weighs a block, from a sample the crystal carries" begin
    # One weight per block (the trace of the block's reference Gram over its columns),
    # and the reference ensemble drawn in the crystal frame: the turned crystal's metric
    # is the same numbers, not a second Monte-Carlo estimate. Before the change the
    # per-column weights of the turned fcc cell differed from the aligned cell's by the
    # seed noise (a few per cent), and a metric Ridge moved the bond matrices by 9e-6.
    for (name, A0, F, spec) in _GAUGE_CASES
        cr0, b0 = _gauge_case(A0, F, _GAUGE_M222, I(3), spec)
        cr, b = _gauge_case(A0, F, _GAUGE_M222, _gauge_rot([0.3, -0.7, 0.5], 0.9), spec)
        m0 = penalty_metric(b0; torque_weight = 0.3, nconfig = 256)
        m = penalty_metric(b; torque_weight = 0.3, nconfig = 256)
        @test all(>(0), m0)
        for cols in _gauge_blocks(salcs(b0))
            @test allequal(m0[cols])
        end
        @test maximum(abs, m - m0) <= 1e-12 * maximum(m0)
    end
end

@testset "persistence carries the gauge convention and the re-gauged columns" begin
    # A fresh build is convention 2 and records its re-gauged columns; both survive a
    # save/load. Documents before schema 8 were written under the Cartesian-frame,
    # column-by-column convention: they load verbatim (same tensors, same predictions)
    # but carry the legacy tag, so their fingerprint never matches a fresh build of the
    # same crystal and no fingerprint check pairs their coefficients with new columns.
    name, A0, F, spec = _GAUGE_CASES[3]
    cr, b = _gauge_case(A0, F, _GAUGE_M222, I(3), spec)
    @test b.salc_basis.gauge == SCEFitting.SALC_GAUGE_CONVENTION
    @test !isempty(regauged_columns(b))
    doc = SCEFitting._to_doc(b)
    @test doc["regauged_columns"] == regauged_columns(b)
    now = _gquiet(() -> SCEFitting._basis_from_doc(doc))
    @test now.salc_basis.fingerprint == b.salc_basis.fingerprint
    @test regauged_columns(now) == regauged_columns(b)
    old = copy(doc)
    old["schema_version"] = 7
    delete!(old, "regauged_columns")
    legacy = _gquiet(() -> SCEFitting._basis_from_doc(old))
    @test legacy.salc_basis.gauge == 1
    @test isempty(regauged_columns(legacy))
    @test legacy.salc_basis.keys == b.salc_basis.keys
    @test legacy.salc_basis.fingerprint != b.salc_basis.fingerprint
    e = randcfg(MersenneTwister(5), n_atoms(cr))
    @test [evaluate_salc(s, e) for s in salcs(legacy)] ==
          [evaluate_salc(s, e) for s in salcs(b)]
    # a schema-8 document must say which columns were re-gauged
    bad = copy(doc)
    delete!(bad, "regauged_columns")
    @test_throws ArgumentError SCEFitting._basis_from_doc(bad)
    bad["regauged_columns"] = [n_salcs(b) + 1]
    @test_throws ArgumentError SCEFitting._basis_from_doc(bad)
end

@testset "pointed basis: the vanishing subspace of a block is frozen" begin
    # hcp (2,2,2), pair stars, mark and environment l ≤ 1, all `L_S`. Each block (one
    # mark class) must vanish on the sampled pointed design exactly as far as
    # `moment_resolvability` freezes it, `regauged_columns` must name exactly the blocks
    # of which a proper part vanishes, and the columns must turn with the crystal.
    cr, frac = _gauge_supercell(_GAUGE_HCP_A0, _GAUGE_HCP_F, _GAUGE_M222, I(3))
    spec = MomentSpec(; lmax_env = [1], sampled = [true], lmax_mark = 1, nbody = 2,
                      cutoff_pair = 4.2, isotropy = false)
    build(c) = _gquiet(() -> MomentBasis(c, spec;
                                         backend = _gauge_ops(_GAUGE_HCP_A0, _GAUGE_M222,
                                                              frac)))
    mb = build(cr)
    vanishing = Set(moment_resolvability(mb).vanishing)
    @test !isempty(vanishing)
    rng = MersenneTwister(53)
    nat = n_atoms(cr)
    es = [randcfg(rng, nat) for _ = 1:30]
    X = SCEFitting._design_moment(mb, es, [copy(e) for e in es])
    ref = maximum(norm(view(X, :, j)) for j in axes(X, 2))
    groups = Dict{Any,Vector{Int}}()
    for (j, s) in enumerate(salcs(mb))
        push!(get!(groups, (SCEFitting._block_key(s), SCEFitting._mark_class(s)), Int[]), j)
    end
    flagged = Set(regauged_columns(mb))
    partial = 0
    for cols in values(groups)
        nullity = length(cols) - count(>(1e-8 * ref), svdvals(X[:, cols]))
        @test nullity == count(in(vanishing), cols)
        @test all(in(flagged), cols) == (0 < nullity < length(cols))
        partial += 0 < nullity < length(cols)
    end
    @test partial > 0                       # a block of this cell is partly vanishing
    @test all(j -> norm(view(X, :, j)) <= 1e-12 * ref, vanishing)
    # the pointed metric weighs a block, frozen columns exempt
    m = penalty_metric(mb; nconfig = 256)
    @test all(iszero, m[collect(vanishing)])
    for cols in values(groups)
        live = setdiff(cols, vanishing, SCEFitting._intercept_columns(mb))
        @test allequal(m[live])
    end
    # turned, the pointed columns are the turned columns (rows: one per marked atom)
    R = _gauge_rot([1, 2, 3], 0.3)
    crR = Crystal(Lattice(R * cr.lattice.vectors), frac, fill(1, nat), ["X"])
    mbR = build(crR)
    @test mbR.salc_basis.keys == mb.salc_basis.keys
    XR = SCEFitting._design_moment(mbR, [R * e for e in es], [R * e for e in es])
    @test maximum(abs, XR - X) <= 1e-10 * ref
    @test regauged_columns(mbR) == regauged_columns(mb)
end
