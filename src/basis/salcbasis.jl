# SALC construction (arbitrary body order). For each cluster orbit, project the
# representative's coupled coefficient tensor onto the trivial irrep of its site
# stabilizer (orbit–stabilizer theorem), gauge-fix a deterministic orthonormal
# basis, then transport to all orbit members.
#
# At N ≥ 3 a stabilizer operation can permute symmetry-equivalent sites. That
# permutation (a) mixes coupling paths sharing the same Lf and (b) — when the
# permuted sites carry unequal l — mixes l-orderings. The projection therefore
# runs over the COMBINED space of (ordering o, coupling path p, Mf), with each
# operation acting by rotating every site axis (wignerD_real(o_i, R)) and then
# permuting the site axes by π; the action matrix is read off by contracting
# against the package's own orthonormal coupled tensors (no 6j/9j). A surviving
# invariant can span several orderings, so a SALC carries one `SALCTerm` per
# ordering (`basis/salc.jl`). For 1/2-body and equal-l channels this reduces to a
# single term and a single path, matching the earlier construction.
#
# Entries of a coefficient tensor that are round-off against its largest entry are set
# to exactly zero (`_ROUNDOFF_RTOL`), and the group average is checked to be idempotent
# (eigenvalues 0 or 1 to 1e-6). Both decisions are relative, and neither depends on how
# the crystal is turned against the Cartesian axes; the invariance tests and the
# orientation tests (`test/unit/test_orientation.jl`) hold them. (Until 2026-11-15 the
# drops were absolute, 1e-10 / 1e-12, which deleted genuine entries of size θ in a
# crystal turned by a small angle θ. [Backported from SLCE.jl 2f7d2d5.])
#
# The gauge inside a block — which orthonormal columns span the invariant space of one
# `(l-multiset, Lf)` — is chosen along axes fixed to the crystal (`_crystal_frame`), so
# the columns of a turned crystal are the turned columns, and the per-orbit reduction at
# the end of this file decides what a reference cell cannot see block by block
# (2026-11-16; `test/unit/test_gauge.jl`).

"""
Sentinel meaning "no `Σl` cap" for a body order (in `BasisSpec.lsum` and
`build_salc_basis`'s `lsum_by_body`). Only ever compared against (the per-site
enumeration cap subtracts at most `nbody − 1` from it, which cannot overflow).
"""
const LSUM_UNCAPPED = typemax(Int)

# Per-site l-tuples to enumerate for an orbit, one canonical representative per
# orbit of the site-permutation group `perms` (Σl even = time-reversal even,
# Σl ≤ lsumN = the body order's angular budget). Each per-site range is also
# tightened to `lsumN − (N − 1)` (every other site carries l ≥ 1), so an
# uncapped species lmax does not blow up the product before the filter.
function _enumerate_ls(N::Int, species::Vector{Int}, lmax::Vector{Int}, lsumN::Int,
                       perms::Vector{Vector{Int}})
    cap = lsumN - (N - 1)
    ranges = ntuple(i -> 1:min(lmax[species[i]], cap), N)
    seen = Set{Vector{Int}}()
    out = Vector{Int}[]
    for tt in Iterators.product(ranges...)
        t = collect(Int, tt)
        s = sum(t)
        (iseven(s) && s <= lsumN) || continue
        t in seen && continue
        orbit = unique([t[p] for p in perms])
        for o in orbit
            push!(seen, o)
        end
        push!(out, minimum(orbit))   # lex-min canonical of the permutation orbit
    end
    return out
end

# Site permutation aligning `src` (= g·rep images) to `dst` up to a common lattice
# translation: `perm[i] = j` iff `g·(rep siteᵢ)` lands on `dst[j]`. `nothing` if no
# alignment exists.
function _align(src::Vector{Tuple{Int,SVector{3,Int}}},
                dst::Vector{Tuple{Int,SVector{3,Int}}})
    N = length(src)
    @inbounds for j = 1:N
        src[1][1] == dst[j][1] || continue
        t = src[1][2] - dst[j][2]               # translation aligning src[1] → dst[j]
        perm = zeros(Int, N)
        used = falses(N)
        ok = true
        for i = 1:N
            atom = src[i][1]
            shift = src[i][2] - t
            found = 0
            for k = 1:N
                used[k] && continue
                if dst[k][1] == atom && dst[k][2] == shift
                    found = k
                    break
                end
            end
            found == 0 && (ok = false; break)
            perm[i] = found
            used[found] = true
        end
        ok && return perm
    end
    return nothing
end

_sites(m::ClusterMember) = Tuple{Int,SVector{3,Int}}[(m.atoms[i], m.shifts[i]) for i in eachindex(m.atoms)]

_op_images(crystal::Crystal, sg::SpaceGroup, g::Int, rep::ClusterMember) =
    Tuple{Int,SVector{3,Int}}[_site_image(crystal, sg, g, rep.atoms[i], rep.shifts[i])
                             for i in eachindex(rep.atoms)]

# Site stabilizer: ops mapping the (unordered) representative to itself, each with
# its induced site permutation.
function _stabilizer(crystal::Crystal, sg::SpaceGroup, rep::ClusterMember)
    dst = _sites(rep)
    out = Tuple{Int,Vector{Int}}[]
    for g = 1:n_ops(sg)
        perm = _align(_op_images(crystal, sg, g, rep), dst)
        perm === nothing || push!(out, (g, perm))
    end
    return out
end

# Translation signature of a raw site list (same normalization as `_member_sig`).
_sig_of_sites(sites::Vector{Tuple{Int,SVector{3,Int}}}) =
    _sig_of_sites(sites, Val(length(sites)))
@inline function _sig_of_sites(sites, ::Val{N}) where {N}
    t = ntuple(k -> (sites[k][1], sites[k][2][1], sites[k][2][2], sites[k][2][3]), Val(N))
    return _normalize_sites(t)
end

# Op + induced permutation connecting the representative to every orbit member
# (`g·rep` aligns to `member`), in a single O(n_ops) sweep: apply each op to `rep`
# once, match its image (by translation signature) to the member(s) of that class
# — the orbit lists every anchor-variant, so one class can hold several — and record
# the first-ascending connecting op for each. The transported term is
# stabilizer-invariant, so which connecting op is chosen is moot.
function _connect_all(crystal::Crystal, sg::SpaceGroup, rep::ClusterMember,
                      members::Vector{ClusterMember})
    M = length(members)
    sig2members = Dict{Any,Vector{Int}}()
    for j = 1:M
        push!(get!(sig2members, _member_sig(members[j]), Int[]), j)
    end
    conns = Vector{Tuple{Int,Vector{Int}}}(undef, M)
    filled = falses(M)
    remaining = M
    for g = 1:n_ops(sg)
        remaining == 0 && break
        imgs = _op_images(crystal, sg, g, rep)
        js = get(sig2members, _sig_of_sites(imgs), nothing)
        js === nothing && continue
        for j in js
            filled[j] && continue
            perm = _align(imgs, _sites(members[j]))
            perm === nothing && continue
            conns[j] = (g, perm)
            filled[j] = true
            remaining -= 1
        end
    end
    remaining == 0 || error("could not connect all orbit members to the representative")
    return conns
end

# The Mf-th multiplet slice of a coupled tensor (rank N+1 → rank N), as a view: the
# last axis is column-major-contiguous, so this is a strided view with no copy. Every
# consumer (`nmode_mul`, `_frob`, the fold broadcast) reads it without mutating.
_mfslice(T::AbstractArray, Mf::Int) = selectdim(T, ndims(T), Mf)

_frob(A::AbstractArray, B::AbstractArray) = sum(A .* B)

# Real Wigner-D matrix `D^l(R_g)`, keyed by `(l, g)` for the whole basis build (it
# depends only on the operation and `l`, but is needed once per stabilizer column and
# once per transported member/term — recomputing it dominated the cost). The full
# `(l, g)` grid is bounded (`l ≤ lmax`, `g ≤ n_ops`) and cheap, so it is precomputed
# **serially** up front; the cache is then **read-only** during the threaded orbit
# loop (a concurrent `get!` on a shared `Dict` would race and corrupt it).
const _WigCache = Dict{Tuple{Int,Int},Matrix{Float64}}
function _build_wig_cache(sg::SpaceGroup, maxl::Int)::_WigCache
    cache = _WigCache()
    for l = 1:maxl, g = 1:n_ops(sg)
        cache[(l, g)] = AngularMomentum.wignerD_real(l, sg.ops[g].rotation_cart)
    end
    return cache
end
_wig(cache::_WigCache, l::Int, g::Int) = cache[(l, g)]   # read-only lookup

# The Wigner matrices of ONE orthogonal matrix, the crystal frame (`_crystal_frame`),
# in the same table shape as the operations' cache: entry `(l, 1)` is `D^l(Q)`. Built
# serially next to the operations' cache and read-only below, for the same reason.
function _build_frame_cache(frame::AbstractMatrix{<:Real}, maxl::Int)::_WigCache
    cache = _WigCache()
    for l = 1:maxl
        cache[(l, 1)] = AngularMomentum.wignerD_real(l, Matrix{Float64}(frame))
    end
    return cache
end
_build_frame_cache(crystal::Crystal, spacegroup::SpaceGroup, maxl::Int) =
    _build_frame_cache(_crystal_frame(crystal, spacegroup), maxl)

# Project the combined (ordering, path, Mf) space onto stabilizer invariants for a
# fixed final `Lf`, gauge-fix, and fold each invariant into per-ordering tensors.
# Returns a list of blocks; each block is `Vector{(ls, folded)}` (the rep terms).
#
# The invariants of a block span a space that the symmetry fixes; which vectors of
# that space become the columns (the gauge) it does not. `_canonical_basis` decides by
# looking at components, so the answer depends on the axes the components are taken
# along. They are taken along the axes of the crystal (`_crystal_frame`; `fcache` is
# the Wigner table of that one matrix), not along the Cartesian axes: the columns of a
# crystal turned as a whole are then the turned columns, and nothing that picks out
# single columns of a block depends on the orientation (ported from SLCE.jl a248d4b,
# itself from SphericalTensorFC c7b9957). Every site axis is a spin axis here and
# turns with the same matrix.
function _project_and_fold(stab::Vector{Tuple{Int,Vector{Int}}},
                           orderings::Vector{Vector{Int}}, cbs::Vector, Lf::Int,
                           wcache::_WigCache, fcache::_WigCache)
    N = length(orderings[1])
    # coupled tensors of this Lf, per ordering — selected from the prebuilt `cbs`
    # (built once per ordering in `_orbit_salcs`, not recomputed for each Lf).
    tens = [Array{Float64}[] for _ in orderings]
    for oi in eachindex(orderings)
        for cb in cbs[oi]
            cb.Lf == Lf && push!(tens[oi], cb.tensor)
        end
    end
    cols = Tuple{Int,Int,Int}[]            # (ordering index, path index, Mf index)
    for oi in eachindex(orderings), p in eachindex(tens[oi]), Mf = 1:(2Lf + 1)
        push!(cols, (oi, p, Mf))
    end
    D = length(cols)
    D == 0 && return Vector{Tuple{Vector{Int},Array{Float64}}}[]
    colidx = Dict(cols[k] => k for k = 1:D)

    P = zeros(Float64, D, D)
    for (g, perm) in stab
        for k = 1:D
            (oi, p, Mf) = cols[k]
            o = orderings[oi]
            v = _mfslice(tens[oi][p], Mf)
            for i = 1:N
                v = AngularMomentum.nmode_mul(v, _wig(wcache, o[i], g), i)
            end
            # U_g sends rep axis i → position π(i)=perm[i], i.e. permutedims by π⁻¹.
            q = invperm(perm)
            v = permutedims(v, q)
            oprime = o[q]
            oi2 = findfirst(==(oprime), orderings)
            oi2 === nothing && error("ordering set not closed under the stabilizer")
            for p2 in eachindex(tens[oi2]), Mf2 = 1:(2Lf + 1)
                c = _frob(_mfslice(tens[oi2][p2], Mf2), v)
                abs(c) < _ROUNDOFF_RTOL && continue    # the tensors are orthonormal
                P[colidx[(oi2, p2, Mf2)], k] += c
            end
        end
    end
    P ./= length(stab)
    P .= (P .+ P') ./ 2
    # `D` is tiny (a few dozen), so OpenBLAS runs this `eigen` serially even though
    # we are inside `Threads.@threads`; no Julia↔BLAS thread oversubscription in practice.
    F = eigen(Symmetric(P))
    all(λ -> abs(λ) < 1e-6 || abs(λ - 1) < 1e-6, F.values) ||
        error("SALC projector is not idempotent (eigenvalues $(F.values)); convention bug")
    V1 = F.vectors[:, findall(>(0.5), F.values)]
    # Turning a tensor from the axes of the crystal to the Cartesian axes is an
    # orthogonal matrix in these coordinates; the sites stay where they are.
    to_cartesian = zeros(Float64, D, D)
    for k = 1:D
        (oi, p, Mf) = cols[k]
        o = orderings[oi]
        v = _mfslice(tens[oi][p], Mf)
        for i = 1:N
            v = AngularMomentum.nmode_mul(v, _wig(fcache, o[i], 1), i)
        end
        for p2 in eachindex(tens[oi]), Mf2 = 1:(2Lf + 1)
            c = _frob(_mfslice(tens[oi][p2], Mf2), v)
            abs(c) < _ROUNDOFF_RTOL && continue
            to_cartesian[colidx[(oi, p2, Mf2)], k] = c
        end
    end
    W = to_cartesian * _canonical_basis(to_cartesian' * V1)

    blocks = Vector{Tuple{Vector{Int},Array{Float64}}}[]
    for b = 1:size(W, 2)
        c = view(W, :, b)
        terms = Tuple{Vector{Int},Array{Float64}}[]
        for (oi, o) in enumerate(orderings)
            Fo = zeros(Float64, ntuple(i -> 2o[i] + 1, N)...)
            for p in eachindex(tens[oi]), Mf = 1:(2Lf + 1)
                coef = c[colidx[(oi, p, Mf)]]
                coef == 0.0 && continue
                Fo .+= coef .* _mfslice(tens[oi][p], Mf)
            end
            _zero_roundoff!(Fo)
            norm(Fo) > _ROUNDOFF_RTOL && push!(terms, (collect(o), Fo))
        end
        push!(blocks, terms)
    end
    return blocks
end

"""
Relative size below which an entry of a coefficient tensor is round-off: an entry
smaller than this times the largest entry of its tensor is set to exactly zero
(`_zero_roundoff!`).

Exact zeros matter downstream — a term whose tensor is zero is dropped, and so is a
member without terms. The cut must sit at round-off and be relative. A crystal turned
by a small angle `θ` against the Cartesian axes has genuine entries of size `θ`, and
an absolute cut of `1e-10`, which stood here before 2026-11-15, kept them for `θ`
above it and deleted them below: the basis then jumped with the orientation (a
tetragonal pair basis emitted 5–6 SALCs instead of 4 for θ between 1e-8 and 1e-5).
[Backported from SLCE.jl 2f7d2d5, itself ported from SphericalTensorFC 311e08f.]
"""
const _ROUNDOFF_RTOL = 1e-13

# Set the entries of `tensor` that are round-off against its largest entry to zero.
function _zero_roundoff!(tensor::AbstractArray{Float64})
    cut = _ROUNDOFF_RTOL * maximum(abs, tensor; init = 0.0)
    @inbounds for index in eachindex(tensor)
        abs(tensor[index]) < cut && (tensor[index] = 0.0)
    end
    return tensor
end

# Deterministic, BLAS-independent gauge: axis-pivoted modified Gram–Schmidt on the
# block projector, with a first-significant-component sign fix.
#
# A column of the projector is taken as the next basis vector when what is left of it
# after removing the vectors already chosen has norm above `_GAUGE_PIVOT_MIN`. The
# threshold is that high on purpose: normalizing a remainder of norm `ν` multiplies
# its round-off by `1/ν`, and in a crystal turned by a small angle `θ` remainders of
# norm about `θ` occur (at 1e-8, which stood here before 2026-11-15, the amplified
# round-off changed the SALC count for θ between 1e-8 and 1e-5). It can never leave
# the basis short: the squared remainders of all columns add up to the number of
# vectors still missing, so while one is missing some remainder exceeds `1/√n` for an
# `n`-dimensional block. [Backported from SLCE.jl 2f7d2d5.]
const _GAUGE_PIVOT_MIN = 1e-4

_canonical_basis(V1::AbstractMatrix{Float64}) = first(_canonical_basis_pivots(V1))

# `_canonical_basis` and, for each of its vectors, the coordinate axis it was pivoted
# on (`W[:, k]` comes from column `pivots[k]` of the projector). The per-orbit
# reduction keys a recombined column by its pivot, so a subspace spanned by single
# coordinate axes comes back as exactly those axes under their own keys.
function _canonical_basis_pivots(V1::AbstractMatrix{Float64})
    n, d = size(V1)
    d == 0 && return zeros(Float64, n, 0), Int[]
    Q = V1 * V1'
    W = zeros(Float64, n, d)
    pivots = Int[]
    k = 0
    for j = 1:n
        u = Q[:, j]
        for i = 1:k
            u .-= dot(view(W, :, i), u) .* view(W, :, i)
        end
        nu = norm(u)
        if nu > _GAUGE_PIVOT_MIN
            k += 1
            W[:, k] = u ./ nu
            _sign_canon!(view(W, :, k))
            push!(pivots, j)
        end
        k == d && break
    end
    k == d || error("canonical gauge failed to span the invariant subspace")
    return W, pivots
end

function _sign_canon!(v)
    for x in v
        if abs(x) > 1e-8
            x < 0 && (v .*= -1)
            return
        end
    end
end

# Transport a rep term `(o, F)` to a member connected by `(g, perm)`: rotate each
# site axis by `wignerD_real(o_i, R_g)`, then relabel rep site i → member site
# perm[i]. Returns the member-ordered `(ls, folded)`.
function _transport_term(o::Vector{Int}, F::Array{Float64}, g::Int,
                         perm::Vector{Int}, wcache::_WigCache)
    N = length(o)
    T = F
    for i = 1:N
        T = AngularMomentum.nmode_mul(T, _wig(wcache, o[i], g), i)
    end
    q = invperm(perm)                       # member axis j ← rep axis q[j]
    G = _zero_roundoff!(Array{Float64}(permutedims(T, q)))
    return o[q], G                          # member ls = o[invperm(perm)]
end

# ---------------------------------------------------------------------------
# Mixed-channel (decor) projection engine.
#
# The pure-spin `_project_and_fold`/`_transport_term` above stay the production
# path for spin-only bases (pin- and oracle-checked bitwise); the `_decors`
# engine below generalizes the same construction to per-site `SiteDecor`
# assignments, which is what the pointed site-moment channel needs (its mark is
# a displacement decor). The "engines agree on pure spin" testset in
# `test/unit/test_mixedsalc.jl` is the anti-drift gate coupling the two — change
# one and re-check the other.
#
# Slot conventions (see `basis/salc.jl`): an assignment `a::Vector{SiteDecor}`
# (per site) realizes the canonical slot list SPIN axes before DISP axes, each
# by site order. Coupling runs over the slot `l`s in slot order (spin first),
# so the total spin rank `L_S` is the running coupled momentum after the last
# spin slot — a good quantum number (site permutations act within channels and
# commute with the diagonal rotation), read off each `coupling_paths` entry by
# `_path_LS`. Both channels rotate through the SAME polar Wigner cache: the
# Σl_spin-even screen makes det(R)^{Σl_spin} ≡ +1, so the axial spin action
# equals the polar one; displacement radial factors |u|^{2k} are
# rotation-invariant passive labels.
# ---------------------------------------------------------------------------

# Canonical slot list of a per-site decor assignment.
function _assignment_slots(a::Vector{SiteDecor})::Vector{Slot}
    slots = Slot[]
    for s in eachindex(a)
        has_spin(a[s]) && push!(slots, Slot(s, SiteFactor(SPIN, 0, a[s].spin_l)))
    end
    for s in eachindex(a)
        has_disp(a[s]) &&
            push!(slots, Slot(s, SiteFactor(DISP, a[s].disp_k, a[s].disp_l)))
    end
    return slots
end

_slot_ls(slots::Vector{Slot})::Vector{Int} = Int[sl.factor.l for sl in slots]
_n_spin_slots(slots::Vector{Slot})::Int = count(sl -> sl.factor.channel == SPIN, slots)

# Total spin rank of a left-coupling path over a spin-first slot list: the
# running coupled momentum after the last spin slot (0 / l₁ / Lseq / Lf edges).
function _path_LS(ls::Vector{Int}, Lseq::Vector{Int}, Lf::Int, n_spin::Int)::Int
    n_spin == 0 && return 0
    n_spin == 1 && return ls[1]
    n_spin == length(ls) && return Lf
    return Lseq[n_spin - 1]
end

# All distinct arrangements of a decor multiset over the cluster sites (swap
# recursion + dedup; orbit body orders are tiny). Sorted for determinism.
function _multiset_arrangements(label::Vector{SiteDecor})::Vector{Vector{SiteDecor}}
    n = length(label)
    out = Set{Vector{SiteDecor}}()
    perm = collect(1:n)
    function rec(k::Int)
        if k >= n
            push!(out, label[perm])
            return
        end
        for i = k:n
            perm[k], perm[i] = perm[i], perm[k]
            rec(k + 1)
            perm[k], perm[i] = perm[i], perm[k]
        end
    end
    n == 0 && return [SiteDecor[]]
    rec(1)
    return sort!(collect(out))
end

# Image of an assignment under a stabilizer site permutation (site s → perm[s]).
function _assignment_image(a::Vector{SiteDecor}, perm::Vector{Int})::Vector{SiteDecor}
    a2 = Vector{SiteDecor}(undef, length(a))
    for s in eachindex(a)
        a2[perm[s]] = a[s]
    end
    return a2
end

# Axis permutation induced on the canonical slot lists by a site permutation:
# σ[j] = position of the image of slot j (site perm[s], same factor) in the
# image assignment's canonical slot list. Unique by the (site, channel) slot
# invariant.
function _slot_sigma(slots::Vector{Slot}, slots2::Vector{Slot},
                     perm::Vector{Int})::Vector{Int}
    σ = Vector{Int}(undef, length(slots))
    for j in eachindex(slots)
        target = Slot(perm[slots[j].site], slots[j].factor)
        jp = findfirst(==(target), slots2)
        jp === nothing && error("slot lists not closed under the site permutation")
        σ[j] = jp
    end
    return σ
end

# Coupled bases of one assignment, tagged with L_S: (L_S, Lf, tensor) triples
# over the slot-order `l`s (spin first ⇒ L_S well-defined per path).
#
# `isotropy` is the decor engine's spelling of the pure-spin engine's screen:
# there it keeps `Lf == 0`, here `L_S == 0`. On a pure-spin label the two agree
# identically (every slot is a spin slot, so `_path_LS` returns `Lf`), which is
# what the cross-engine gate checks. On a mixed label they do NOT agree and
# neither implies the other (see `build_real_bases`' `keep`), so the screen is
# handed down as a path predicate and rejected paths never build a tensor.
function _decor_coupled_bases(slots::Vector{Slot}, isotropy::Bool)
    ls = _slot_ls(slots)
    n_spin = _n_spin_slots(slots)
    out = Tuple{Int,Int,Array{Float64}}[]
    keep = (Lseq, Lf) -> !isotropy || _path_LS(ls, Lseq, Lf, n_spin) == 0
    for (Lseq, Lf, tensor) in AngularMomentum.build_real_bases(ls; keep = keep)
        push!(out, (_path_LS(ls, Lseq, Lf, n_spin), Lf, tensor))
    end
    return out
end

# Decor-general `_project_and_fold`: project the combined (assignment, path, Mf)
# space onto stabilizer invariants for a fixed (L_S, Lf) block. Returns a list
# of invariant blocks; each is `Vector{(assignment index, folded)}`.
#
# The gauge inside a block is chosen along the crystal axes exactly as in
# `_project_and_fold`: every slot, spin and displacement alike, turns with the frame
# `fcache`. On an `L_S = 0` path the spin slots couple to a scalar, which no rotation
# moves, so including them changes nothing there (gated in `test_gauge.jl`).
function _project_and_fold_decors(stab::Vector{Tuple{Int,Vector{Int}}},
                                  assignments::Vector{Vector{SiteDecor}},
                                  slotlists::Vector{Vector{Slot}},
                                  cbs::Vector, L_S::Int, Lf::Int,
                                  wcache::_WigCache, fcache::_WigCache)
    tens = [Array{Float64}[] for _ in assignments]
    for ai in eachindex(assignments)
        for (ls_path, lf_path, tensor) in cbs[ai]
            (ls_path == L_S && lf_path == Lf) && push!(tens[ai], tensor)
        end
    end
    cols = Tuple{Int,Int,Int}[]            # (assignment index, path index, Mf)
    for ai in eachindex(assignments), p in eachindex(tens[ai]), Mf = 1:(2Lf + 1)
        push!(cols, (ai, p, Mf))
    end
    D = length(cols)
    D == 0 && return Vector{Tuple{Int,Array{Float64}}}[]
    colidx = Dict(cols[k] => k for k = 1:D)

    P = zeros(Float64, D, D)
    for (g, perm) in stab
        for k = 1:D
            (ai, p, Mf) = cols[k]
            slots = slotlists[ai]
            v = _mfslice(tens[ai][p], Mf)
            for j in eachindex(slots)
                slots[j].factor.l == 0 && continue          # D⁰ = 1 (trace axes)
                v = AngularMomentum.nmode_mul(v, _wig(wcache, slots[j].factor.l, g), j)
            end
            a2 = _assignment_image(assignments[ai], perm)
            ai2 = findfirst(==(a2), assignments)
            ai2 === nothing && error("assignment set not closed under the stabilizer")
            σ = _slot_sigma(slots, slotlists[ai2], perm)
            v = permutedims(v, invperm(σ))
            for p2 in eachindex(tens[ai2]), Mf2 = 1:(2Lf + 1)
                c = _frob(_mfslice(tens[ai2][p2], Mf2), v)
                abs(c) < _ROUNDOFF_RTOL && continue    # the tensors are orthonormal
                P[colidx[(ai2, p2, Mf2)], k] += c
            end
        end
    end
    P ./= length(stab)
    P .= (P .+ P') ./ 2
    F = eigen(Symmetric(P))
    all(λ -> abs(λ) < 1e-6 || abs(λ - 1) < 1e-6, F.values) ||
        error("decor SALC projector is not idempotent (eigenvalues $(F.values))")
    V1 = F.vectors[:, findall(>(0.5), F.values)]
    to_cartesian = zeros(Float64, D, D)
    for k = 1:D
        (ai, p, Mf) = cols[k]
        slots = slotlists[ai]
        v = _mfslice(tens[ai][p], Mf)
        for j in eachindex(slots)
            slots[j].factor.l == 0 && continue          # D⁰ = 1 (trace axes)
            v = AngularMomentum.nmode_mul(v, _wig(fcache, slots[j].factor.l, 1), j)
        end
        for p2 in eachindex(tens[ai]), Mf2 = 1:(2Lf + 1)
            c = _frob(_mfslice(tens[ai][p2], Mf2), v)
            abs(c) < _ROUNDOFF_RTOL && continue
            to_cartesian[colidx[(ai, p2, Mf2)], k] = c
        end
    end
    W = to_cartesian * _canonical_basis(to_cartesian' * V1)

    blocks = Vector{Tuple{Int,Array{Float64}}}[]
    for b = 1:size(W, 2)
        c = view(W, :, b)
        terms = Tuple{Int,Array{Float64}}[]
        for ai in eachindex(assignments)
            isempty(tens[ai]) && continue
            Fo = zeros(Float64, size(tens[ai][1])[1:(end - 1)]...)
            for p in eachindex(tens[ai]), Mf = 1:(2Lf + 1)
                coef = c[colidx[(ai, p, Mf)]]
                coef == 0.0 && continue
                Fo .+= coef .* _mfslice(tens[ai][p], Mf)
            end
            _zero_roundoff!(Fo)
            # `c` and the coupled tensors are orthonormal, so Σ_ai ‖Fo‖² = 1 and a
            # term is round-off when its norm is round-off against 1.
            norm(Fo) > _ROUNDOFF_RTOL && push!(terms, (ai, Fo))
        end
        push!(blocks, terms)
    end
    return blocks
end

# Decor-general `_transport_term`: rotate every slot axis, relabel sites through
# the connecting permutation, and return the member term in its canonical slot
# order.
function _transport_term_decors(a::Vector{SiteDecor}, slots::Vector{Slot},
                                F::Array{Float64}, g::Int, perm::Vector{Int},
                                wcache::_WigCache)::SALCTerm
    T = F
    for j in eachindex(slots)
        slots[j].factor.l == 0 && continue                  # D⁰ = 1 (trace axes)
        T = AngularMomentum.nmode_mul(T, _wig(wcache, slots[j].factor.l, g), j)
    end
    a2 = _assignment_image(a, perm)
    slots2 = _assignment_slots(a2)
    σ = _slot_sigma(slots, slots2, perm)
    G = _zero_roundoff!(Array{Float64}(permutedims(T, invperm(σ))))
    return SALCTerm(slots2, G)
end

# All SALCs of one cluster orbit for explicitly-given decoration labels (sorted
# `SiteDecor` multisets). Callers hand the labels in; nothing enumerates them
# yet. `isotropy = true` keeps only L_S = 0.
#
# `isotropy` is a REQUIRED KEYWORD on purpose. Upstream SLCE.jl's engine takes
# `soc::Bool` in this positional slot with the opposite meaning (`soc = true`
# keeps every L_S; `isotropy = true` keeps L_S = 0 only), so a verbatim copy of
# an upstream call would compile here and silently invert the screen. As a
# keyword with no default, a positional copy is a MethodError and an omitted
# screen an UndefKeywordError — both loud.
function _orbit_salcs_decors(crystal::Crystal, spacegroup::SpaceGroup, N::Int,
                             orbit_id::Int, O::ClusterOrbit,
                             labels::Vector{Vector{SiteDecor}}, wcache::_WigCache,
                             fcache::_WigCache; isotropy::Bool,
                             admit::Union{Nothing,Function} = nothing)::Vector{SALC}
    # `admit(t::Vector{SiteDecor})::Bool` screens a whole site-permutation orbit of
    # assignments, at the canonical representative and before any block index is
    # consumed. A label is a sorted multiset, so it cannot say WHICH site may carry
    # which rank; the pure-spin engine expresses that through its per-species `lmax`
    # (`_enumerate_ls`), and a caller that needs the same (or the pointed basis'
    # mark-aware variant) hands it in here. Correctness requirement: the verdict must
    # be a permutation-orbit invariant — stabilizer permutations preserve species and
    # edge geometry, so any rule built from (decor, species, edges-from-site) data is,
    # while one keyed on a bare site index is not.
    out = SALC[]
    rep = O.representative
    stab = _stabilizer(crystal, spacegroup, rep)
    perms = unique([perm for (_, perm) in stab])
    conns = _connect_all(crystal, spacegroup, rep, O.members)
    allunique(labels) || throw(ArgumentError(
        "duplicate decoration labels: each label projects to the same SALCs " *
        "twice (exactly collinear design columns)"))
    blockcount = Dict{Tuple{Vector{SiteDecor},Int,Int},Int}()
    for label in labels
        length(label) == N ||
            throw(ArgumentError("decor label has $(length(label)) sites for an " *
                                "$N-body orbit"))
        issorted(label) ||
            throw(ArgumentError("decor label must be sorted (canonical multiset)"))
        iseven(sum(d.spin_l for d in label)) ||
            throw(ArgumentError("Σl_spin must be even (time-reversal screen); " *
                                "got label $label"))
        # Distinct site assignments of the multiset, one canonical representative
        # per orbit of the site-permutation group — in EXACTLY `_enumerate_ls`'s
        # order (upstream review blocker): orbits are discovered in colex
        # (site-1-fastest, i.e. `Iterators.product`) order, the emitted
        # representative is the lex-min of the orbit, and the assignment list is
        # `unique(rep[p])`. Any other order relabels `block` indices / permutes
        # the gauge columns on pure-spin labels, silently breaking key-addressed
        # coefficient re-pairing the moment this engine backs a public builder.
        arrangements = sort(_multiset_arrangements(label);
                            by = a -> reverse!([_decortuple(d) for d in a]))
        seen = Set{Vector{SiteDecor}}()
        for tarr in arrangements
            tarr in seen && continue
            orbit = unique([tarr[p] for p in perms])
            for o in orbit
                push!(seen, o)
            end
            t = minimum(orbit)                   # lex-min canonical representative
            admit === nothing || admit(t)::Bool || continue
            assignments = unique([t[p] for p in perms])
            slotlists = [_assignment_slots(a) for a in assignments]
            cbs = [_decor_coupled_bases(sl, isotropy) for sl in slotlists]
            # `isotropy` already screened the paths inside `_decor_coupled_bases`,
            # before any tensor was built, so `blockset` carries only live blocks.
            blockset = sort(unique((ls, lf) for cbo in cbs for (ls, lf, _) in cbo))
            for (L_S, Lf) in blockset
                blocks = _project_and_fold_decors(stab, assignments, slotlists,
                                                  cbs, L_S, Lf, wcache, fcache)
                for terms_rep in blocks
                    isempty(terms_rep) && continue
                    members = SALCMember[]
                    for (m, (g, perm)) in zip(O.members, conns)
                        mterms = SALCTerm[]
                        for (ai, F) in terms_rep
                            push!(mterms,
                                  _transport_term_decors(assignments[ai],
                                                         slotlists[ai], F, g,
                                                         perm, wcache))
                        end
                        push!(members, SALCMember(m.atoms, m.shifts, mterms))
                    end
                    members = _canonicalize_members(members)
                    isempty(members) && continue
                    block = get(blockcount, (label, L_S, Lf), 0) + 1
                    blockcount[(label, L_S, Lf)] = block
                    key = SALCKey(N, orbit_id, copy(label), L_S, Lf, block)
                    push!(out, SALC(key, N, key.decors, L_S, Lf, members))
                end
            end
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Function-space reduction of one orbit's SALCs.
#
# `evaluate_salc` reads a member's `atoms` but never its `shifts`: under plain
# periodic evaluation every lattice image of an atom carries the same spin, so two
# members differing only in `shifts` contribute to the SAME monomials
# `∏ᵢ Z_{lsᵢμᵢ}(e_{atomᵢ})`. When an orbit folds distinct cluster instances onto one
# atom set — WS-boundary ties kept whole, or a widened `tie_tol` merging a near-tie
# shell — their tensors aggregate, and the aggregate can vanish or become linearly
# dependent across the orbit's SALCs even though every per-instance tensor is
# nonzero and the SALCs are individually space-group invariant. The design matrix
# then silently loses rank: OLS returns huge mutually-cancelling coefficients with
# normal-looking energies and R², and every bond-resolved readout (`bilinear_terms`,
# Sunny export, MC Hamiltonians) stops being unique. Measured on bulk MnTe with SOC
# (3×3×3, 108 atoms, P6_3/mmc): 51 SALCs emitted, 14 of them aggregating to zero
# (relative residual ≤ 5.3e-12), max|coef| ~1e7 from OLS — invisible from the fit
# diagnostics (found 2026-08-12; the "what a tie costs" phenomenon of
# docs/src/theory/resolvability.md reaching the basis builder).
#
# The reduction below is exact linear algebra, not a numerical heuristic: each SALC
# is expanded into its aggregated coefficient vector keyed by
# `(atoms, ls-assignment, tensor index)`. For members with all-distinct atoms these
# monomials are an orthogonal function family (products of orthonormal `Z_{lμ}` on
# distinct unit spheres; the atom set is recoverable from the key because cluster
# sites always carry `l ≥ 1` — `_enumerate_ls` starts every per-site range at 1, so
# no `l = 0` factor can make a pair monomial collide with a single-site one), so
# dependence of the vectors is EQUIVALENT to dependence of the Φ's. With a repeated
# atom in a member (`AllImages` self-pairs) the monomials themselves may be
# dependent, so vector-dependence still implies function-dependence (dropping stays
# sound) but independence is not certified — the OLS rank warning is the backstop
# there. Dependence ACROSS orbits is what this per-orbit pass never sees; the
# proportional case (tied images in distinct orbits) is picked up on the finished
# basis by `column_ties` (`basis/aliases.jl`), which reuses `_function_vector`.
#
# The rank decision is taken per BLOCK, not per column (2026-11-16; ported in spirit from
# SLCE.jl a248d4b, whose build makes the vanishing subspace of a block columns of its
# own and freezes them — this package drops them instead). A block is the set of
# columns that share every label of their key but `block`: one orbit, one `l`-multiset,
# one `Lf`. The symmetry fixes the space a block spans, not the columns inside it, so a
# test that walks the columns one by one, each against its own norm, answers a
# question about the gauge: in a crystal turned by a small angle θ a column whose
# aggregate is θ-small was kept and its coefficient came out about 1/θ too large
# (per-bond `bilinear_terms` off by 2.3e6 at θ = 1e-6 on a tetragonal pair basis),
# and at θ = 1e-8 / 1e-7 the count itself was wrong (6 instead of 4). Now each block's
# aggregated functions, projected off the span already kept in this orbit (blocks in
# emission order, so a lower `Lf` survives as before), are decided by one SVD against
# the block's GROSS mass — the largest unaggregated tensor norm of its columns — at
# `_CANCELLATION_RTOL`, the reading of the pointed gate (`moment_resolvability`) and of
# upstream's per-column freeze. What survives is `canon(N1)`, the complement of the
# vanishing and dependent directions gauged in the coordinates of the block's own
# columns; those columns are already the crystal-frame gauge, so this choice does not
# depend on how the crystal is turned either. A singular value between the cut and
# `_REGAUGE_BAND_TOP` is refused, never guessed (`_block_split`, `basis/regauge.jl`).
#
# Measured margins (2026-11-16, the tied fixtures of `test_gauge.jl` along the axes and
# turned): cancellations ≤ 1.6e-15 of the gross mass, content ≥ 1.0 — the gate in that
# file holds both with three decades of headroom. The zero residual is bounded by the
# entry prune in `_project_and_fold` / `_transport_term` (`_ROUNDOFF_RTOL`, relative
# at 1e-13 per entry per member), so it scales like √(K·M)·1e-13 with the key count K
# and member count M — far below the cut at any orbit size this package builds.
#
# Dict iteration order (insertion order, deterministic in Julia) fixes the
# floating-point summation order in `_dictnorm`/`_dictdot`/the MGS update; only
# threshold comparisons and the rank decision consume it.
#
# `_AGG_DEP_RTOL` is the dependence cut of the cross-orbit span-collapse REPORT
# (`_span_collapsed`, `basis/aliases.jl`), judged against each vector's own norm; it
# no longer decides anything inside an orbit.
const _AGG_DEP_RTOL = 1e-8

# The aggregated (shift-blind) monomial coefficient vector of a SALC, keyed by
# `(member atoms, spin-axis ranks, tensor index)`.
const _FunctionVector = Dict{Tuple{Vector{Int},Vector{Int},Int},Float64}

# Aggregated (shift-blind) monomial coefficient vector of one SALC, plus the
# unaggregated tensor norm used as the zero-test scale.
function _function_vector(s::SALC)
    # Pure spin only. The aggregate key below reads the SPIN-axis ranks and the
    # member atoms, which names a function uniquely ONLY while every slot is a spin
    # slot; on a decorated term two different terms can share a key and the
    # reduction would merge or drop them silently. Refusing here is the gate that
    # keeps the decor engine out of `_reduce_orbit_salcs` until it is wired in on
    # purpose — and the message is the checklist for doing so.
    all(is_pure_spin, s.decors) || throw(ArgumentError(
        "function-space reduction of a displacement-decorated SALC is not wired: " *
        "the aggregate key reads SPIN ranks only. Wiring checklist — (a) key on " *
        "the full slot list (channel, site, k, l) with a content-based " *
        "hash(::Slot) [M1(a)]; (b) the mixed-label 3-cycle invariance gate on the " *
        "C3v triangle [C-1(2)]; (c) the L_S block-diagonality gate over all " *
        "(L_S, Lf) blocks at once [C-2]; (d) an lsum cap for decorated labels and " *
        "its anti-drift run [MIN-7]"))
    v = _FunctionVector()
    raw2 = 0.0
    for m in s.members, t in m.terms
        raw2 += sum(abs2, t.folded)
        # The SPIN-axis ranks stand in for the old per-site `ls`: this engine only
        # ever builds pure-spin terms, whose slots are the identity map, so the
        # partition is the one it always was. Hoisted out of the entry loop —
        # `_term_spin_ls` allocates where the v4 per-term `ls` was a field read.
        # Keying on `t.slots` directly would save that allocation and cost 44%
        # wall time (measured): `hash(::Slot)` goes through `objectid`, and this
        # dict is looked up once per nonzero tensor entry.
        tls = _term_spin_ls(t)
        flat = vec(t.folded)                    # reshape view, no copy
        @inbounds for li in eachindex(flat)
            c = flat[li]
            c == 0.0 && continue
            k = (m.atoms, tls, li)
            v[k] = get(v, k, 0.0) + c
        end
    end
    return v, sqrt(raw2)
end

_dictnorm(v::Dict) = sqrt(sum(abs2, values(v); init = 0.0))
function _dictdot(a::Dict, b::Dict)
    # iterate the smaller dict
    length(a) > length(b) && return _dictdot(b, a)
    acc = 0.0
    for (k, x) in a
        acc += x * get(b, k, 0.0)
    end
    return acc
end

# Modified Gram–Schmidt step: project `v` (in place) against the orthonormalized
# `kept` set; returns `(‖v‖ after, ‖v‖ before)`. Shared by the per-orbit reduction and
# the cross-orbit span-collapse test in `basis/aliases.jl` — one loop, one summation
# order (see the file header on why that order is fixed).
function _mgs_residual!(v::_FunctionVector, kept::Vector{_FunctionVector})
    nv0 = _dictnorm(v)
    for q in kept
        c = _dictdot(q, v)
        c == 0.0 && continue
        for (k, x) in q
            v[k] = get(v, k, 0.0) - c * x
        end
    end
    return _dictnorm(v), nv0
end

# What the reduction did to one block that lost directions: the block's keys as
# emitted, how many directions vanish on their own (`zero`, the nullity of the block's
# aggregate) and how many depend on blocks kept earlier in the orbit (`dependent`), and
# — when the survivors are an exact selection of the block's columns — the dropped
# keys with their reasons (`nothing` when the survivors recombine the block, so no
# single key names a dropped direction).
const _BlockDrop = @NamedTuple{keys::Vector{SALCKey},zero::Int,dependent::Int,
                               dropped::Union{Nothing,Vector{Tuple{SALCKey,Symbol}}}}

# Aggregated function vectors as dense columns over the union of their keys, with the
# row index of each key.
function _dense_block(vecs::Vector{_FunctionVector})
    rows = Dict{Tuple{Vector{Int},Vector{Int},Int},Int}()
    for v in vecs, k in keys(v)
        get!(rows, k, length(rows) + 1)
    end
    S = zeros(Float64, length(rows), length(vecs))
    for (j, v) in enumerate(vecs), (k, x) in v
        S[rows[k], j] = x
    end
    return S, rows
end

"""
    _reduce_orbit_salcs(salcs) -> (kept, drops, regauged, ratios)

The function-space reduction of one orbit, block by block (see the comment above
`_AGG_DEP_RTOL`). Blocks are taken in the emission order of their first column; each
block's aggregated functions are projected off the span already kept and its rank is
decided against the block's gross mass (`_block_split`). A block that loses nothing is
returned untouched; one that loses everything is dropped whole; otherwise its
survivors are `canon(N1)` in the coordinates of its columns (`_canonical_basis` with
the round-off of a unit vector set to zero), each under the key of the column it was
pivoted on — a survivor that is a single column is that column, bit for bit, under its
own key. Returns the kept SALCs, one `_BlockDrop`
per block that lost directions, and the keys of the kept columns of every block a
proper, nonzero part of which was dropped (the blocks whose gauge this cell chose),
and every relative singular value the rank decisions read (`ratios`, for the margin
gate of the test suite).
"""
function _reduce_orbit_salcs(salcs::Vector{SALC})
    isempty(salcs) && return salcs, _BlockDrop[], SALCKey[], Float64[]
    # The pure-spin door comes first, so a decorated orbit is refused whether or not it
    # has a tie (`_function_vector` throws with the wiring checklist).
    for s in salcs
        all(is_pure_spin, s.decors) || _function_vector(s)
    end
    # Without a tie nothing can go: every member of a column reaches its own atom
    # multiset, so the aggregates are the member tensors themselves, orthonormal
    # across the orbit's columns (orthogonal transport of an orthonormal set). The
    # orbit is returned as built, bit for bit, without an expansion.
    _block_has_tie(salcs) || return salcs, _BlockDrop[], SALCKey[], Float64[]
    groups = Dict{Any,Vector{Int}}()
    order = Any[]
    for (j, s) in enumerate(salcs)
        g = _block_key(s)
        haskey(groups, g) || push!(order, g)
        push!(get!(groups, g, Int[]), j)
    end
    fv = [_function_vector(s) for s in salcs]
    kept = SALC[]
    keptvecs = _FunctionVector[]             # orthonormal basis of the span kept so far
    drops = _BlockDrop[]
    regauged = SALCKey[]
    ratios = Float64[]
    for g in order
        index = sort(groups[g]; by = j -> salcs[j].key.block)
        block = salcs[index]
        k = length(index)
        gross = [fv[j][2] for j in index]
        own = [fv[j][1] for j in index]
        projected = [copy(v) for v in own]
        for v in projected
            _mgs_residual!(v, keptvecs)
        end
        # one row space for both, so the counts below compare like with like
        S, rows = _dense_block(vcat(own, projected))
        S0 = S[:, 1:k]
        Sp = S[:, (k + 1):(2k)]
        split = _block_split(Sp, gross)
        _refuse_undecidable(split.ratios, block[1].key)
        append!(ratios, split.ratios)
        r = size(split.kept, 2)
        if r == k                            # nothing lost: the block as projected
            append!(kept, block)
            _extend_kept!(keptvecs, Sp, rows, split)
            continue
        end
        n_zero = min(k - r, count(<=(_CANCELLATION_RTOL), _block_split(S0, gross).ratios))
        n_dep = (k - r) - n_zero
        if r == 0                            # the whole block goes
            push!(drops, (; keys = [s.key for s in block], zero = n_zero,
                            dependent = n_dep,
                            dropped = [(s.key, _drop_reason(fv[j][1], gross))
                                       for (s, j) in zip(block, index)]))
            continue
        end
        W, pivots = _canonical_basis_pivots(split.kept)
        for i in eachindex(W)
            abs(W[i]) < _ROUNDOFF_RTOL && (W[i] = 0.0)
        end
        # A survivor that is a unit vector up to round-off is that column, copied; any
        # other is the combination, keyed by its pivot column.
        unit = zeros(Int, r)
        for b = 1:r
            nz = findall(!iszero, view(W, :, b))
            length(nz) == 1 && (unit[b] = nz[1])
        end
        for b = 1:r
            push!(kept, unit[b] > 0 ? block[unit[b]] :
                        _combine_salcs(block, view(W, :, b), block[pivots[b]].key))
            push!(regauged, kept[end].key)
        end
        _extend_kept!(keptvecs, Sp, rows, split)
        dropped = all(>(0), unit) ?
                  [(block[j].key, _drop_reason(fv[index[j]][1], gross))
                   for j = 1:k if !(j in unit)] : nothing
        push!(drops, (; keys = [s.key for s in block], zero = n_zero, dependent = n_dep,
                        dropped))
    end
    return kept, drops, regauged, ratios
end

# Why a single dropped column went: its own aggregate cancels against the block's
# gross mass, or it is a function the orbit already has.
_drop_reason(v::_FunctionVector, gross::AbstractVector{Float64})::Symbol =
    _dictnorm(v) <= _CANCELLATION_RTOL * maximum(gross) ? :zero : :dependent

# Add an orthonormal basis of the span a block contributes (its projected aggregate
# restricted to the kept right singular vectors) to the orbit's kept span. Thin
# Householder QR of `Sp·V_kept` keeps the stored vectors orthonormal to round-off, so
# the projections of later blocks stay at round-off of their own mass.
function _extend_kept!(keptvecs::Vector{_FunctionVector}, Sp::Matrix{Float64},
                       rows::Dict, split)
    size(split.kept, 2) == 0 && return keptvecs
    Q = Matrix(qr(Sp * split.kept).Q)
    keyof = Vector{Tuple{Vector{Int},Vector{Int},Int}}(undef, length(rows))
    for (key, r) in rows
        keyof[r] = key
    end
    for c in axes(Q, 2)
        v = _FunctionVector()
        for r in axes(Q, 1)
            Q[r, c] == 0.0 || (v[keyof[r]] = Q[r, c])
        end
        push!(keptvecs, v)
    end
    return keptvecs
end

# All SALCs of one cluster orbit. Self-contained (its `blockcount` and output are
# orbit-local; `wcache` is read-only), so orbits are processed independently — the
# unit of parallelism in `build_salc_basis`.
function _orbit_salcs(crystal::Crystal, spacegroup::SpaceGroup, N::Int, orbit_id::Int,
                      O::ClusterOrbit, lmax::Vector{Int}, lsumN::Int, isotropy::Bool,
                      wcache::_WigCache, fcache::_WigCache)::Vector{SALC}
    out = SALC[]
    rep = O.representative
    stab = _stabilizer(crystal, spacegroup, rep)
    perms = unique([perm for (_, perm) in stab])
    conns = _connect_all(crystal, spacegroup, rep, O.members)
    # `block` runs across ALL canonical l-tuples sharing one sorted label so the key
    # stays injective: when the site permutations are a proper subgroup of Sₙ they
    # split a degenerate multiset's arrangements into several orbits (e.g. (1,1,2) and
    # (2,1,1) on a mirror-only triangle put l=2 on inequivalent sites), all sharing
    # `ls = sort(t)` but giving distinct SALCs.
    blockcount = Dict{Tuple{Vector{Int},Int},Int}()
    for t in _enumerate_ls(N, O.species, lmax, lsumN, perms)
        orderings = unique([t[p] for p in perms])
        # Build the coupled bases for each ordering once. They are reused across every
        # `Lf`; rebuilding them inside `_project_and_fold` per `Lf` recomputed the whole
        # chained-CG construction `|Lfset|` times. The `Lf` set is the CG decomposition
        # of the `l`-multiset (permutation-invariant), so the union over orderings equals
        # the previous `coupled_bases(t)` set.
        cbs = [coupled_bases(o; isotropy = isotropy) for o in orderings]
        Lfset = sort(unique(cb.Lf for cbo in cbs for cb in cbo))
        lab = sort(collect(t))
        for Lf in Lfset
            blocks = _project_and_fold(stab, orderings, cbs, Lf, wcache, fcache)
            for terms_rep in blocks
                isempty(terms_rep) && continue
                members = SALCMember[]
                for (m, (g, perm)) in zip(O.members, conns)
                    mterms = SALCTerm[]
                    for (o, F) in terms_rep
                        mls, G = _transport_term(o, F, g, perm, wcache)
                        push!(mterms, SALCTerm(spin_slots(mls), G))
                    end
                    push!(members, SALCMember(m.atoms, m.shifts, mterms))
                end
                # Fold the ordered, anchored images into the canonical duplicate-free
                # form (exact regrouping; empty ⇒ the channel vanished — skip it, as
                # the previous per-term norm check did).
                members = _canonicalize_members(members)
                isempty(members) && continue
                block = get(blockcount, (lab, Lf), 0) + 1
                blockcount[(lab, Lf)] = block
                # Pure-spin construction: decors = spin_decors(lab), L_S = Lf
                # (the p = 0 edge of the joint key layout).
                key = SALCKey(N, orbit_id, spin_decors(lab), Lf, Lf, block)
                push!(out, SALC(key, N, key.decors, Lf, Lf, members))
            end
        end
    end
    return out
end

"""
    build_salc_basis(crystal, spacegroup, clusters; lmax_by_species,
                     lsum_by_body = nothing, isotropy = false) -> SALCBasis

Construct the symmetry-adapted (and time-reversal-even) SCE basis for every cluster
orbit and body order. For each `(orbit, l-multiset, Lf)` the stabilizer-invariant
coefficient subspace (over orderings × coupling paths) is found, gauge-fixed, and
transported to all orbit members.

The columns inside a block (the SALCs sharing `(orbit, l-multiset, Lf)`) are chosen
along axes fixed to the crystal — the vectors of the symmetrized lattice made
orthonormal in order (`symmetrized_lattice(basis).crystal_axes`) — so the columns of a
crystal turned as a whole are the turned columns. They still follow the lattice
vectors: another choice of lattice vectors for the same crystal can give other
columns inside a block (same span, same model).

Each orbit is then **function-space reduced, block by block**: the combinations of a
block's columns whose orbit sum aggregates to the zero function of the cell-periodic
spins, or depends on blocks kept earlier in the orbit (a supercell folding distinct
cluster instances — WS-boundary ties, merged near-tie shells — onto the same atom
sets), are dropped with a warning naming them: by key when the survivors are single
columns, as a block with a count otherwise. The block's rank is decided against its
gross mass, and a value that is neither a cancellation nor content is refused. The
survivors are an orthonormal complement chosen in the coordinates of the block's own
columns; a block that loses nothing keeps its columns bit for bit, and one whose
survivors are single columns keeps exactly those, under their own keys, so gaps in
`block` are legal. The columns of a block that lost a proper part are listed by
[`regauged_columns`](@ref): they belong to this cell, not to the crystal. The
reduction is exact (the spanned model space is unchanged) and, **per orbit and for
members with all-distinct atoms**, guarantees the orbit's surviving functions are
linearly independent. Dependence **across** orbits — tied
images the space group does not relate, placed in different orbits — is not folded
here: `SCEBasis` detects such cross-orbit alias groups on the built basis
(`column_ties`, `basis/aliases.jl`) and ties their proportional columns into one
design column, read back with an equal per-bond split by convention. What neither
pass covers — and where the design matrix can still lose rank, caught only by the
`OLS` rank warning: cross-orbit classes that are dependent without being
proportional (reported as `:span_collapsed`), repeated-atom members (`AllImages`
self-pairs; dropping stays sound but independence is uncertified), and training
sets with fewer rows than columns.

# Keyword arguments
- `lmax_by_species::AbstractVector{<:Integer}`: per-species maximum `l`.
- `lsum_by_body`: per-body-order cap on `Σl` over the cluster sites, indexed by
  body order (`nothing` = no cap; entries of `LSUM_UNCAPPED` mean no cap for
  that order).
- `isotropy::Bool = false`: keep only the scalar `Lf == 0` channel if `true`.

# Status
Arbitrary body order. Isotropic (`Lf == 0`) and anisotropic (`Lf > 0`) channels —
including those that mix coupling paths (`N ≥ 3`) or `l`-orderings on
symmetry-equivalent sites (e.g. `l=(1,1,2)`, `Lf>0` on an equilateral triangle) —
are validated by the ground-truth invariance test `Φ(g·e)=Φ(e)` (non-collinear
spins, all `Lf`) and time-reversal evenness `Φ(−e)=Φ(e)`.

The orbits are built in parallel over `Threads.nthreads()` (set `julia -t` /
`JULIA_NUM_THREADS`); each orbit is independent and the result is sorted by
`SALCKey`, so it is byte-for-byte identical at any thread count.
"""
function build_salc_basis(crystal::Crystal, spacegroup::SpaceGroup, clusters::ClusterSet;
                          lmax_by_species::AbstractVector{<:Integer},
                          lsum_by_body::Union{Nothing,AbstractVector{<:Integer}} = nothing,
                          isotropy::Bool = false)::SALCBasis
    lmax = collect(Int, lmax_by_species)
    maxbody = isempty(clusters.by_body) ? 0 : maximum(keys(clusters.by_body))
    lsum_by_body === nothing || length(lsum_by_body) >= maxbody ||
        throw(ArgumentError("lsum_by_body has $(length(lsum_by_body)) entries; " *
                            "the clusters reach body order $maxbody"))
    lsumN(N) = lsum_by_body === nothing ? LSUM_UNCAPPED : Int(lsum_by_body[N])
    # Precompute the Wigner-D caches (operations and crystal frame) serially so they are
    # read-only (race-free) below.
    maxl = isempty(lmax) ? 0 : maximum(lmax)
    wcache = _build_wig_cache(spacegroup, maxl)
    fcache = _build_frame_cache(crystal, spacegroup, maxl)
    # Flatten the orbits into one work list, keeping `(N, orbit_id)` for the SALC keys.
    # Orbits are independent; the final `sort!` by key makes the result identical at any
    # thread count, so the parallel order is irrelevant.
    work = Tuple{Int,Int,ClusterOrbit}[]
    for N in sort(collect(keys(clusters.by_body)))
        for (orbit_id, O) in enumerate(clusters.by_body[N])
            push!(work, (N, orbit_id, O))
        end
    end
    parts = Vector{Vector{SALC}}(undef, length(work))
    drops = Vector{Vector{_BlockDrop}}(undef, length(work))
    touched = Vector{Vector{SALCKey}}(undef, length(work))
    Threads.@threads for w in eachindex(work)
        (N, orbit_id, O) = work[w]
        raw = _orbit_salcs(crystal, spacegroup, N, orbit_id, O, lmax, lsumN(N),
                           isotropy, wcache, fcache)
        parts[w], drops[w], touched[w], _ = _reduce_orbit_salcs(raw)
    end
    salcs = isempty(parts) ? SALC[] : reduce(vcat, parts)
    sort!(salcs; by = s -> s.key)
    dropped = isempty(drops) ? _BlockDrop[] : reduce(vcat, drops)
    isempty(dropped) || _warn_reduction(sort!(dropped; by = d -> first(d.keys)))
    keyvec = SALCKey[s.key for s in salcs]
    regauged_keys = Set(isempty(touched) ? SALCKey[] : reduce(vcat, touched))
    regauged = [j for j in eachindex(keyvec) if keyvec[j] in regauged_keys]
    return SALCBasis(salcs, keyvec; regauged)
end

# The reduction's report: every dropped direction by name — the key of a dropped
# column when the survivors of its block are an exact selection, otherwise the block
# and a count (the survivors then recombine the block, and no single key names what
# went).
function _warn_reduction(dropped::Vector{_BlockDrop})
    nz = sum(d.zero for d in dropped)
    nd = sum(d.dependent for d in dropped)
    lines = String[]
    for d in dropped
        k = first(d.keys)
        head = "  $(k.body)-body orbit $(k.orbit_id), ls = $(spin_ls(k)), Lf = $(k.Lf)"
        if d.dropped === nothing
            blocks = join((key.block for key in d.keys), ", ")
            push!(lines, head * ", blocks $blocks: $(d.zero + d.dependent) " *
                         "combination(s) dropped ($(d.zero) zero, $(d.dependent) " *
                         "dependent); the $(length(d.keys) - d.zero - d.dependent) " *
                         "kept column(s) recombine the block")
        else
            for (key, why) in d.dropped
                push!(lines, head * ", block $(key.block): $why")
            end
        end
    end
    @warn "SALC reduction: dropped $(nz + nd) redundant basis " *
          "function(s) ($nz aggregating to identically zero, $nd linearly " *
          "dependent as functions of the spins). The supercell folds distinct " *
          "cluster instances of the affected orbit(s) onto the same atom sets " *
          "(WS-boundary ties, or a merged near-tie shell), so their tensors " *
          "aggregate — the dropped combinations are not resolvable from this " *
          "supercell, and keeping them would make the design matrix rank " *
          "deficient (non-unique coefficients). The reduction is exact: the " *
          "spanned model space and all energy/torque predictions are " *
          "unchanged. For a `dependent` drop, note the surviving channel's " *
          "fitted coefficient absorbs the dropped one's — the bond-level " *
          "attribution within the affected orbit remains unresolvable from " *
          "this supercell, with or without the drop. A block that lost part of " *
          "itself is listed by `regauged_columns` (and flagged in `coeftable`): its " *
          "columns belong to this cell, not to the crystal.\n" * join(lines, "\n")
end
