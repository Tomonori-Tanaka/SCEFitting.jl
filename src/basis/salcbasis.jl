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
# The absolute drop thresholds below (1e-10/1e-12) and the idempotency guard (1e-6)
# are sized for the validated regime (per-site l ≤ 2, small body order), where the
# CG/recoupling products stay O(1e-3)–O(1). The ground-truth invariance test and the
# exact-0/1 eigenvalue assertion are the in-band gates; at much higher l these
# thresholds would want to scale with the coupled-tensor magnitude.

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

# Project the combined (ordering, path, Mf) space onto stabilizer invariants for a
# fixed final `Lf`, gauge-fix, and fold each invariant into per-ordering tensors.
# Returns a list of blocks; each block is `Vector{(ls, folded)}` (the rep terms).
function _project_and_fold(stab::Vector{Tuple{Int,Vector{Int}}},
                           orderings::Vector{Vector{Int}}, cbs::Vector, Lf::Int,
                           wcache::_WigCache)
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
                abs(c) < 1e-12 && continue
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
    W = _canonical_basis(V1)

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
            @inbounds for idx in eachindex(Fo)
                abs(Fo[idx]) < 1e-10 && (Fo[idx] = 0.0)
            end
            norm(Fo) > 1e-10 && push!(terms, (collect(o), Fo))
        end
        push!(blocks, terms)
    end
    return blocks
end

# Deterministic, BLAS-independent gauge: axis-pivoted modified Gram–Schmidt on the
# block projector, with a first-significant-component sign fix.
function _canonical_basis(V1::AbstractMatrix{Float64})
    n, d = size(V1)
    d == 0 && return zeros(Float64, n, 0)
    Q = V1 * V1'
    W = zeros(Float64, n, d)
    k = 0
    for j = 1:n
        u = Q[:, j]
        for i = 1:k
            u .-= dot(view(W, :, i), u) .* view(W, :, i)
        end
        nu = norm(u)
        if nu > 1e-8
            k += 1
            W[:, k] = u ./ nu
            _sign_canon!(view(W, :, k))
        end
        k == d && break
    end
    k == d || error("canonical gauge failed to span the invariant subspace")
    return W
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
    G = Array{Float64}(permutedims(T, q))
    @inbounds for idx in eachindex(G)
        abs(G[idx]) < 1e-10 && (G[idx] = 0.0)
    end
    return o[q], G                          # member ls = o[invperm(perm)]
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
# there, as it is for dependence ACROSS orbits, which this per-orbit pass never
# sees.
#
# Thresholds: an aggregated function is dropped as zero when its vector norm is
# ≤ `_AGG_ZERO_RTOL ×` its unaggregated tensor norm, and as dependent when its
# MGS residual against the kept set is ≤ `_AGG_DEP_RTOL ×` its own aggregated
# norm. Measured gap on the MnTe case: true zeros cancel to ≤ 5.3e-12 relative
# while surviving SALCs sit at O(1e-2)–O(1) — so 1e-8 has ~3 orders of headroom
# above the measured zeros and ~6 below the smallest genuine survivor. The zero
# residual is bounded by the 1e-10 absolute entry prune in `_project_and_fold` /
# `_transport_term` *per entry per member*, so it scales like √(K·M)·1e-10 with
# the key count K and member count M — like the prune constants themselves (see
# the file header), sized for the validated regime (per-site l ≤ 2, small body
# order) and to be revisited at much higher l or very large orbits. A miss is in
# the safe direction (a true zero kept, then caught by the OLS rank warning).
# Dict iteration order (insertion order, deterministic in Julia) fixes the
# floating-point summation order in `_dictnorm`/`_dictdot`/the MGS update; it is
# the one place in this file where that order follows hash-map layout rather than
# an explicit sort, which is fine because only threshold comparisons consume it.
const _AGG_ZERO_RTOL = 1e-8
const _AGG_DEP_RTOL = 1e-8

# Aggregated (shift-blind) monomial coefficient vector of one SALC, plus the
# unaggregated tensor norm used as the zero-test scale.
function _function_vector(s::SALC)
    v = Dict{Tuple{Vector{Int},Vector{Int},Int},Float64}()
    raw2 = 0.0
    for m in s.members, t in m.terms
        raw2 += sum(abs2, t.folded)
        flat = vec(t.folded)                    # reshape view, no copy
        @inbounds for li in eachindex(flat)
            c = flat[li]
            c == 0.0 && continue
            k = (m.atoms, t.ls, li)
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

# Keep the maximal leading (in emission order) subset of `salcs` whose aggregated
# functions are independent; return `(kept, dropped)` where each dropped entry is
# `(key, :zero | :dependent)`. Emission order puts lower `Lf` first within an
# `l`-tuple, so the scalar (Heisenberg-like) channel survives and the exotic
# aggregate-degenerate ones are the ones named as dropped. Orbit-local and
# deterministic — safe under the threaded orbit loop.
function _reduce_orbit_salcs(salcs::Vector{SALC})
    isempty(salcs) && return salcs, Tuple{SALCKey,Symbol}[]
    kept = SALC[]
    keptvecs = Dict{Tuple{Vector{Int},Vector{Int},Int},Float64}[]
    dropped = Tuple{SALCKey,Symbol}[]
    for s in salcs
        v, raw = _function_vector(s)
        nv0 = _dictnorm(v)
        if nv0 <= _AGG_ZERO_RTOL * raw
            push!(dropped, (s.key, :zero))
            continue
        end
        # modified Gram–Schmidt against the kept (orthonormalized) set
        for q in keptvecs
            c = _dictdot(q, v)
            c == 0.0 && continue
            for (k, x) in q
                v[k] = get(v, k, 0.0) - c * x
            end
        end
        nv = _dictnorm(v)
        if nv <= _AGG_DEP_RTOL * nv0
            push!(dropped, (s.key, :dependent))
            continue
        end
        for k in keys(v)
            v[k] /= nv
        end
        push!(keptvecs, v)
        push!(kept, s)
    end
    return kept, dropped
end

# All SALCs of one cluster orbit. Self-contained (its `blockcount` and output are
# orbit-local; `wcache` is read-only), so orbits are processed independently — the
# unit of parallelism in `build_salc_basis`.
function _orbit_salcs(crystal::Crystal, spacegroup::SpaceGroup, N::Int, orbit_id::Int,
                      O::ClusterOrbit, lmax::Vector{Int}, lsumN::Int, isotropy::Bool,
                      wcache::_WigCache)::Vector{SALC}
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
            blocks = _project_and_fold(stab, orderings, cbs, Lf, wcache)
            for terms_rep in blocks
                isempty(terms_rep) && continue
                members = SALCMember[]
                for (m, (g, perm)) in zip(O.members, conns)
                    mterms = SALCTerm[]
                    for (o, F) in terms_rep
                        mls, G = _transport_term(o, F, g, perm, wcache)
                        push!(mterms, SALCTerm(mls, G))
                    end
                    push!(members, SALCMember(m.atoms, m.shifts, mterms))
                end
                # Fold the ordered, anchored images into the canonical duplicate-free
                # form (exact regrouping; empty ⇒ the channel vanished — skip it, as
                # the previous per-term `norm > 1e-10` check did).
                members = _canonicalize_members(members)
                isempty(members) && continue
                block = get(blockcount, (lab, Lf), 0) + 1
                blockcount[(lab, Lf)] = block
                key = SALCKey(N, orbit_id, lab, Lf, block)
                push!(out, SALC(key, N, lab, Lf, members))
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

Each orbit is then **function-space reduced**: SALCs whose orbit sum aggregates to
the zero function of the cell-periodic spins, or becomes linearly dependent within
the orbit (a supercell folding distinct cluster instances — WS-boundary ties, merged
near-tie shells — onto the same atom sets), are dropped with a warning naming the
orbit, channel, and reason. The reduction is exact (the spanned model space is
unchanged) and, **per orbit and for members with all-distinct atoms**, guarantees
the orbit's surviving functions are linearly independent; surviving keys keep their
`block` numbers, so gaps in `block` are legal. What it does NOT cover — and where
the design matrix can still lose rank, caught only by the `OLS` rank warning:
dependence **across** orbits (a trivial/underreported space group can place tied
images in different orbits; more generally cross-orbit supercell aliasing, which is
unresolvable from the cell, not mergeable), repeated-atom members (`AllImages`
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
    # Precompute the Wigner-D cache serially so it is read-only (race-free) below.
    wcache = _build_wig_cache(spacegroup, isempty(lmax) ? 0 : maximum(lmax))
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
    drops = Vector{Vector{Tuple{SALCKey,Symbol}}}(undef, length(work))
    Threads.@threads for w in eachindex(work)
        (N, orbit_id, O) = work[w]
        raw = _orbit_salcs(crystal, spacegroup, N, orbit_id, O, lmax, lsumN(N),
                           isotropy, wcache)
        parts[w], drops[w] = _reduce_orbit_salcs(raw)
    end
    salcs = isempty(parts) ? SALC[] : reduce(vcat, parts)
    sort!(salcs; by = s -> s.key)
    dropped = isempty(drops) ? Tuple{SALCKey,Symbol}[] : reduce(vcat, drops)
    if !isempty(dropped)
        sort!(dropped; by = first)
        nz = count(d -> d[2] === :zero, dropped)
        nd = length(dropped) - nz
        detail = join(("  $(k.body)-body orbit $(k.orbit_id), ls = $(k.ls), " *
                       "Lf = $(k.Lf), block $(k.block): $why" for (k, why) in dropped),
                      "\n")
        @warn "SALC reduction: dropped $(length(dropped)) redundant basis " *
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
              "this supercell, with or without the drop.\n" * detail
    end
    return SALCBasis(salcs, SALCKey[s.key for s in salcs])
end
