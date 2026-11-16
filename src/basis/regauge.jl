# The gauge inside a block, decided by the reference cell where it has to be.
#
# A BLOCK is the set of columns that share every label of their key but `block`: one
# orbit, one decoration multiset, one `(L_S, Lf)` (and, in the pointed basis, one mark
# class). The symmetry fixes the space the block spans, not the columns inside it; the
# columns are a gauge, chosen along the crystal axes (`_project_and_fold`).
#
# On a finite reference cell a Wigner–Seitz boundary tie can make a SUBSPACE of a block
# vanish identically on every cell-periodic configuration. Whether that subspace happens
# to be spanned by single columns depends on the gauge, so every decision about it is
# taken over the block as a whole, against the block's gross mass:
#
# - the energy basis DROPS it (`_reduce_orbit_salcs`, `basis/salcbasis.jl`): the block
#   keeps `canon(N1)`, the complement of what vanishes or depends on blocks kept
#   earlier in the orbit;
# - the pointed basis FREEZES it, as upstream does (ported from SLCE.jl a248d4b): the
#   block is re-gauged as `[canon(N1) canon(N0)]` so that the vanishing part is columns
#   of its own, and `moment_resolvability`'s per-column test then names exactly those
#   columns, which every pointed fit holds at zero.
#
# In both cases the new columns are chosen in the coordinates of the block's old
# columns, which are already the crystal-frame gauge, so the choice does not depend on
# how the crystal is turned either.

# The relative band in which a singular value of a block's expansion would be neither a
# cancellation nor content. A cancellation lands at round-off of the gross mass (below
# `_CANCELLATION_RTOL`); genuine content sits at O(1) of it, whatever the orientation,
# because turning the crystal changes the expansion by an invertible map of bounded
# condition. A value between the two means the rank cannot be decided, and the build
# refuses instead of guessing. The test suite holds the measured margins on both sides
# (`test/unit/test_gauge.jl`).
const _REGAUGE_BAND_TOP = 1e-6

"""
    _block_split(S, gross) -> (; null, kept, ratios)

The rank decision of one block: `S` is the block's expansion (rows × block columns) and
`gross` the gross mass of each column. Returns the right singular vectors spanning the
vanishing subspace (`null`) and its complement (`kept`), and the singular values divided
by the block's gross mass (`ratios`, descending, padded with zeros to the block size).
A singular value is a cancellation when it is at most `_CANCELLATION_RTOL` times the
largest gross mass of the block.
"""
function _block_split(S::AbstractMatrix{Float64}, gross::AbstractVector{Float64})
    k = size(S, 2)
    mass = maximum(gross; init = 0.0)
    # QR first: the singular values and right vectors of `S` are those of its `R`
    # factor, and `R` is at most `k × k` whatever the number of rows.
    R = size(S, 1) == 0 ? zeros(Float64, 0, k) : Matrix(qr(S).R)
    F = svd(R; full = true)
    sv = vcat(F.S, zeros(Float64, k - length(F.S)))
    ratios = mass > 0.0 ? sv ./ mass : zeros(Float64, k)
    vanish = [r <= _CANCELLATION_RTOL for r in ratios]
    return (; null = F.V[:, vanish], kept = F.V[:, .!vanish], ratios)
end

# The gate on a block's rank decision: refuse a singular value that is neither a
# cancellation nor content (see `_REGAUGE_BAND_TOP`).
function _refuse_undecidable(ratios::AbstractVector{Float64}, what)
    undecided = findall(r -> _CANCELLATION_RTOL < r <= _REGAUGE_BAND_TOP, ratios)
    isempty(undecided) || error(
        "cannot decide which combinations of the columns of block $what vanish on " *
        "this reference cell: relative singular values $(ratios[undecided]) lie " *
        "between the cancellation cut $(_CANCELLATION_RTOL) and $(_REGAUGE_BAND_TOP). " *
        "This is a bug or a cell far outside what the gauge was measured on; please " *
        "report the crystal and the spec")
    return nothing
end

# The new gauge of a pointed block in the coordinates of its old columns, or `nothing`
# when the block keeps its columns: no vanishing direction, or nothing else. Entries
# that are round-off against a unit vector are set to zero, so a block whose vanishing
# subspace was already spanned by single columns comes back as an exact permutation and
# its tensors are copied, not recombined.
function _vanishing_gauge(S::AbstractMatrix{Float64}, gross::AbstractVector{Float64},
                          what)::Union{Nothing,Matrix{Float64}}
    split = _block_split(S, gross)
    _refuse_undecidable(split.ratios, what)
    k = size(S, 2)
    n0 = size(split.null, 2)
    (n0 == 0 || n0 == k) && return nothing
    W = hcat(_canonical_basis(split.kept), _canonical_basis(split.null))
    for i in eachindex(W)
        abs(W[i]) < _ROUNDOFF_RTOL && (W[i] = 0.0)
    end
    W == Matrix{Float64}(I, k, k) && return nothing
    return W
end

# The SALC `Σᵢ w[i]·columns[i]` under `key`. The columns of one block share their labels,
# so the combination is the term-by-term sum of their member tensors; the canonical
# member fold does the summing. Round-off of the sums is set to zero the way the
# projection does it: an entry against its own tensor, a whole tensor against the
# column (whose coefficient vector is a unit vector).
function _combine_salcs(columns::Vector{SALC}, w::AbstractVector{Float64},
                        key::SALCKey)::SALC
    raw = SALCMember[]
    for (i, s) in enumerate(columns)
        w[i] == 0.0 && continue
        for m in s.members
            push!(raw, SALCMember(m.atoms, m.shifts,
                                  [SALCTerm(t.slots, w[i] .* t.folded) for t in m.terms]))
        end
    end
    folded = _canonicalize_members(raw)
    big = maximum((maximum(abs, t.folded) for m in folded for t in m.terms); init = 0.0)
    members = SALCMember[]
    for m in folded
        terms = SALCTerm[]
        for t in m.terms
            maximum(abs, t.folded) > _ROUNDOFF_RTOL * big || continue
            push!(terms, SALCTerm(t.slots, _zero_roundoff!(copy(t.folded))))
        end
        isempty(terms) || push!(members, SALCMember(m.atoms, m.shifts, terms))
    end
    return SALC(key, key.body, key.decors, key.L_S, key.Lf, members)
end

# A tie is necessary for a vanishing direction: rows are keyed by per-atom content, so
# a block whose members all reach different atom multisets cannot cancel.
function _block_has_tie(columns::Vector{SALC})::Bool
    seen = Set{Vector{Int}}()
    for s in columns
        empty!(seen)
        for m in s.members
            atoms = sort(m.atoms)
            atoms in seen && return true
            push!(seen, atoms)
        end
    end
    return false
end

# The block of a column: its key without `block`.
_block_key(s::SALC) = (s.key.body, s.key.orbit_id, s.key.decors, s.key.L_S, s.key.Lf)

# No refinement of the block: the energy basis.
_no_class(::SALC) = nothing

# The mark class of a pointed column (the `salc_groups(::MomentBasis)` fingerprint): the
# atoms and the sites that carry the mark in the representative member. Columns of one
# energy-side block with different mark classes predict different atoms, so they are
# never mixed.
function _mark_class(s::SALC)
    mem = s.members[1]
    marks = Int[]
    sites = Int[]
    for t in mem.terms
        mark = findfirst(sl -> sl.factor.channel == DISP, t.slots)
        mark === nothing && continue
        push!(sites, t.slots[mark].site)
        push!(marks, mem.atoms[t.slots[mark].site])
    end
    return (sort!(unique!(marks)), sort!(unique!(sites)))
end

# The pointed signature expansion of a list of SALCs, densely: `S[r, j]` is the weight
# of row function `r` (mark atom, the mark's own `(l, m)`, the sorted environment
# symbols — `moment_resolvability`'s rows, with the mark atom that names its blocks
# folded into the key) in column `j`, `gross[j]` the absolute mass accumulated into it.
# The basis build calls it on the columns of one block (`_regauge_vanishing`).
function _moment_signature_of(sal::AbstractVector{SALC})
    rows = Dict{Tuple{Int,_MomentRowKey},Int}()
    entries = Vector{Vector{Tuple{Int,Float64}}}()   # per row: (col, weight)
    gross = zeros(length(sal))
    for (j, s) in enumerate(sal)
        scale = (4π)^(count(has_spin, s.decors) / 2)
        for mem in s.members, t in mem.terms
            slots = t.slots
            mslot = findfirst(sl -> sl.factor.channel == DISP, slots)
            mslot === nothing && error("pointed SALC without a mark slot")
            mark_site = slots[mslot].site
            mark_atom = mem.atoms[mark_site]
            for index in CartesianIndices(t.folded)
                w = t.folded[index] * scale
                w == 0.0 && continue
                mark_lm = (0, 0)
                env = Tuple{Int,Int,Int}[]
                for i in eachindex(slots)
                    sl = slots[i]
                    if sl.factor.channel == DISP
                        continue                     # the mark's |u|²R₀₀ ≡ const
                    elseif sl.site == mark_site
                        mark_lm = (sl.factor.l, index[i] - sl.factor.l - 1)
                    else
                        push!(env, (mem.atoms[sl.site], sl.factor.l,
                                    index[i] - sl.factor.l - 1))
                    end
                end
                sort!(env)
                key = (mark_atom, _MomentRowKey(mark_lm, env))
                r = get!(rows, key) do
                    push!(entries, Tuple{Int,Float64}[])
                    length(entries)
                end
                push!(entries[r], (j, w))
                gross[j] += abs(w)
            end
        end
    end
    S = zeros(length(entries), length(sal))
    for (r, es) in enumerate(entries), (j, w) in es
        S[r, j] += w
    end
    return S, gross
end

# The pointed expansion as the build calls it, or `nothing` for a block it cannot
# classify — the same refusal as `moment_resolvability`'s: two spin factors on one
# reference-cell atom (two images of one neighbour, or an environment on an image of
# the mark). The build leaves such a block alone, and the door still refuses it by name.
function _moment_block_signature(columns::Vector{SALC})
    for s in columns, m in s.members, t in m.terms
        mark = findfirst(sl -> sl.factor.channel == DISP, t.slots)
        mark === nothing && return nothing
        mark_site = t.slots[mark].site
        atoms = Int[m.atoms[mark_site]]
        for sl in t.slots
            (sl.factor.channel == SPIN && sl.site != mark_site) || continue
            push!(atoms, m.atoms[sl.site])
        end
        allunique(atoms) || return nothing
    end
    return _moment_signature_of(columns)
end

"""
    _regauge_vanishing(columns, signature, class) -> Vector{SALC}

The columns of one orbit with every block whose vanishing subspace is a proper,
nonzero part of it re-gauged as `[canon(N1) canon(N0)]` (see the head of this file).
`signature` expands a list of columns (`_moment_block_signature`), `class` refines the
block (`_mark_class`). The re-gauged columns of a block take the block's `block`
indices in ascending order, so the keys of the orbit do not change; only the functions
behind the keys of a re-gauged block do. Blocks without a tie are not expanded at all,
so a cell with unique minimum images builds exactly the basis it built without this
step.
"""
function _regauge_vanishing(columns::Vector{SALC}, signature::Function,
                            class::Function)::Vector{SALC}
    length(columns) < 2 && return columns
    groups = Dict{Any,Vector{Int}}()
    order = Any[]
    for (j, s) in enumerate(columns)
        g = (_block_key(s), class(s))
        haskey(groups, g) || push!(order, g)
        push!(get!(groups, g, Int[]), j)
    end
    out = copy(columns)
    for g in order
        index = sort(groups[g]; by = j -> columns[j].key.block)
        length(index) < 2 && continue
        block = columns[index]
        _block_has_tie(block) || continue
        expansion = signature(block)
        expansion === nothing && continue
        S, gross = expansion
        W = _vanishing_gauge(S, gross, block[1].key)
        W === nothing && continue
        for (b, j) in enumerate(index)
            out[j] = _combine_salcs(block, view(W, :, b), columns[j].key)
        end
    end
    return out
end

"""
    regauged_columns(basis::SCEBasis) -> Vector{Int}
    regauged_columns(mb::MomentBasis) -> Vector{Int}

The columns, in ascending order, of every block whose gauge this reference cell chose:
a block (the columns sharing every label of their [`SALCKey`](@ref) but `block`; in a
pointed basis also one mark class) of which a proper, nonzero part vanishes
identically on every configuration the cell can express, or depends on functions the
orbit already has.

- In an [`SCEBasis`](@ref) that part is DROPPED at the build (`build_salc_basis`
  names it in its warning) and the block keeps the complement, chosen in the
  coordinates of its own columns along the crystal axes; when the complement is not
  spanned by single columns, the kept columns are combinations of them. The list is
  recorded at the build and persisted with the basis.
- In a [`MomentBasis`](@ref) the part is kept as columns of its own, which
  [`moment_resolvability`](@ref) names and every pointed fit freezes; the remaining
  columns of the block are the complement. A basis the classification cannot read
  reports none.

These columns depend on the cell — the same crystal on a (3,3,3) supercell, or with
other lattice vectors, has other columns in that block — so their coefficients are not
comparable column by column across cells. Compare models through readouts that do not
depend on the gauge (`bilinear_terms`, predicted energies and torques).
[`coeftable`](@ref) carries the same flag as its `regauged` column.
"""
regauged_columns(basis::SCEBasis)::Vector{Int} = copy(basis.salc_basis.regauged)

function regauged_columns(mb::MomentBasis)::Vector{Int}
    vanishing = try
        moment_resolvability(mb).vanishing
    catch err
        err isa UnclassifiableBasis || rethrow()
        return Int[]
    end
    return _regauged(salcs(mb), vanishing, _mark_class)
end

function _regauged(sal::AbstractVector{SALC}, vanishing::AbstractVector{Int},
                   class::Function)::Vector{Int}
    isempty(vanishing) && return Int[]
    groups = Dict{Any,Vector{Int}}()
    for (j, s) in enumerate(sal)
        push!(get!(groups, (_block_key(s), class(s)), Int[]), j)
    end
    vset = Set(vanishing)
    out = Int[]
    for cols in values(groups)
        n0 = count(in(vset), cols)
        0 < n0 < length(cols) && append!(out, cols)
    end
    return sort!(out)
end
