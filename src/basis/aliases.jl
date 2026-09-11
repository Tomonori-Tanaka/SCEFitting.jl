# Cross-orbit alias groups: distinct cluster orbits whose members join the SAME
# reference-cell atom sets through periodic images the space group does not relate
# (a Wigner–Seitz boundary tie the point group fails to fuse — the second face of
# periodic resolvability, see docs/src/theory/resolvability.md). Under cell-periodic
# evaluation every image of an atom carries the same spin, so such orbits evaluate to
# proportional design columns: the data determine only the SUM of their couplings.
#
# This package ties the columns of an alias group into ONE design column and reads
# the fitted coefficient back to every member orbit with equal per-bond weight — a
# recorded CONVENTION, not a measurement (the cell samples J(q) only on its
# commensurate mesh, and the tied images alias exactly there). Detection is
# structural (no training data involved) and runs once, at basis construction, on
# the same aggregated shift-blind monomial vectors `_reduce_orbit_salcs` uses within
# an orbit. Nothing is merged or dropped in the basis itself: the orbits, their
# `SALCKey`s and their members stay distinct, so a future datum type that breaks the
# periodicity (a spin spiral) could separate them on the same basis object.

"""
Relative residual below which two orbits' aggregated function vectors count as
proportional (`‖v_b − c·v_a‖ / ‖v_b‖`, `c` the least-squares factor). Measured on
Nd₂Fe₁₄B 1×1×1 (`P4₂/mnm`, 68 atoms, isotropic pairs to the WS boundary): the ten
genuine alias pairs sit at exactly `0` (six) and `1.5e-8` (four — gauge-fixing
round-off of the two independently projected representatives), while the smallest
NON-alias same-channel residual is `O(1e-2)` (the MnTe survivors of the per-orbit
reduction). `1e-6` has two decades above the measured aliases and four below the
first survivor. It is deliberately looser than `_AGG_DEP_RTOL` (`1e-8`), which
would miss the four round-off-limited pairs.
"""
const _ALIAS_RTOL = 1e-6

# Hard cap on `alias_rtol`, the analogue of `_TIE_TOL_MAX`: past 1e-2 the band reaches
# the smallest genuine same-channel residual and would impose an equal split between
# couplings the cell CAN tell apart.
const _ALIAS_RTOL_MAX = 1e-2

# Snap band for the per-bond weights: transported copies of one representative have
# identical per-member tensor norms up to round-off, and the tie must then be an
# EXACT ±1 (so the folded column is the exact sum and the read-out the exact
# coefficient), not 1 ± 1e-16.
const _ALIAS_WEIGHT_SNAP = 1e-10

# The admissible values of `SCEPredictor.split` / the persisted `split` field.
const _SPLIT_KINDS = (:free, :convention, :legacy)

"""
    AliasGroup

One cross-orbit alias group of an [`SCEBasis`](@ref): SALCs (design-order indices
`salcs`, from ≥ 2 distinct `orbit_ids`) whose aggregated functions are proportional
on the reference cell, so that cell-periodic data determine only the sum of their
couplings. Read them from [`alias_groups`](@ref).

- `body`, `atoms`: body order and the representative member's atoms (of `salcs[1]`).
- `delta_shifts[k]`: per-site lattice-vector offset of orbit `k`'s matching member
  relative to orbit 1's — the images that alias. An axis with a nonzero component in
  any of these is one the training cell must be doubled along to break the tie.
- `distance`: the pair separation in Å (`NaN` for `body ≥ 3`).
- `weights[k]`: the per-bond weight of `salcs[k]` in the tied column and in the
  read-out (`jϕ_salc = weights[k] · jϕ_column`); `±1` for transported copies of one
  representative tensor, `NaN` when `kind == :span_collapsed`.
- `kind`: `:proportional` (tied as one column); `:span_collapsed` (the orbits'
  functions are linearly dependent without being pairwise proportional — an
  anisotropic channel whose member tensors carry their own bond geometry); or
  `:unequal_norm` (proportional, but the per-member tensor norms differ, so the
  transported-copy premise behind the per-bond weights fails). The last two are
  reported only, never tied: the design stays rank deficient there and the `OLS`
  rank warning applies.
"""
struct AliasGroup
    salcs::Vector{Int}
    orbit_ids::Vector{Int}
    body::Int
    atoms::Vector{Int}
    delta_shifts::Vector{Vector{SVector{3,Int}}}
    distance::Float64
    weights::Vector{Float64}
    kind::Symbol
end

# The SALC → design-column map of a basis. `col_of[j]` / `weight[j]` place SALC `j`
# into its column with its per-bond weight; `columns[t]` lists the SALCs of column
# `t` (a singleton for every untied SALC). `trivial` ⇔ no group is tied, in which
# case every consumer takes the identity fast path and the design is byte-identical
# to a basis built without alias detection.
struct _ColumnTies
    col_of::Vector{Int}
    weight::Vector{Float64}
    columns::Vector{Vector{Int}}
    groups::Vector{AliasGroup}
    trivial::Bool
end

_trivial_ties(p::Int, groups::Vector{AliasGroup} = AliasGroup[]) =
    _ColumnTies(collect(1:p), ones(Float64, p), [[j] for j = 1:p], groups, true)

# Minimal union–find over 1:n (path halving, union by index).
function _uf_find(parent::Vector{Int}, i::Int)::Int
    while parent[i] != i
        parent[i] = parent[parent[i]]
        i = parent[i]
    end
    return i
end
function _uf_union!(parent::Vector{Int}, a::Int, b::Int)
    ra = _uf_find(parent, a)
    rb = _uf_find(parent, b)
    ra == rb && return nothing
    parent[max(ra, rb)] = min(ra, rb)
    return nothing
end

# The structural signature two SALCs must share before their functions are compared:
# same channel `(body, decors, L_S, Lf)` and the same set of member atom multisets.
# Distinct orbits with equal signatures are exactly the "one atom multiset reached
# from two orbits" condition of a non-fused tie; equal signatures within ONE orbit
# are the ordinary `block` / `Lf` siblings the per-orbit reduction already made
# independent, and are never compared here.
function _alias_signature(s::SALC)
    atomsets = unique!(sort!([sort(m.atoms) for m in s.members]))
    return (s.body, s.decors, s.L_S, s.Lf, atomsets)
end

# Signature classes with ≥ 2 distinct orbits, each as a sorted SALC index list.
# Classes containing a displacement-decorated SALC are skipped (the aggregate key of
# `_function_vector` names a function only on pure-spin SALCs), and so are classes
# with a repeated-atom member (`AllImages` self-images: a tiling template that is
# never fitted on this cell, whose distinct images are the whole point).
function _alias_classes(salcs::Vector{SALC})::Vector{Vector{Int}}
    classes = Dict{Any,Vector{Int}}()
    for j in eachindex(salcs)
        push!(get!(classes, _alias_signature(salcs[j]), Int[]), j)
    end
    out = Vector{Int}[]
    for members in values(classes)
        length(unique(salcs[j].key.orbit_id for j in members)) >= 2 || continue
        all(j -> all(is_pure_spin, salcs[j].decors), members) || continue
        all(j -> all(m -> allunique(m.atoms), salcs[j].members), members) || continue
        push!(out, sort!(members))
    end
    sort!(out; by = first)
    return out
end

# `v_b − c·v_a` for the least-squares factor `c = ⟨v_a, v_b⟩ / ‖v_a‖²`, as a relative
# residual, plus `c`. Computed by explicit vector subtraction, not `1 − cos²`: at the
# measured residual scale (1e-8) the cosine form would cancel to round-off.
function _proportional_residual(va::_FunctionVector, na::Float64,
                                vb::_FunctionVector, nb::Float64)
    c = _dictdot(va, vb) / (na * na)
    r2 = 0.0
    for (k, x) in vb
        r2 += (x - c * get(va, k, 0.0))^2
    end
    for (k, x) in va
        haskey(vb, k) || (r2 += (c * x)^2)
    end
    return sqrt(r2) / nb, c
end

# The aggregated function vector, its norm, and the unaggregated tensor norm of
# every SALC in `members`.
function _class_vectors(salcs::Vector{SALC}, members::Vector{Int})
    vecs = Dict{Int,Tuple{_FunctionVector,Float64,Float64}}()
    for j in members
        v, raw = _function_vector(salcs[j])
        vecs[j] = (v, _dictnorm(v), raw)
    end
    return vecs
end

# Union–find over one class: SALCs of DIFFERENT orbits whose vectors are proportional
# within `rtol` end up in one set. Returns the sets (sorted index lists, ≥ 1 each).
function _proportional_sets(salcs::Vector{SALC}, members::Vector{Int}, vecs, rtol::Float64)
    slot = Dict(j => i for (i, j) in enumerate(members))
    parent = collect(1:length(members))
    for a in members, b in members
        a < b || continue
        salcs[a].key.orbit_id == salcs[b].key.orbit_id && continue
        va, na, _ = vecs[a]
        vb, nb, _ = vecs[b]
        (na > 0 && nb > 0) || continue
        residual, c = _proportional_residual(va, na, vb, nb)
        (residual <= rtol && c != 0.0) && _uf_union!(parent, slot[a], slot[b])
    end
    sets = Dict{Int,Vector{Int}}()
    for j in members
        push!(get!(sets, _uf_find(parent, slot[j]), Int[]), j)
    end
    return sort!([sort!(s) for s in values(sets)]; by = first)
end

# Per-bond weights of a proportional set relative to its lowest-index member:
# `w_j = sign(c_j) · ν_ref / ν_j` with `ν_j = raw_j / √n_members(j)` the per-member
# tensor norm (`raw` is `_function_vector`'s unaggregated norm, members counted
# after canonicalization). Transported copies of one representative give `|w| = 1`
# up to round-off, snapped to an exact ±1; a weight that does not snap is evidence
# that the equal-norm premise failed — the caller then reports the set as
# `:unequal_norm` and leaves it untied.
function _set_weights(salcs::Vector{SALC}, members::Vector{Int}, vecs)::Vector{Float64}
    ref = members[1]
    vref, nref, rawref = vecs[ref]
    nuref = rawref / sqrt(length(salcs[ref].members))
    weights = Float64[]
    for j in members
        v, n, raw = vecs[j]
        _, c = _proportional_residual(vref, nref, v, n)
        w = sign(c) * nuref / (raw / sqrt(length(salcs[j].members)))
        abs(abs(w) - 1) < _ALIAS_WEIGHT_SNAP && (w = sign(c))
        push!(weights, w)
    end
    return weights
end

# Whether the representatives of a class (one SALC per proportional set) are linearly
# dependent as functions — the span-collapse case, detected by the same modified
# Gram–Schmidt pass `_reduce_orbit_salcs` runs within an orbit, here across orbits.
function _span_collapsed(reps::Vector{Int}, vecs, rtol::Float64)::Bool
    length(reps) >= 2 || return false
    kept = _FunctionVector[]
    for j in reps
        v = copy(vecs[j][1])
        nv, nv0 = _mgs_residual!(v, kept)
        nv <= max(rtol, _AGG_DEP_RTOL) * nv0 && return true
        for k in keys(v)
            v[k] /= nv
        end
        push!(kept, v)
    end
    return false
end

# Per-orbit lattice-vector offsets of the matching members: for each orbit in the
# group, the member with the same atoms as orbit 1's first member, and the per-site
# difference of its (re-anchored) shifts. Two orbits with the same atom sets always
# have such a member (the signature class guarantees it).
function _delta_shifts(salcs::Vector{SALC}, members::Vector{Int})
    ref = salcs[members[1]].members[1]
    out = Vector{SVector{3,Int}}[]
    for j in members
        k = findfirst(mem -> mem.atoms == ref.atoms, salcs[j].members)
        k === nothing && error("internal: alias-group SALC $j has no member on the " *
                               "atoms $(ref.atoms) its signature class shares")
        mem = salcs[j].members[k]
        push!(out, [mem.shifts[i] - ref.shifts[i] for i in eachindex(ref.atoms)])
    end
    return out
end

function _pair_distance(crystal::Crystal, m::SALCMember)::Float64
    length(m.atoms) == 2 || return NaN
    A = crystal.lattice.vectors
    fr = crystal.frac_positions
    d = A * (SVector{3,Float64}(fr[:, m.atoms[2]]) + m.shifts[2] -
             SVector{3,Float64}(fr[:, m.atoms[1]]) - m.shifts[1])
    return norm(d)
end

function _alias_group(salcs::Vector{SALC}, crystal::Crystal, members::Vector{Int},
                      weights::Vector{Float64}, kind::Symbol)::AliasGroup
    lead = salcs[members[1]]
    return AliasGroup(copy(members), [salcs[j].key.orbit_id for j in members], lead.body,
                      copy(lead.members[1].atoms), _delta_shifts(salcs, members),
                      _pair_distance(crystal, lead.members[1]), weights, kind)
end

# The SALC → column map for the tied groups: columns in order of their first SALC,
# untied SALCs as singletons.
function _column_map(p::Int, tied::Vector{AliasGroup})
    col_of = zeros(Int, p)
    weight = ones(Float64, p)
    group_of = Dict{Int,Int}()
    for (g, grp) in enumerate(tied), (k, j) in enumerate(grp.salcs)
        group_of[j] = g
        weight[j] = grp.weights[k]
    end
    columns = Vector{Int}[]
    placed = falses(length(tied))
    for j = 1:p
        g = get(group_of, j, 0)
        if g == 0
            push!(columns, [j])
            col_of[j] = length(columns)
        elseif !placed[g]
            push!(columns, copy(tied[g].salcs))
            placed[g] = true
            for jj in tied[g].salcs
                col_of[jj] = length(columns)
            end
        end
    end
    return col_of, weight, columns
end

"""
    _column_ties(salc_basis, crystal; alias_rtol = _ALIAS_RTOL) -> _ColumnTies

Detect the cross-orbit alias groups of `salc_basis` and build its SALC → design-column
map. Within each structural signature class (same channel, same member atom sets,
≥ 2 orbits) the aggregated shift-blind function vectors (`_function_vector`) are
compared pairwise across orbits; pairs whose relative residual is `≤ alias_rtol` are
united into a `:proportional` group and tied into one column with per-bond weights
(`_set_weights`). A class whose representatives are linearly dependent without being
pairwise proportional is reported as one `:span_collapsed` group and left untied.
`alias_rtol = nothing` disables detection (every SALC its own column).
"""
function _column_ties(sb::SALCBasis, crystal::Crystal;
                      alias_rtol::Union{Nothing,Real} = _ALIAS_RTOL)::_ColumnTies
    salcs = sb.salcs
    p = length(salcs)
    alias_rtol === nothing && return _trivial_ties(p)
    rtol = Float64(alias_rtol)
    (isfinite(rtol) && rtol >= 0) || throw(ArgumentError(
        "alias_rtol must be finite and ≥ 0 (or nothing); got $alias_rtol"))
    rtol < _ALIAS_RTOL_MAX || throw(ArgumentError(
        "alias_rtol = $alias_rtol is at or above the hard cap $_ALIAS_RTOL_MAX: a " *
        "band that wide reaches the smallest genuine same-channel residual and would " *
        "tie couplings the cell can tell apart"))
    groups = AliasGroup[]
    for members in _alias_classes(salcs)
        vecs = _class_vectors(salcs, members)
        sets = _proportional_sets(salcs, members, vecs, rtol)
        for s in sets
            length(s) >= 2 || continue
            w = _set_weights(salcs, s, vecs)
            kind = all(x -> abs(x) == 1.0, w) ? :proportional : :unequal_norm
            push!(groups, _alias_group(salcs, crystal, s, w, kind))
        end
        if _span_collapsed([s[1] for s in sets], vecs, rtol)
            push!(groups, _alias_group(salcs, crystal, members, fill(NaN, length(members)),
                                       :span_collapsed))
        end
    end
    sort!(groups; by = g -> g.salcs[1])
    tied = filter(g -> g.kind === :proportional, groups)
    isempty(tied) && return _trivial_ties(p, groups)
    col_of, weight, columns = _column_map(p, tied)
    return _ColumnTies(col_of, weight, columns, groups, false)
end

# Fold SALC columns into design columns: `X_t[:, t] = Σ_{j ∈ columns[t]} w_j · X[:, j]`.
# Identity (the same matrix object) when no group is tied, so an alias-free basis
# keeps its design byte-identical. For the `±1` weights of transported copies the
# fold is an exact sum.
function _fold_columns(X::Matrix{Float64}, ties::_ColumnTies)::Matrix{Float64}
    ties.trivial && return X
    n = size(X, 1)
    Xt = zeros(Float64, n, length(ties.columns))
    @inbounds for j = 1:size(X, 2)
        t = ties.col_of[j]
        w = ties.weight[j]
        if w == 1.0
            @views Xt[:, t] .+= X[:, j]
        else
            @views Xt[:, t] .+= w .* X[:, j]
        end
    end
    return Xt
end

# Expand a column-space coefficient vector to SALC space:
# `jϕ_salc[j] = w_j · jϕ_col[col_of[j]]`.
function _expand_coefficients(jphi_col::AbstractVector{<:Real},
                              ties::_ColumnTies)::Vector{Float64}
    length(jphi_col) == length(ties.columns) || throw(DimensionMismatch(
        "got $(length(jphi_col)) column coefficients for $(length(ties.columns)) " *
        "design columns"))
    ties.trivial && return collect(Float64, jphi_col)
    return Float64[ties.weight[j] * jphi_col[ties.col_of[j]] for j in eachindex(ties.col_of)]
end

# Per-SALC split status of a fitted model: `:convention` inside a tied group, `:free`
# elsewhere (see `SCEPredictor.split`).
function _split_status(ties::_ColumnTies)::Vector{Symbol}
    status = fill(:free, length(ties.col_of))
    for g in ties.groups
        g.kind === :proportional || continue
        for j in g.salcs
            status[j] = :convention
        end
    end
    return status
end

# Per-SALC index into `ties.groups` (0 = not in any alias group), the `alias_group`
# column of `coeftable`.
function _alias_index(ties::_ColumnTies)::Vector{Int}
    index = zeros(Int, length(ties.col_of))
    for (g, grp) in enumerate(ties.groups), j in grp.salcs
        # a tied SALC keeps pointing at its tie when the same class is also
        # reported as span-collapsed
        (index[j] != 0 && grp.kind !== :proportional) && continue
        index[j] = g
    end
    return index
end

# One line per group for the construction-time report: the offsets of every orbit
# after the first, one bracket per site.
function _alias_group_line(crystal::Crystal, g::AliasGroup)::String
    labels = crystal.species_labels[crystal.species]
    atoms = join(("$(labels[a])$a" for a in g.atoms), "–")
    orbits = join(string.(g.orbit_ids), "/")
    site(s) = "(" * join(string.(s), ",") * ")"
    offsets = join((join((site(s) for s in d), "") for d in g.delta_shifts[2:end]), " ")
    dist = isnan(g.distance) ? "" : ", d = $(round(g.distance; digits = 3)) Å"
    tag = g.kind === :proportional ? "tied" :
          g.kind === :span_collapsed ? "span collapsed, NOT tied" :
          "unequal per-member norms $(round.(g.weights; sigdigits = 4)), NOT tied"
    return "  $(g.body)-body orbits $orbits: $atoms$dist, image offsets $offsets [$tag]"
end

function _report_alias_groups(crystal::Crystal, ties::_ColumnTies)
    isempty(ties.groups) && return nothing
    tied = [g for g in ties.groups if g.kind === :proportional]
    collapsed = [g for g in ties.groups if g.kind !== :proportional]
    if !isempty(tied)
        tie_axes = Set{Int}()
        for g in tied, d in g.delta_shifts, s in d, k = 1:3
            s[k] == 0 || push!(tie_axes, k)
        end
        axis_names = join((("a", "b", "c")[k] for k in sort!(collect(tie_axes))), ", ")
        @info "cross-orbit alias groups: $(length(tied)) group(s) of distinct cluster " *
              "orbits join the same reference-cell atoms through periodic images the " *
              "space group does not relate (Wigner–Seitz boundary ties). Cell-periodic " *
              "data determine only the SUM of each group's couplings: the columns are " *
              "tied into one design column and the fitted coefficient is read back to " *
              "every orbit with equal per-bond weight — a recorded convention " *
              "(`coeftable` column `split = :convention`), not a measurement. It leaves " *
              "training-cell energies, torques and every q = 0 quantity unchanged, but " *
              "J(q) off Γ, the tiled-supercell energy of configurations that are not " *
              "periodic in this cell, and any bond-resolved J(R) table depend on it. To " *
              "MEASURE the split, break the tie: double the cell along $axis_names " *
              "(every axis with a nonzero image offset below), or add spin-spiral " *
              "data.\n" * join((_alias_group_line(crystal, g) for g in tied), "\n")
    end
    if !isempty(collapsed)
        @warn "cross-orbit alias groups NOT tied: $(length(collapsed)) class(es) of " *
              "distinct orbits share their atom sets but either span a linearly " *
              "dependent function space without being pairwise proportional (an " *
              "anisotropic channel whose member tensors carry their own bond geometry) " *
              "or have unequal per-member tensor norms. No equal per-bond split is " *
              "defined for them: the design stays rank deficient there and the OLS rank " *
              "warning applies. Break the tie with a doubled cell.\n" *
              join((_alias_group_line(crystal, g) for g in collapsed), "\n")
    end
    return nothing
end
