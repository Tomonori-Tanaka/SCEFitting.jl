"""
    SALCKey

Canonical structural address of a SALC (design-matrix column): a tuple of
integers that identifies the column independent of construction order. The
`(body, orbit_id, ls, Lf)` part is gauge-invariant; `block` indexes within a
degenerate invariant subspace (made reproducible by the canonical gauge).
"""
struct SALCKey
    body::Int
    orbit_id::Int
    ls::Vector{Int}
    Lf::Int
    block::Int
end

_keytuple(k::SALCKey) = (k.body, k.orbit_id, k.ls, k.Lf, k.block)
Base.isless(a::SALCKey, b::SALCKey) = _keytuple(a) < _keytuple(b)
Base.:(==)(a::SALCKey, b::SALCKey) = _keytuple(a) == _keytuple(b)
Base.hash(k::SALCKey, h::UInt) = hash(_keytuple(k), h)

"""
    SALCTerm

One `l`-ordering contributing to a (member of a) SALC: per-site angular momenta
`ls` and the transported real coefficient tensor `folded` (rank `N`, `Mf` axis
already contracted against the SALC coefficient). A SALC built from an `l`-multiset
with unequal `l`'s on symmetry-equivalent sites carries one term per ordering; the
common case (1/2-body, or equal `l`'s) is a single term.
"""
struct SALCTerm
    ls::Vector{Int}
    folded::Array{Float64}
end

"""
    SALCMember

One concrete cluster instance contributing to a SALC's orbit sum: `atoms` and
per-site lattice `shifts` (in the representative's site order), with this member's
list of transported [`SALCTerm`](@ref)s (one per `l`-ordering).
"""
struct SALCMember
    atoms::Vector{Int}
    shifts::Vector{SVector{3,Int}}
    terms::Vector{SALCTerm}
end

"""
    SALC

One symmetry-adapted, time-reversal-even scalar basis function over a cluster
orbit. The function is the orbit sum `Φ(e) = (4π)^(N/2) Σ_members Σ_terms Σ_μ
folded[μ] ∏ᵢ Z_{lsᵢ,μᵢ}(e_{site})`, where each member contributes one or more
`l`-ordering terms. `SALC.ls` is the sorted `l`-multiset label (the per-term `ls`
are its orderings). See [`build_salc_basis`](@ref).
"""
struct SALC
    key::SALCKey
    body::Int
    ls::Vector{Int}
    Lf::Int
    members::Vector{SALCMember}
end

"""
    SALCBasis

An ordered list of [`SALC`](@ref)s (sorted by key), the matching `keys` vector
(used to address design-matrix columns), and a structural `fingerprint`
(`hash` of the sorted keys; gauge-independent).
"""
struct SALCBasis
    salcs::Vector{SALC}
    keys::Vector{SALCKey}
    fingerprint::UInt64
end

Base.length(b::SALCBasis) = length(b.salcs)

"""
    evaluate_salc(salc, e[, cache]) -> Float64

Evaluate `Φ(e)` for a single (super)cell spin configuration `e` (`3 × n_atoms`,
unit columns). Periodic with the cell, so lattice shifts map back to the same
column (`site(a, R) → a`). Hot loops may pass `cache` (a reusable
`Vector{Float64}`, grown on demand) — the dnPl recursion workspace of the
cache-threaded `Harmonics.Zlm_unsafe`; values are identical with or without it.
"""
function evaluate_salc(salc::SALC, e::AbstractMatrix{<:Real},
                       cache::Vector{Float64} = Vector{Float64}(undef, 4))::Float64
    scale = (4π)^(salc.body / 2)
    total = 0.0
    @inbounds for m in salc.members
        atoms = m.atoms
        for t in m.terms
            total += _eval_term(t.folded, t.ls, atoms, e, cache)
        end
    end
    return scale * total
end

# Per-term kernel behind a function barrier: `SALCTerm.folded` is stored as the
# abstract `Array{Float64}` (rank not in the type), so the multi-index loop is
# type-unstable when written inline. Dispatching here specializes on the concrete
# rank `D` (== body order) — same values and the same multiply/accumulate order, so
# the result is bit-identical to the inlined loop.
@inline function _eval_term(folded::Array{Float64,D}, ls, atoms,
                            e::AbstractMatrix{<:Real}, cache::Vector{Float64}) where {D}
    # Per-site harmonic tables. The site direction `u_i` is fixed for this term, so
    # `Z_{lsᵢ,μ}` depends only on `(i, μ)`; tabulate it once over `μ ∈ -lsᵢ:lsᵢ` (the
    # tensor axis index is `idx[i] = μ+lsᵢ+1`) instead of recomputing inside the
    # nonzero multi-index loop, where each call hit `dnPl` and dominated the design-
    # matrix cost. Same values, same multiply order ⇒ bit-identical. `u_i` is a unit
    # column by the config contract, so the unchecked harmonic is safe.
    tabs = _site_ztables(Val(D), ls, atoms, e, cache)
    acc = 0.0
    @inbounds for idx in CartesianIndices(folded)
        w = folded[idx]
        w == 0.0 && continue
        for i = 1:D
            w *= tabs[i][idx[i]]
        end
        acc += w
    end
    return acc
end

# Tuple of per-site harmonic tables `tabs[i][μ+lsᵢ+1] = Z_{lsᵢ,μ}(u_i)`. `Val(D)`
# keeps the tuple length a compile-time constant (type-stable). `cache` is the
# reusable dnPl recursion workspace of the cache-threaded `Zlm_unsafe` (grown on
# demand; values are bit-identical to the cache-less call).
@inline function _site_ztables(::Val{D}, ls, atoms, e::AbstractMatrix{<:Real},
                               cache::Vector{Float64}) where {D}
    ntuple(Val(D)) do i
        a = atoms[i]
        u = SVector{3,Float64}(e[1, a], e[2, a], e[3, a])
        li = ls[i]
        length(cache) < li + 1 && resize!(cache, li + 1)
        Float64[Harmonics.Zlm_unsafe(li, μ, u, cache) for μ = -li:li]
    end
end

"""
    accumulate_grad!(G, salc, e, weight[, cache]) -> G

Accumulate `weight · ∇Φ(e)` into `G` (a `3 × n_atoms` buffer), where `∇Φ` is the
per-site direction gradient of the SALC orbit sum: column `a` of `G` receives
`weight · ∂Φ/∂e_a`. Summing this over a model's SALCs (with `weight = jϕ`) yields
`Σ_ϕ jϕ ∇Φ_ϕ`, whose cross product with the spins gives the torque
`τ_a = −e_a × ∂E/∂e_a = ∂E/∂e_a × e_a` (the physical / Landau–Lifshitz sign).

The gradient distributes over the product rule: for each cluster member, each
folded-tensor multi-index `μ`, and each site `i` of the member,
`∂/∂e_{aᵢ} ∏ₖ Zₗₖμₖ = (∏_{k≠i} Zₗₖμₖ) ∇Zₗᵢμᵢ`, landing on column `aᵢ`. A site that
appears more than once in a member contributes once per occurrence (the chain
rule), so repeated sites are handled correctly. `∇Zₗₘ` is the tangent-projected
gradient; the radial part it drops would cancel in `e × ∇Φ` anyway.
"""
function accumulate_grad!(G::AbstractMatrix{Float64}, salc::SALC,
                          e::AbstractMatrix{<:Real}, weight::Real,
                          cache::Vector{Float64} = Vector{Float64}(undef, 4))
    weight == 0.0 && return G
    scale = weight * (4π)^(salc.body / 2)
    @inbounds for m in salc.members
        atoms = m.atoms
        for t in m.terms
            _accum_grad_term!(G, t.folded, t.ls, atoms, e, scale, cache)
        end
    end
    return G
end

# Gradient counterpart of `_eval_term`: a function barrier specializing on the
# concrete rank `D`. The site harmonics `Z` and gradients `∇Z` are tabulated per
# site (same rationale as `_eval_term`: `u_i` is fixed for the term, so each entry
# depends only on `(i, μ)`), then the product-rule expansion reads them back. Same
# expansion, same evaluation order ⇒ the accumulated gradient is bit-identical.
@inline function _accum_grad_term!(G::AbstractMatrix{Float64}, folded::Array{Float64,D},
                                   ls, atoms, e::AbstractMatrix{<:Real}, scale::Float64,
                                   cache::Vector{Float64}) where {D}
    ztab = _site_ztables(Val(D), ls, atoms, e, cache)
    gtab = _site_gtables(Val(D), ls, atoms, e, cache)
    @inbounds for idx in CartesianIndices(folded)
        w = scale * folded[idx]
        w == 0.0 && continue
        for i = 1:D
            # leave-one-out product of the other sites' harmonics
            p = 1.0
            for k = 1:D
                k == i && continue
                p *= ztab[k][idx[k]]
            end
            p == 0.0 && continue
            gi = gtab[i][idx[i]]
            c = w * p
            a = atoms[i]
            G[1, a] += c * gi[1]
            G[2, a] += c * gi[2]
            G[3, a] += c * gi[3]
        end
    end
    return G
end

# Tuple of per-site gradient tables `gtab[i][μ+lsᵢ+1] = ∇Z_{lsᵢ,μ}(u_i)`.
@inline function _site_gtables(::Val{D}, ls, atoms, e::AbstractMatrix{<:Real},
                               cache::Vector{Float64}) where {D}
    ntuple(Val(D)) do i
        a = atoms[i]
        u = SVector{3,Float64}(e[1, a], e[2, a], e[3, a])
        li = ls[i]
        length(cache) < li + 1 && resize!(cache, li + 1)
        SVector{3,Float64}[Harmonics.grad_Zlm_unsafe(li, μ, u, cache) for μ = -li:li]
    end
end
