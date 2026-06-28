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
    evaluate_salc(salc, e) -> Float64

Evaluate `Φ(e)` for a single (super)cell spin configuration `e` (`3 × n_atoms`,
unit columns). Periodic with the cell, so lattice shifts map back to the same
column (`site(a, R) → a`).
"""
function evaluate_salc(salc::SALC, e::AbstractMatrix{<:Real})::Float64
    N = salc.body
    scale = (4π)^(N / 2)
    total = 0.0
    @inbounds for m in salc.members
        atoms = m.atoms
        for t in m.terms
            ls = t.ls
            for idx in CartesianIndices(t.folded)
                w = t.folded[idx]
                w == 0.0 && continue
                for i = 1:N
                    μ = idx[i] - ls[i] - 1
                    u = SVector{3,Float64}(e[1, atoms[i]], e[2, atoms[i]], e[3, atoms[i]])
                    # μ is in range by the tensor bounds and u is a unit column by the
                    # config contract, so the unchecked variant is safe here.
                    w *= Harmonics.Zlm_unsafe(ls[i], μ, u)
                end
                total += w
            end
        end
    end
    return scale * total
end

"""
    accumulate_grad!(G, salc, e, weight) -> G

Accumulate `weight · ∇Φ(e)` into `G` (a `3 × n_atoms` buffer), where `∇Φ` is the
per-site direction gradient of the SALC orbit sum: column `a` of `G` receives
`weight · ∂Φ/∂e_a`. Summing this over a model's SALCs (with `weight = jϕ`) yields
`Σ_ϕ jϕ ∇Φ_ϕ`, whose cross product with the spins gives the torque
`τ_a = e_a × ∂E/∂e_a`.

The gradient distributes over the product rule: for each cluster member, each
folded-tensor multi-index `μ`, and each site `i` of the member,
`∂/∂e_{aᵢ} ∏ₖ Zₗₖμₖ = (∏_{k≠i} Zₗₖμₖ) ∇Zₗᵢμᵢ`, landing on column `aᵢ`. A site that
appears more than once in a member contributes once per occurrence (the chain
rule), so repeated sites are handled correctly. `∇Zₗₘ` is the tangent-projected
gradient; the radial part it drops would cancel in `e × ∇Φ` anyway.
"""
function accumulate_grad!(G::AbstractMatrix{Float64}, salc::SALC,
                          e::AbstractMatrix{<:Real}, weight::Real)
    weight == 0.0 && return G
    N = salc.body
    scale = weight * (4π)^(N / 2)
    @inbounds for m in salc.members
        atoms = m.atoms
        for t in m.terms
            ls = t.ls
            for idx in CartesianIndices(t.folded)
                w = scale * t.folded[idx]
                w == 0.0 && continue
                for i = 1:N
                    # leave-one-out product of the other sites' harmonics
                    p = 1.0
                    for k = 1:N
                        k == i && continue
                        μk = idx[k] - ls[k] - 1
                        uk = SVector{3,Float64}(e[1, atoms[k]], e[2, atoms[k]], e[3, atoms[k]])
                        p *= Harmonics.Zlm_unsafe(ls[k], μk, uk)
                    end
                    p == 0.0 && continue
                    μi = idx[i] - ls[i] - 1
                    ui = SVector{3,Float64}(e[1, atoms[i]], e[2, atoms[i]], e[3, atoms[i]])
                    gi = Harmonics.grad_Zlm_unsafe(ls[i], μi, ui)
                    c = w * p
                    a = atoms[i]
                    G[1, a] += c * gi[1]
                    G[2, a] += c * gi[2]
                    G[3, a] += c * gi[3]
                end
            end
        end
    end
    return G
end
