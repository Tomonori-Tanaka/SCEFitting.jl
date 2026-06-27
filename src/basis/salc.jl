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
    SALCMember

One concrete cluster instance contributing to a SALC's orbit sum: `atoms` and
per-site lattice `shifts` (in the representative's site order), with this member's
transported real coefficient tensor `folded` (rank `N`, `Mf` axis already
contracted against the SALC coefficient).
"""
struct SALCMember
    atoms::Vector{Int}
    shifts::Vector{SVector{3,Int}}
    folded::Array{Float64}
end

"""
    SALC

One symmetry-adapted, time-reversal-even scalar basis function over a cluster
orbit. The function is the orbit sum `Φ(e) = (4π)^(N/2) Σ_members Σ_μ
folded[μ] ∏ᵢ Zₗᵢ,μᵢ(e_{site})`. See [`build_salc_basis`](@ref).
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
    evaluate(salc, e) -> Float64

Evaluate `Φ(e)` for a single (super)cell spin configuration `e` (`3 × n_atoms`,
unit columns). Periodic with the cell, so lattice shifts map back to the same
column (`site(a, R) → a`).
"""
function evaluate(salc::SALC, e::AbstractMatrix{<:Real})::Float64
    N = salc.body
    scale = (4π)^(N / 2)
    total = 0.0
    @inbounds for m in salc.members
        acc = 0.0
        for idx in CartesianIndices(m.folded)
            w = m.folded[idx]
            w == 0.0 && continue
            for i = 1:N
                μ = idx[i] - salc.ls[i] - 1
                u = SVector{3,Float64}(e[1, m.atoms[i]], e[2, m.atoms[i]], e[3, m.atoms[i]])
                # μ is in range by the tensor bounds and u is a unit column by the
                # config contract, so the unchecked variant is safe here.
                w *= Harmonics.Zlm_unsafe(salc.ls[i], μ, u)
            end
            acc += w
        end
        total += acc
    end
    return scale * total
end
