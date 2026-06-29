# Fitted-model introspection — a stable, code-neutral view of the multipole terms of a
# fitted SCE. Downstream packages (e.g. the mean-field samplers in `SCETools.jl`) consume
# the fitted Hamiltonian through this surface instead of reaching into the SALC-basis
# internals (`model.basis.salc_basis.salcs`, `SALCMember` / `SALCTerm` fields), so the basis
# representation can evolve without breaking them. The general per-term dump is
# `multipole_terms`; the bilinear / single-ion `3×3` extraction (reusing the validated
# Sunny conversion in `interop/sunny.jl`) is `bilinear_terms`.

"""
    MultipoleTerm

One multipole term of a fitted [`SCEPredictor`](@ref): the **raw** fitted coefficient
`coef = jϕ_k` (the per-N scale `(4π)^(body/2)` is *not* applied — apply it once, at the
consumer, from `body`), the cluster `body` order, the member `atoms` and their periodic
`shifts`, the per-site angular momenta `ls`, and the rank-`body` real coefficient tensor
`folded`. The term's energy contribution on a configuration `e` (`3 × n_atoms`, unit
columns) is

```
coef · (4π)^(body/2) · Σ_μ folded[μ] ∏ᵢ Z_{ls[i], μ_i}(e_{atoms[i]}),   μ_i = idx_i − ls[i] − 1
```

(the same kernel as [`evaluate_salc`](@ref) for a single SALC term). See [`multipole_terms`](@ref).
"""
struct MultipoleTerm
    coef::Float64
    body::Int
    atoms::Vector{Int}
    shifts::Vector{SVector{3,Int}}
    ls::Vector{Int}
    folded::Array{Float64}
end

# The number of atoms in the (training-cell) crystal the model was fitted on (a method of
# the exported `n_atoms`, so downstream consumers need not reach into `model.basis.crystal`).
n_atoms(model::SCEPredictor)::Int = n_atoms(model.basis.crystal)

"""
    multipole_terms(model::SCEPredictor) -> Vector{MultipoleTerm}

A code-neutral, flat view of a fitted SCE: one [`MultipoleTerm`](@ref) per cluster member
and `l`-ordering of every SALC with a nonzero coefficient. This is the stable public
contract a downstream consumer (a mean-field sampler, an energy evaluator, …) reads
*instead of* the SALC-basis internals; the per-N scale `(4π)^(body/2)` is left for the
consumer to apply once, so the returned `coef` is exactly the fitted `jϕ`.
"""
function multipole_terms(model::SCEPredictor)::Vector{MultipoleTerm}
    salcs = model.basis.salc_basis.salcs
    out = MultipoleTerm[]
    @inbounds for k in eachindex(model.jphi)
        j = model.jphi[k]
        j == 0.0 && continue
        salc = salcs[k]
        for mem in salc.members
            for t in mem.terms
                push!(out, MultipoleTerm(j, salc.body, mem.atoms, mem.shifts, t.ls, t.folded))
            end
        end
    end
    return out
end

"""
    bilinear_terms(model::SCEPredictor)

The bilinear pair (`ls=[1,1]`: Heisenberg + Dzyaloshinskii–Moriya + symmetric-anisotropic)
and single-ion (`ls=[2]`) channels of a fitted SCE, extracted as Cartesian `3×3` matrices
and folded onto the training supercell. Returns a named tuple `(; pairs, onsites, skipped)`:

- `pairs::Dict{Tuple{Int,Int,SVector{3,Int}}, SMatrix{3,3}}` — one matrix `M` per directed
  bond `(a, b, R)` (canonicalized to `a ≤ b`), energy `Σ eₐ'·M·e_b`;
- `onsites::Dict{Int, SMatrix{3,3}}` — one single-ion matrix `A` per atom, energy `eₐ'·A·eₐ`;
- `skipped::Vector{String}` — the higher-order / higher-`l` SALCs not representable as a
  bilinear pair or single-ion term.

Reuses the same validated tesseral→Cartesian conversion as the Sunny export
([`to_sunny`](@ref)); the higher-order channels it lists are kept by the full
[`multipole_terms`](@ref) view.
"""
function bilinear_terms(model::SCEPredictor)
    t = _bilinear_terms(model)
    return (; pairs = t.pairs, onsites = t.onsites, skipped = t.skipped)
end
