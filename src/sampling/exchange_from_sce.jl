# Mean-field sampler — extracting an `ExchangeModel` / `MultipoleField` from a fitted SCE
# (see `docs/specs/mfa-sampling.md`). Separated from `sampling/exchange.jl` because it
# depends on the Sunny bilinear/single-ion extraction `_sunny_supercell_terms` (loaded after
# the SCE model and the Sunny conversion core).

# Build the directed bilinear tensor `bilinear[a,b] = S_ab` and single-ion `onsite[a] = A_a`
# from the Sunny extraction (no warnings; callers decide what to report). Returns
# `(bilinear, onsite, nselfbond)`; `nselfbond` counts on-site (a == image-of-a) bilinear
# terms, which the rigid-axis mean field does not represent (only reachable via AllImages).
function _extract_bilinear_onsite(model::SCEModel)
    terms = _sunny_supercell_terms(model)
    n = num_atoms(model.basis.crystal)
    bilinear = fill(zero(SMatrix{3,3,Float64,9}), n, n)
    nselfbond = 0
    for ((a, b, _), M) in terms.pairs
        if a == b
            nselfbond += 1
            continue
        end
        bilinear[a, b] += M
        bilinear[b, a] += SMatrix{3,3,Float64}(transpose(M))
    end
    onsite = fill(zero(SMatrix{3,3,Float64,9}), n)
    for (a, A) in terms.onsites
        onsite[a] += A
    end
    return bilinear, onsite, nselfbond
end

"""
    ExchangeModel(model::SCEModel) -> ExchangeModel

Extract the full bilinear exchange (`ls=[1,1]`: Heisenberg + DMI + anisotropic) and the
single-ion anisotropy (`ls=[2]`) of a fitted [`SCEModel`](@ref) into an
[`ExchangeModel`](@ref), reusing the Sunny per-bond `3×3` matrix and per-atom single-ion
extraction. The bond matrices are placed directionally (`bilinear[a,b] = S_ab`, the reverse
member transposed), so the molecular field is `g_a = Σ_b S_ab ⟨e_b⟩`. Only the
higher-order / higher-`l` SALCs (3-body and up) are dropped — captured instead by the full
[`MultipoleField`](@ref) path — and reported via `@warn`.
"""
function ExchangeModel(model::SCEModel)
    bilinear, onsite, nselfbond = _extract_bilinear_onsite(model)
    nskipped = length(_sunny_supercell_terms(model).skipped)
    nskipped > 0 && @warn "ExchangeModel keeps the bilinear (Heisenberg + DMI + anisotropic) " *
        "and single-ion channels; dropped $nskipped higher-order / higher-l SALC(s) " *
        "(use the full SCE `MFASampler(model; reference)` to keep them)."
    nselfbond > 0 && @warn "ExchangeModel: skipping $nselfbond on-site (a == image-of-a) " *
        "bilinear term(s); the rigid-axis mean field does not represent them (only reachable " *
        "via AllImages). The default MinimumImage selection drops such self-pairs."
    return ExchangeModel(bilinear; onsite = onsite)
end

"""
    MultipoleField(model::SCEModel) -> MultipoleField

Digest a fitted [`SCEModel`](@ref) into the full-multipole mean field (P4): one mean-field
term per cluster member / `l`-ordering (carrying `jϕ·(4π)^(N/2)`, the member atoms, the
per-site `ls`, and the folded coefficient tensor), the model `lmax`, and the bilinear
[`ExchangeModel`](@ref) used only for the `l=1` temperature scale. This keeps **all**
channels — bilinear, single-ion, and higher-order / many-body — so the mean field iterates
the full multipole averages `⟨Z_lm⟩`.
"""
function MultipoleField(model::SCEModel)
    n = num_atoms(model.basis.crystal)
    salcs = model.basis.salcs.salcs
    terms = _MFATerm[]
    lmax = 0
    @inbounds for k in eachindex(model.jphi)
        j = model.jphi[k]
        j == 0.0 && continue
        salc = salcs[k]
        scale = j * (4π)^(salc.body / 2)
        for mem in salc.members
            # The mean-field contraction decouples every site in a cluster (⟨∏Z⟩ → ∏⟨Z⟩),
            # which assumes distinct sites; the cluster enumeration drops reused-atom
            # clusters, so this holds, but assert it (a repeated site would need CG
            # recoupling, not a self-⟨Z⟩ factor).
            allunique(mem.atoms) || throw(ArgumentError(
                "MultipoleField: cluster member with a repeated atom $(mem.atoms); the " *
                "mean-field factorization assumes distinct sites"))
            for t in mem.terms
                push!(terms, _MFATerm(scale, mem.atoms, t.ls, t.folded))
                lmax = max(lmax, maximum(t.ls))
            end
        end
    end
    isempty(terms) && throw(ArgumentError(
        "the model has no spin-dependent SALCs with a nonzero coefficient"))
    bilinear, onsite, _ = _extract_bilinear_onsite(model)
    return MultipoleField(n, lmax, terms, ExchangeModel(bilinear; onsite = onsite))
end

"""
    MFASampler(model::SCEModel; reference) -> MFASampler

The full-multipole mean-field sampler (P4): build a [`MultipoleField`](@ref) from the
fitted SCE (keeping every channel — bilinear, single-ion, and higher-order / many-body) and
sample about `reference`. The `l=1` temperature scale `T_MF = ρ/3` comes from the bilinear
part; the single-site distribution (a Bingham / higher-multipole shape) is drawn with the
Metropolis engine. See [`MFASampler(::ExchangeModel)`](@ref) for the bilinear-only path.
"""
MFASampler(model::SCEModel; reference::AbstractMatrix{<:Real}) =
    MFASampler(MultipoleField(model); reference = reference)
