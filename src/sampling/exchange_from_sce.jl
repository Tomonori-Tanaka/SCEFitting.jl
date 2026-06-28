# Mean-field sampler — extracting an `ExchangeModel` from a fitted SCE (see
# `docs/specs/mfa-sampling.md`). Separated from `sampling/exchange.jl` because it depends on
# the Sunny bilinear/single-ion extraction `_sunny_supercell_terms` (loaded after the SCE
# model and the Sunny conversion core).

"""
    ExchangeModel(model::SCEModel) -> ExchangeModel

Extract the full bilinear exchange (`ls=[1,1]`: Heisenberg + DMI + anisotropic) and the
single-ion anisotropy (`ls=[2]`) of a fitted [`SCEModel`](@ref) into an
[`ExchangeModel`](@ref), reusing the Sunny per-bond `3×3` matrix and per-atom single-ion
extraction. The bond matrices are placed directionally (`bilinear[a,b] = S_ab`, the reverse
member transposed), so the molecular field is `g_a = Σ_b S_ab ⟨e_b⟩`. Only the
higher-order / higher-`l` SALCs (3-body and up) are dropped — a P4 extension — and reported
via `@warn`.
"""
function ExchangeModel(model::SCEModel)
    terms = _sunny_supercell_terms(model)
    n = num_atoms(model.basis.crystal)
    bilinear = fill(zero(SMatrix{3,3,Float64,9}), n, n)
    selfbond = false
    for ((a, b, _), M) in terms.pairs
        if a == b
            selfbond = true                          # on-site (a–image-of-a) bilinear bond
            continue
        end
        bilinear[a, b] += M
        bilinear[b, a] += SMatrix{3,3,Float64}(transpose(M))
    end
    selfbond && @warn "ExchangeModel: skipping on-site (a == image-of-a) bilinear term(s); " *
        "the rigid-axis mean field does not represent them (only reachable via AllImages, a " *
        "future spin-spiral path). The default MinimumImage selection drops such self-pairs."
    onsite = fill(zero(SMatrix{3,3,Float64,9}), n)
    for (a, A) in terms.onsites
        onsite[a] += A
    end
    isempty(terms.skipped) || @warn "ExchangeModel keeps the bilinear (Heisenberg + DMI + " *
        "anisotropic) and single-ion channels; dropped $(length(terms.skipped)) higher-order " *
        "/ higher-l SALC(s) (a P4 extension)."
    return ExchangeModel(bilinear; onsite = onsite)
end
