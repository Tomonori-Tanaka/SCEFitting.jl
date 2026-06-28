"""
Sunny export — the code-agnostic core.

A fitted [`SCEModel`](@ref) is a polynomial in unit spin directions of arbitrary
body order and `l`. Sunny.jl represents only **bilinear pair exchange** (a 3×3
matrix per bond) and **single-ion anisotropy** (a quadratic form per site), so the
export covers exactly the two SCE channels that map onto those:

- `ls = [1, 1]` (2-body, the `Lf = 0,1,2` parts fold into one 3×3 matrix that
  carries Heisenberg, Dzyaloshinskii–Moriya, and symmetric-anisotropic exchange),
- `ls = [2]` (1-body, a traceless symmetric single-ion tensor).

`ls = [0…]` channels are constants (folded into `j0`, dropped); every other SALC
(3-body and up, higher `l`) is **skipped** and reported.

The numerically delicate part — turning a fitted coefficient and its folded
tesseral tensor into a Cartesian matrix with the right normalization, and summing
the directed cluster members into one matrix per undirected bond — lives here, in
the core, and is validated by reconstructing the SCE energy
(`_reconstruct_energy ≈ predict_energy − j0`) **without** Sunny. The actual
`Sunny.System` assembly is the thin `MagestyRebuildSunnyExt` extension (loaded by
`using Sunny`), which consumes these matrices.
"""

# Cartesian axis (x,y,z) = (1,2,3) → the m-index of the folded tensor, with index
# (1,2,3) ↔ m = (−1,0,+1): real `Z_{1,m}` is `√(3/4π)·(y,z,x)` for m = (−1,0,+1),
# so x↔m=+1 (index 3), y↔m=−1 (index 1), z↔m=0 (index 2).
const _L1_AXIS_TO_MIDX = (3, 1, 2)

"""
    _l1_pair_matrix(folded) -> SMatrix{3,3}

The Cartesian bilinear matrix `M` with `eₐ'·M·e_b = Σ_{m₁,m₂} folded[m₁,m₂]
Z_{1,m₁}(eₐ) Z_{1,m₂}(e_b)`, i.e. `(3/4π)` times the folded `(m₁,m₂)` tensor
reindexed to Cartesian axes. `folded` is the rank-2 term tensor of an `ls = [1,1]`
SALC (already `Lf`-contracted), so its symmetric part is Heisenberg + Γ and its
antisymmetric part is the DM vector.
"""
function _l1_pair_matrix(folded::AbstractArray)::SMatrix{3,3,Float64,9}
    pref = 3 / (4π)
    p = _L1_AXIS_TO_MIDX
    return SMatrix{3,3,Float64}(ntuple(9) do k
        r = (k - 1) % 3 + 1
        c = (k - 1) ÷ 3 + 1
        pref * folded[p[r], p[c]]
    end...)
end

"""
    _l2_onsite_matrix(folded) -> SMatrix{3,3}

The traceless symmetric single-ion matrix `A` with `e'·A·e = Σ_m folded[m]
Z_{2,m}(e)`, from the rank-1 term tensor of an `ls = [2]` SALC (`folded[1..5]` ↔
`m = −2..+2`). Constants `√(15/16π)`, `√(5/16π)` carry the `Z_{2,m}`
normalization.
"""
function _l2_onsite_matrix(folded::AbstractArray)::SMatrix{3,3,Float64,9}
    a = sqrt(15 / (16π))
    b = sqrt(5 / (16π))
    Bxy = a * folded[1]      # m = −2
    Byz = a * folded[2]      # m = −1
    d0 = folded[3]           # m =  0
    Bxz = a * folded[4]      # m = +1
    d2 = folded[5]           # m = +2
    Bxx = -b * d0 + a * d2
    Byy = -b * d0 - a * d2
    Bzz = 2b * d0
    return SMatrix{3,3,Float64}(Bxx, Bxy, Bxz,
                                Bxy, Byy, Byz,
                                Bxz, Byz, Bzz)
end

# Which Sunny channel an SALC maps to, from its sorted l-multiset label.
function _classify_salc(ls::AbstractVector{<:Integer})::Symbol
    all(==(0), ls) && return :drop          # spin-independent constant (in j0)
    ls == [1, 1] && return :pair            # bilinear exchange
    ls == [2] && return :onsite             # single-ion anisotropy
    return :unsupported                      # 3-body+, higher-l pairs
end

_unsupported_salc_string(k::SALCKey) =
    "SALC(body=$(k.body), ls=$(k.ls), Lf=$(k.Lf)): unsupported by Sunny export " *
    "(only ls=[1,1] pairs and ls=[2] single-ion are exported)"

"""
    SunnyTerms

The supercell-route intermediate of a Sunny export: one bilinear matrix per
directed supercell bond `(a, b, R)` (canonicalized to `a ≤ b`), one single-ion
matrix per atom, and the human-readable list of skipped SALCs. Built by
[`_sunny_supercell_terms`](@ref); consumed by the `Sunny` extension and by the
Sunny-free energy reconstruction [`_reconstruct_energy`](@ref).
"""
struct SunnyTerms
    pairs::Dict{Tuple{Int,Int,SVector{3,Int}},SMatrix{3,3,Float64,9}}
    onsites::Dict{Int,SMatrix{3,3,Float64,9}}
    skipped::Vector{String}
end

# Accumulate `m` into `d[key]` (zero-initialized).
_accum!(d, key, m::SMatrix{3,3,Float64,9}) =
    (d[key] = get(d, key, zero(SMatrix{3,3,Float64,9})) + m)

"""
    _sunny_supercell_terms(model) -> SunnyTerms

Walk the fitted SALCs and fold each Sunny-representable channel onto the training
supercell: an `ls=[1,1]` member `(a, b, R)` contributes `jϕ·(4π)^(N/2)·_l1_pair_matrix`
to the bond, summed over the two directed members `(a,b,R)`/`(b,a,−R)` into one
matrix on the canonical `a ≤ b` bond (the reverse member transposed); an `ls=[2]`
member contributes `jϕ·(4π)^(N/2)·_l2_onsite_matrix` to its atom. The energy is then
`Σ eₐ'·M·e_b + Σ eₐ'·A·eₐ`, matching `predict_energy − j0` for a model with only
these channels.
"""
function _sunny_supercell_terms(model::SCEModel)::SunnyTerms
    pairs = Dict{Tuple{Int,Int,SVector{3,Int}},SMatrix{3,3,Float64,9}}()
    onsites = Dict{Int,SMatrix{3,3,Float64,9}}()
    skipped = String[]
    salcs = model.basis.salcs.salcs
    @inbounds for k in eachindex(model.jphi)
        salc = salcs[k]
        j = model.jphi[k]
        cls = _classify_salc(salc.ls)
        cls === :drop && continue
        if cls === :unsupported
            push!(skipped, _unsupported_salc_string(salc.key))
            continue
        end
        scale = j * (4π)^(salc.body / 2)
        scale == 0.0 && continue
        for m in salc.members
            for t in m.terms
                if cls === :pair
                    a, b = m.atoms[1], m.atoms[2]
                    C = scale * _l1_pair_matrix(t.folded)
                    R = m.shifts[2] - m.shifts[1]
                    if a <= b
                        _accum!(pairs, (a, b, R), C)
                    else
                        _accum!(pairs, (b, a, -R), SMatrix{3,3,Float64}(transpose(C)))
                    end
                else # :onsite
                    _accum!(onsites, m.atoms[1], scale * _l2_onsite_matrix(t.folded))
                end
            end
        end
    end
    return SunnyTerms(pairs, onsites, skipped)
end

"""
    _reconstruct_energy(terms, e) -> Float64

The classical energy `Σ_{(a,b,R)} eₐ'·M·e_b + Σ_a eₐ'·A·eₐ` of a [`SunnyTerms`](@ref)
on a spin configuration `e` (`3 × n_atoms`, unit columns). For a model whose only
channels are `ls=[1,1]`/`ls=[2]` this equals `predict_energy(model, e) − j0` — the
Sunny-free gate on the conversion math.
"""
function _reconstruct_energy(terms::SunnyTerms, e::AbstractMatrix{<:Real})::Float64
    tot = 0.0
    @inbounds for ((a, b, _), M) in terms.pairs
        ea = SVector{3,Float64}(e[1, a], e[2, a], e[3, a])
        eb = SVector{3,Float64}(e[1, b], e[2, b], e[3, b])
        tot += dot(ea, M * eb)
    end
    @inbounds for (a, A) in terms.onsites
        ea = SVector{3,Float64}(e[1, a], e[2, a], e[3, a])
        tot += dot(ea, A * ea)
    end
    return tot
end

"""
    to_sunny(model; spins, g = 2, mode = :auto, placement = :auto) -> Sunny.System

Export a fitted [`SCEModel`](@ref) to a Sunny.jl `System`. **Requires `using Sunny`**
(the export is a package extension). Only the Sunny-representable channels are
exported — bilinear pair (`ls=[1,1]`) exchange and single-ion (`ls=[2]`) anisotropy;
higher-order / higher-`l` SALCs are skipped and reported via `@warn`.

`spins` gives the effective spin length `S` (a number, or a per-species
`label => S` mapping); the exchange matrices are rescaled by `1/(SₐS_b)` so the
Sunny system's classical energy equals `predict_energy − j0`. `placement = :primitive`
unfolds onto the chemical primitive cell (physical dispersion) when the model maps
cleanly, else falls back to the training supercell (`:explicit`); `:auto` picks
primitive when clean.
"""
function to_sunny(model::SCEModel, args...; kwargs...)
    error("to_sunny requires Sunny.jl; run `using Sunny` first " *
          "(it activates the MagestyRebuild Sunny export extension).")
end
