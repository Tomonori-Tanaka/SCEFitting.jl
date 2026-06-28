module MagestyRebuildSunnyExt

# Assembles a Sunny.jl `System` from a fitted SCEModel. All the conversion math
# (folded tesseral tensor → Cartesian matrix, directed-member folding, the supercell
# → primitive unfold, the energy gate) lives in the core (`src/sce/sunny.jl`); this
# extension only places the core-computed matrices into Sunny constructs, so it stays
# thin and the numerics are validated without Sunny.

using MagestyRebuild: SCEModel, SunnyTerms, SunnyPrimitive, _sunny_supercell_terms,
    _sunny_primitive, num_atoms
import MagestyRebuild: to_sunny
using Sunny
using StaticArrays
using LinearAlgebra: dot

# A `label -> S` effective-spin lookup from the user's `spins` (scalar or mapping).
function _spin_lookup(spins)
    _checked(s) = (s > 0 ? Float64(s) :
                   error("effective spin S must be positive; got $s"))
    spins isa Real && (s = _checked(spins); return _ -> s)
    if spins isa AbstractDict
        return label -> (haskey(spins, label) ? _checked(spins[label]) :
                         error("`spins` has no entry for species \"$label\""))
    end
    error("`spins` must be a number or a `species_label => S` Dict")
end

# Classical (`:dipole_uncorrected`) vs quantum (`:dipole`) single-ion rescaling. The
# bilinear exchange is mode-independent (`Sₐ·J·S_b`); only the rank-2 single-ion
# term is renormalized in the corrected dipole mode (where spin-1/2 has no rank-2
# observable, so the factor is undefined).
function _onsite_factor(s::Float64, mode::Symbol)::Float64
    mode === :dipole_uncorrected && return 1 / s^2
    if mode === :dipole
        s > 0.5 ||
            error("mode = :dipole has no rank-2 single-ion anisotropy for s ≤ 1/2 " *
                  "(got s = $s); use :dipole_uncorrected or drop the ls=[2] channel.")
        return 2 / (s * (2s - 1))
    end
    error("unsupported mode $mode for single-ion anisotropy (use :dipole or :dipole_uncorrected)")
end

# The Sunny anisotropy operator `factor · Σ_ij A[i,j] (Sᵢ Sⱼ + Sⱼ Sᵢ)/2` from the
# traceless symmetric Cartesian matrix `A` and the site spin operators `S`.
_onsite_operator(A::SMatrix{3,3,Float64,9}, S, factor::Float64) =
    factor * sum(A[i, j] * (S[i] * S[j] + S[j] * S[i]) / 2 for i = 1:3, j = 1:3)

# The training supercell as a P1 Sunny crystal (every atom its own sublattice).
function _build_supercell(model::SCEModel, lookup, g::Real, mode::Symbol)
    cr = model.basis.crystal
    nat = num_atoms(cr)
    latvecs = Matrix{Float64}(cr.lattice.vectors)
    positions = [Vector{Float64}(cr.frac_positions[:, i]) for i = 1:nat]
    types = [cr.species_labels[cr.species[i]] for i = 1:nat]
    spinvec = [lookup(types[i]) for i = 1:nat]
    cryst = Sunny.Crystal(latvecs, positions, 1; types = types)
    moments = [i => Sunny.Moment(s = spinvec[i], g = Float64(g)) for i = 1:nat]
    sys = Sunny.to_inhomogeneous(Sunny.System(cryst, moments, mode))

    terms = _sunny_supercell_terms(model)
    for ((a, b, R), M) in terms.pairs
        J = M / (spinvec[a] * spinvec[b])
        Sunny.set_exchange_at!(sys, Matrix(J), (1, 1, 1, a), (1, 1, 1, b);
                               offset = (R[1], R[2], R[3]))
    end
    for (a, A) in terms.onsites
        f = _onsite_factor(spinvec[a], mode)
        Sunny.set_onsite_coupling_at!(sys, S -> _onsite_operator(A, S, f), (1, 1, 1, a))
    end
    return sys, terms.skipped
end

# The chemical primitive cell, with one Sunny bond per primitive bond (Sunny's
# periodic replication restores the supercell multiplicity).
function _build_primitive(prim::SunnyPrimitive, lookup, g::Real, mode::Symbol)
    latvecs = Matrix{Float64}(prim.latvecs)
    positions = [Vector{Float64}(p) for p in prim.positions]
    types = prim.types
    spinvec = [lookup(types[i]) for i in eachindex(types)]
    # symbol = 1 forces P1: bonds are placed explicitly from the SCE model, so Sunny
    # must not re-derive (and symmetrize over) the crystal symmetry.
    cryst = Sunny.Crystal(latvecs, positions, 1; types = types)
    moments = [i => Sunny.Moment(s = spinvec[i], g = Float64(g)) for i in eachindex(types)]
    sys = Sunny.System(cryst, moments, mode)
    for ((i, j, n), M) in prim.bonds
        J = M / (spinvec[i] * spinvec[j])
        Sunny.set_exchange!(sys, Matrix(J), Sunny.Bond(i, j, [n[1], n[2], n[3]]))
    end
    for (i, A) in prim.onsites
        f = _onsite_factor(spinvec[i], mode)
        Sunny.set_onsite_coupling!(sys, S -> _onsite_operator(A, S, f), i)
    end
    return sys
end

function to_sunny(model::SCEModel; spins, g::Real = 2, mode::Symbol = :dipole,
                  placement::Symbol = :auto)
    placement in (:auto, :primitive, :explicit) ||
        error("placement must be :auto, :primitive, or :explicit; got $placement")
    lookup = _spin_lookup(spins)
    local sys, skipped
    if placement === :explicit
        sys, skipped = _build_supercell(model, lookup, g, mode)
    else
        prim = _sunny_primitive(model)
        if prim.clean
            sys = _build_primitive(prim, lookup, g, mode)
            skipped = prim.skipped
        else
            msg = "to_sunny: model does not unfold cleanly onto the primitive cell; " *
                  "using the exact supercell route instead (dispersion will be folded)."
            placement === :primitive ? @warn(msg) : @info(msg)
            sys, skipped = _build_supercell(model, lookup, g, mode)
        end
    end
    isempty(skipped) ||
        @warn "to_sunny: $(length(skipped)) SALC channel(s) are not Sunny-representable " *
              "and were skipped:\n  " * join(skipped, "\n  ")
    return sys
end

end # module MagestyRebuildSunnyExt
