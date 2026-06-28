module MagestyRebuildSunnyExt

# Assembles a Sunny.jl `System` from a fitted SCEModel. All the conversion math
# (folded tesseral tensor → Cartesian matrix, directed-member folding, the energy
# gate) lives in the core (`src/sce/sunny.jl`); this extension only places the
# core-computed matrices into Sunny constructs, so it stays thin and the numerics
# are validated without Sunny.

using MagestyRebuild: SCEModel, SunnyTerms, _sunny_supercell_terms, num_atoms
import MagestyRebuild: to_sunny
using Sunny
using StaticArrays
using LinearAlgebra: dot

# Per-atom effective spin S: a scalar (uniform) or a `species_label => S` mapping.
function _resolve_spins(model::SCEModel, spins)::Vector{Float64}
    cr = model.basis.crystal
    nat = num_atoms(cr)
    if spins isa Real
        return fill(Float64(spins), nat)
    elseif spins isa AbstractDict
        return [Float64(spins[cr.species_labels[cr.species[i]]]) for i = 1:nat]
    end
    error("`spins` must be a number or a `species_label => S` Dict")
end

# Classical (`:dipole_uncorrected`) vs quantum (`:dipole`) single-ion rescaling. The
# bilinear exchange is mode-independent (`Sₐ·J·S_b`); only the rank-2 single-ion
# term is renormalized in the corrected dipole mode.
function _onsite_factor(s::Float64, mode::Symbol)::Float64
    mode === :dipole_uncorrected && return 1 / s^2
    mode === :dipole && return 2 / (s * (2s - 1))
    error("unsupported mode $mode for single-ion anisotropy (use :dipole or :dipole_uncorrected)")
end

# Build the Sunny anisotropy operator `factor · Σ_ij A[i,j] (Sᵢ Sⱼ + Sⱼ Sᵢ)/2` from
# the traceless symmetric Cartesian matrix `A` and the site spin operators `S`.
_onsite_operator(A::SMatrix{3,3,Float64,9}, S, factor::Float64) =
    factor * sum(A[i, j] * (S[i] * S[j] + S[j] * S[i]) / 2 for i = 1:3, j = 1:3)

function _build_supercell(model::SCEModel, spinvec::Vector{Float64}, g::Real,
                          mode::Symbol)
    cr = model.basis.crystal
    nat = num_atoms(cr)
    latvecs = Matrix{Float64}(cr.lattice.vectors)
    positions = [Vector{Float64}(cr.frac_positions[:, i]) for i = 1:nat]
    types = [cr.species_labels[cr.species[i]] for i = 1:nat]
    # symbol = 1 forces P1: every supercell atom is its own sublattice.
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

function to_sunny(model::SCEModel; spins, g::Real = 2, mode::Symbol = :dipole,
                  placement::Symbol = :auto)
    spinvec = _resolve_spins(model, spins)
    # The supercell (explicit) route is always exact. The primitive unfold is added
    # on top; for now route everything through the supercell.
    sys, skipped = _build_supercell(model, spinvec, g, mode)
    isempty(skipped) ||
        @warn "to_sunny: $(length(skipped)) SALC channel(s) are not Sunny-representable " *
              "and were skipped:\n  " * join(skipped, "\n  ")
    return sys
end

end # module MagestyRebuildSunnyExt
