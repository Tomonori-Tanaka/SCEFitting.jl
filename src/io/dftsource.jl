"""
    AbstractDFTSource

A source of DFT training data: `read_configs(src::AbstractDFTSource) ->
Vector{<:AbstractTrainingDatum}` turns some DFT output into fit-ready configurations,
and `SCEDataset(basis, src)` goes straight from a source to a dataset. Concrete sources
(e.g. `SCETools.VASP.Oszicar`, in the SCETools.jl package) live alongside their format
reader, so the SCE pipeline consumes only [`SpinDatum`](@ref) / [`SCEDataset`](@ref) and
never depends on the originating DFT code.
"""
abstract type AbstractDFTSource end

"""
    AbstractTrainingDatum

One DFT configuration's observables, as consumed by the SCE pipeline. The concrete type
is [`SpinDatum`](@ref) (energy, per-atom spin directions, magnetic moments, constraining
field, and the derived torque target).
"""
abstract type AbstractTrainingDatum end

"""
    SpinDatum

One DFT spin configuration used for fitting:

- `energy::Float64` — total energy (eV).
- `directions::Matrix{Float64}` — `3 × n_atoms` unit spin directions `e_a`.
- `magmoms::Vector{Float64}` — per-atom moment magnitude `m_a` (μ_B).
- `field::Matrix{Float64}` — `3 × n_atoms` constraining field `B_a` (eV/μ_B), from a
  constrained-noncollinear calculation (zero where unconstrained).
- `torques::Matrix{Float64}` — `3 × n_atoms` per-atom torque target
  `τ_a = m_a × B_a = ‖m_a‖ (e_a × B_a)` (eV), the physical / Landau–Lifshitz
  torque, the observable that the SCE torque `τ_a = −e_a × ∂E/∂e_a` is fit to.

Build it from raw per-atom moment vectors and the constraining field with
`SpinDatum(energy, moments, field)` (which derives directions, magnitudes, and
torques), or construct the fields directly.
"""
struct SpinDatum <: AbstractTrainingDatum
    energy::Float64
    directions::Matrix{Float64}
    magmoms::Vector{Float64}
    field::Matrix{Float64}
    torques::Matrix{Float64}
end

"""
    SpinDatum(energy, moments, field; zero_moment_atol = 1e-10) -> SpinDatum

Build a [`SpinDatum`](@ref) from the per-atom magnetic moment vectors `moments`
(`3 × n_atoms`, μ_B) and the per-atom constraining field `field` (`3 × n_atoms`,
eV/μ_B). The spin direction is `e_a = m_a / ‖m_a‖` (a near-zero moment, below
`zero_moment_atol`, gets the placeholder `ẑ` and a zero torque), the magnitude is
`‖m_a‖`, and the torque target is `τ_a = m_a × B_a` (eV) — the physical /
Landau–Lifshitz torque, matching the SCE model torque `−e_a × ∂E/∂e_a`.

The target carries the per-config moment magnitude `‖m_a‖`, while the SCE model torque
depends on directions only — so a co-fit assumes the moment magnitudes
are roughly constant across configurations (large longitudinal variation would bias it).
The zero-moment placeholder direction is a divide-by-zero guard; a *magnetic* site that
quenches to `‖m_a‖ ≈ 0` in some configuration therefore enters with a fictitious
direction and a zero torque (a small bias) — prefer dropping such configurations.
"""
function SpinDatum(energy::Real, moments::AbstractMatrix{<:Real},
                   field::AbstractMatrix{<:Real}; zero_moment_atol::Real = 1e-10)::SpinDatum
    size(moments, 1) == 3 || throw(ArgumentError("`moments` must be 3 × n_atoms"))
    size(field) == size(moments) ||
        throw(ArgumentError("`field` $(size(field)) must match `moments` $(size(moments))"))
    n = size(moments, 2)
    dirs = Matrix{Float64}(undef, 3, n)
    mags = Vector{Float64}(undef, n)
    torq = Matrix{Float64}(undef, 3, n)
    @inbounds for i = 1:n
        mi = SVector{3,Float64}(moments[1, i], moments[2, i], moments[3, i])
        Bi = SVector{3,Float64}(field[1, i], field[2, i], field[3, i])
        mag = norm(mi)
        mags[i] = mag
        ei = mag <= zero_moment_atol ? SVector{3,Float64}(0, 0, 1) : mi / mag
        ti = cross(mi, Bi)                  # τ = m × B  (physical / LL torque)  [eV]
        for k = 1:3
            dirs[k, i] = ei[k]
            torq[k, i] = ti[k]
        end
    end
    return SpinDatum(Float64(energy), dirs, mags, Matrix{Float64}(field), torq)
end

"""
    read_configs(src::AbstractDFTSource) -> Vector{SpinDatum}

Read all training configurations from a DFT source. Implemented per source type
(e.g. `SCETools.VASP.Oszicar` in the SCETools.jl package).
"""
read_configs(src::AbstractDFTSource) =
    throw(ArgumentError("read_configs is not implemented for $(typeof(src))"))

"""
    SCEDataset(basis, data::AbstractVector{SpinDatum}; use_torque = true) -> SCEDataset
    SCEDataset(basis, src::AbstractDFTSource; use_torque = true) -> SCEDataset

Build a fit-ready [`SCEDataset`](@ref) from training data (or directly from a DFT
source, which is read first). The spin directions become the configurations and the
energies the energy targets; with `use_torque = true` the per-atom torque targets
are included as well (for an energy+torque co-fit — see [`fit`](@ref)). Pass
`use_torque = false` for an energy-only dataset (e.g. unconstrained data with no
meaningful field).
"""
function SCEDataset(basis::SCEBasis, data::AbstractVector{SpinDatum};
                    use_torque::Bool = true)::SCEDataset
    isempty(data) && throw(ArgumentError("no training data"))
    configs = [d.directions for d in data]
    energies = Float64[d.energy for d in data]
    if use_torque
        torques = [d.torques for d in data]
        any(t -> any(!iszero, t), torques) ||
            throw(ArgumentError("use_torque = true but every torque target is zero — no " *
                                "constraining field was found in the data; pass " *
                                "use_torque = false for an energy-only dataset"))
        return SCEDataset(basis, configs, energies, torques)
    else
        return SCEDataset(basis, configs, energies)
    end
end

SCEDataset(basis::SCEBasis, src::AbstractDFTSource; use_torque::Bool = true)::SCEDataset =
    SCEDataset(basis, read_configs(src); use_torque = use_torque)
