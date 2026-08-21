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

The adiabatic site-moment channel's trio, each `nothing` when absent (the
energy / torque paths never read them, and a datum without them is exactly what it
was before the channel existed):

- `moments_bare` — `3 × n_atoms` **bare** (non-smoothed) magnetic moment vectors
  `M_a` (μ_B; VASP `M_int`), the projection target of the adiabatic site-moment
  channel `y_a = ê_a · M_a`. Kept as a raw vector — components may be negative and
  the magnitude may pass through zero (only finiteness is validated), because the
  signed readout is exactly what makes the target analytic where `‖M‖ → 0`.
  Distinct from the smoothed decomposition `magmoms[a] * directions[:, a]` (VASP
  `MW_int`):
  the constraining penalty acts on `MW`, so `MW` supplies the configuration
  coordinates `e` and the torque channel, while `M_int` supplies the moment
  target — their ratio is configuration-dependent, so neither substitutes for the
  other and both are stored.
- `constraint_axes` — `3 × n_atoms` constraint axes `ê_a^c` (unit columns; an
  exactly-zero column is this package's convention for "no axis for this atom" —
  the generator emits zeros for unconstrained atoms). For a transverse-penalty
  constraint (`constraint_mode = 1`) this is the axis the constrained DFT run
  evaluated the adiabatic map along (VASP `M_CONSTR`), and it is **required** —
  the matrix must be present and carry at least one axis; an atom whose column is
  zero under mode 1 contributes **no** moment row and must never fall back to
  `directions` (the gates skip it, a consumer must too). For a direction-pinning
  constraint (`constraint_mode = 4`) it is optional but recommended — it feeds
  the axis-consistency gates at the dataset boundary. A mode-1 datum without
  `moments_bare` is legal (it carries the axis without the target) and is
  skipped by the gates.
- `constraint_mode` — which class of constrained-DFT scheme produced this datum:
  `1` = transverse-penalty type (the axis is prescribed, the sign of the moment
  along it is free) or `4` = direction-pinning type (the full direction is
  pinned). The numbers follow VASP's `I_CONSTRAINED_M`, but the key is the
  physical class — another code's scheme maps onto one of the two. `nothing`
  means "no constraint information": the datum can never feed the moment
  channel's evaluation-axis rule (mode 4 reads `ê` from `directions`, mode 1
  from `constraint_axes` — keyed here, deliberately never by which fields happen
  to be present).

Build it from raw per-atom moment vectors and the constraining field with
`SpinDatum(energy, moments, field)` (which derives directions, magnitudes, and
torques; the trio passes through as keywords), or construct the fields directly
— the five-argument form leaves the trio absent.
"""
struct SpinDatum <: AbstractTrainingDatum
    energy::Float64
    directions::Matrix{Float64}
    magmoms::Vector{Float64}
    field::Matrix{Float64}
    torques::Matrix{Float64}
    moments_bare::Union{Matrix{Float64},Nothing}
    constraint_axes::Union{Matrix{Float64},Nothing}
    constraint_mode::Union{Int,Nothing}

    function SpinDatum(energy::Float64, directions::Matrix{Float64},
                       magmoms::Vector{Float64}, field::Matrix{Float64},
                       torques::Matrix{Float64},
                       moments_bare::Union{Matrix{Float64},Nothing},
                       constraint_axes::Union{Matrix{Float64},Nothing},
                       constraint_mode::Union{Int,Nothing})
        nat = size(directions, 2)
        if moments_bare !== nothing
            size(moments_bare) == (3, nat) || throw(ArgumentError(
                "`moments_bare` must be 3 × $nat (got $(size(moments_bare)))"))
            # finiteness is deliberately the ONLY value constraint: the bare moment
            # is a signed vector whose magnitude legitimately passes through zero —
            # that analyticity is the whole point of the projection target y = ê·M
            all(isfinite, moments_bare) ||
                throw(ArgumentError("`moments_bare` contains non-finite entries"))
        end
        if constraint_axes !== nothing
            size(constraint_axes) == (3, nat) || throw(ArgumentError(
                "`constraint_axes` must be 3 × $nat (got $(size(constraint_axes)))"))
            @inbounds for a = 1:nat
                u = SVector{3,Float64}(constraint_axes[1, a], constraint_axes[2, a],
                                       constraint_axes[3, a])
                # An exactly-zero column asserts "no axis for this atom" (the VASP
                # M_CONSTR convention); anything else must be a unit axis. Near-zero
                # noise is neither — it is refused rather than silently normalized.
                u == SVector{3,Float64}(0, 0, 0) && continue
                all(isfinite, u) || throw(ArgumentError(
                    "`constraint_axes` column $a is not finite ($(Tuple(u)))"))
                abs(norm(u) - 1) <= _DIRECTION_ATOL || throw(ArgumentError(
                    "`constraint_axes` column $a is not a unit vector or exactly " *
                    "zero (‖u‖ = $(norm(u)))"))
                # The same component bound every direction door of this package
                # carries (`_validate_config`, `Harmonics._validate_unit`): an axis
                # is a direction the moment channel will hand to the harmonic
                # kernels, whose `dnPl` domain is |z| ≤ 1 — a near-pole column 5e-9
                # off unit clears any norm band and would throw a bare DomainError
                # from inside an accumulation. Not redundant with the band.
                maximum(abs, u) <= 1 || throw(ArgumentError(
                    "`constraint_axes` column $a has a component of magnitude " *
                    "$(maximum(abs, u)) > 1 (a unit axis cannot; the harmonic " *
                    "kernels' domain is |z| ≤ 1)"))
            end
        end
        if constraint_mode !== nothing
            constraint_mode in (1, 4) ||
                throw(ArgumentError("`constraint_mode` must be 1 (transverse-penalty " *
                                    "type) or 4 (direction-pinning type); got " *
                                    "$constraint_mode"))
            constraint_mode == 1 && constraint_axes === nothing &&
                throw(ArgumentError("`constraint_mode = 1` (transverse-penalty type) " *
                                    "requires `constraint_axes`: the moment readout " *
                                    "axis cannot be reconstructed from the converged " *
                                    "moment direction where ‖M‖ → 0, so the " *
                                    "constraint axis must be carried explicitly"))
            constraint_mode == 1 && all(iszero, constraint_axes) &&
                throw(ArgumentError("`constraint_mode = 1` with an all-zero " *
                                    "`constraint_axes`: every column says \"no " *
                                    "axis\", so the matrix carries none — a " *
                                    "transverse-penalty run constrains at least " *
                                    "one atom"))
        elseif constraint_axes !== nothing
            throw(ArgumentError("`constraint_axes` without `constraint_mode`: the " *
                                "evaluation-axis rule is keyed by the constraint " *
                                "class (1 or 4), deliberately never by which fields " *
                                "happen to be present — declare the mode"))
        end
        return new(energy, directions, magmoms, field, torques, moments_bare,
                   constraint_axes, constraint_mode)
    end
end

# The five-field form: the adiabatic-moment trio absent. Kept so every existing
# direct construction (and the v4-era docs) reads unchanged — converting, as the
# default constructor it replaces was.
SpinDatum(energy::Real, directions::AbstractMatrix{<:Real},
          magmoms::AbstractVector{<:Real}, field::AbstractMatrix{<:Real},
          torques::AbstractMatrix{<:Real}) =
    SpinDatum(Float64(energy), Matrix{Float64}(directions), Vector{Float64}(magmoms),
              Matrix{Float64}(field), Matrix{Float64}(torques), nothing, nothing,
              nothing)

"""
    SpinDatum(energy, moments, field; zero_moment_atol = 1e-10,
              moments_bare = nothing, constraint_axes = nothing,
              constraint_mode = nothing) -> SpinDatum

Build a [`SpinDatum`](@ref) from the per-atom magnetic moment vectors `moments`
(`3 × n_atoms`, μ_B) and the per-atom constraining field `field` (`3 × n_atoms`,
eV/μ_B). The spin direction is `e_a = m_a / ‖m_a‖` (a near-zero moment, below
`zero_moment_atol`, gets the placeholder `ẑ` and a zero torque), the magnitude is
`‖m_a‖`, and the torque target is `τ_a = m_a × B_a` (eV) — the physical /
Landau–Lifshitz torque, matching the SCE model torque `−e_a × ∂E/∂e_a`.

`moments_bare` / `constraint_axes` / `constraint_mode` pass through to
[`SpinDatum`](@ref) unchanged (see there): `moments` here is the smoothed
decomposition source (VASP `MW_int` — the quantity the constraint acts on), while
`moments_bare` is the bare `M_int` the adiabatic moment channel targets.

The target carries the per-config moment magnitude `‖m_a‖`, while the SCE model torque
depends on directions only — so a co-fit assumes the moment magnitudes
are roughly constant across configurations (large longitudinal variation would bias it).
The zero-moment placeholder direction is a divide-by-zero guard; a *magnetic* site that
quenches to `‖m_a‖ ≈ 0` in some configuration therefore enters with a fictitious
direction and a zero torque (a small bias) — prefer dropping such configurations.
(`SCEDataset` rejects a placeholder on a basis-referenced atom; if you change this
tolerance, pass the same value to its `zero_moment_atol` so the guard stays aligned.)
"""
function SpinDatum(energy::Real, moments::AbstractMatrix{<:Real},
                   field::AbstractMatrix{<:Real}; zero_moment_atol::Real = 1e-10,
                   moments_bare::Union{AbstractMatrix{<:Real},Nothing} = nothing,
                   constraint_axes::Union{AbstractMatrix{<:Real},Nothing} = nothing,
                   constraint_mode::Union{Integer,Nothing} = nothing)::SpinDatum
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
    _m(x) = x === nothing ? nothing : Matrix{Float64}(x)
    constraint_mode isa Bool &&
        throw(ArgumentError("`constraint_mode` is a class (1 or 4), not a flag"))
    return SpinDatum(Float64(energy), dirs, mags, Matrix{Float64}(field), torq,
                     _m(moments_bare), _m(constraint_axes),
                     constraint_mode === nothing ? nothing : Int(constraint_mode))
end

"""
    read_configs(src::AbstractDFTSource) -> Vector{SpinDatum}

Read all training configurations from a DFT source. Implemented per source type
(e.g. `SCETools.VASP.Oszicar` in the SCETools.jl package).
"""
read_configs(src::AbstractDFTSource) =
    throw(ArgumentError("read_configs is not implemented for $(typeof(src))"))

# Every atom the SALC basis references must carry a nonzero magnetic moment in every
# configuration: a quenched (‖m‖ ≈ 0) moment on a referenced atom enters the design
# matrix through the ẑ placeholder direction of `SpinDatum` and silently biases the
# fit. Unreferenced atoms (species removed with `lmax = 0`, or sites outside every
# admitted cluster) may be non-magnetic — their moments are never consulted.
function _check_referenced_moments(basis::SCEBasis, data::AbstractVector{SpinDatum};
                                   atol::Real)
    ref = _referenced_atoms(basis)
    nat = length(ref)
    labels = basis.crystal.species_labels
    for (i, d) in enumerate(data)
        length(d.magmoms) == nat ||
            throw(DimensionMismatch("config $i has $(length(d.magmoms)) atoms, " *
                                    "basis expects $nat"))
        for a = 1:nat
            (ref[a] && d.magmoms[a] <= atol) || continue
            lab = labels[basis.crystal.species[a]]
            throw(ArgumentError(
                "config $i: atom $a ($lab) has a zero magnetic moment " *
                "(‖m‖ = $(d.magmoms[a]) ≤ $atol) but is referenced by the SALC basis " *
                "— its placeholder ẑ direction would silently bias the fit. Drop the " *
                "configuration, or, if the species is non-magnetic, remove it from " *
                "the basis with lmax = 0"))
        end
    end
    return nothing
end

"""
    SCEDataset(basis, data::AbstractVector{SpinDatum}; use_torque = true,
               zero_moment_atol = 1e-10) -> SCEDataset
    SCEDataset(basis, src::AbstractDFTSource; use_torque = true,
               zero_moment_atol = 1e-10) -> SCEDataset

Build a fit-ready [`SCEDataset`](@ref) from training data (or directly from a DFT
source, which is read first). The spin directions become the configurations and the
energies the energy targets; with `use_torque = true` the per-atom torque targets
are included as well (for an energy+torque co-fit — see [`fit`](@ref)). Pass
`use_torque = false` for an energy-only dataset (e.g. unconstrained data with no
meaningful field).

Every atom the SALC basis references must carry a nonzero moment (`> zero_moment_atol`,
μ_B) in every configuration — a quenched moment would otherwise enter the fit through
the `ẑ` placeholder direction of [`SpinDatum`](@ref) and silently bias it, so it is an
error. Atoms the basis never reads (a species removed with `lmax = 0`, or sites outside
every admitted cluster) are exempt. The guard re-derives quenched atoms from the stored
magnitudes, so if the `SpinDatum`s were built with a custom `zero_moment_atol`, pass the
same value here — a looser build tolerance with the default guard tolerance would let a
placeholder direction through (both default to `1e-10`).
"""
function SCEDataset(basis::SCEBasis, data::AbstractVector{SpinDatum};
                    use_torque::Bool = true,
                    zero_moment_atol::Real = 1e-10)::SCEDataset
    isempty(data) && throw(ArgumentError("no training data"))
    _check_referenced_moments(basis, data; atol = zero_moment_atol)
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

# `zero_moment_atol` is forwarded: the docstring tells a custom-atol adapter to pass
# the same value here, and this path used to silently drop it.
# [Backported from SLCE.jl 011e3c9 (the spin-applicable hunk of the co-fit commit).]
SCEDataset(basis::SCEBasis, src::AbstractDFTSource; use_torque::Bool = true,
           zero_moment_atol::Real = 1e-10)::SCEDataset =
    SCEDataset(basis, read_configs(src); use_torque = use_torque,
               zero_moment_atol = zero_moment_atol)
