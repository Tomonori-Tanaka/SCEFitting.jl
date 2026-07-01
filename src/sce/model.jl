"""
    BasisSpec(; nbody, pair_cutoff, lmax, isotropy = false)

An SCE basis **specification** — the knobs that define the basis, not an interaction
term itself: maximum body order, the 2-body cutoff (Å), per-species maximum `l`, and
whether to keep only the isotropic (`Lf = 0`) channel.

`pair_cutoff` may be `Inf`, meaning "every resolvable pair" — the whole
Wigner–Seitz cell of the (super)cell — under the default [`MinimumImage`](@ref)
selection (cf. Magesty's `-1` sentinel).

!!! note
    Renamed from `Interaction` (the old name wrongly suggested a fitted coupling
    term).
"""
struct BasisSpec
    nbody::Int
    pair_cutoff::Float64
    lmax::Vector{Int}
    isotropy::Bool

    # Inner constructor only: validation can't be bypassed via the field constructor.
    function BasisSpec(; nbody::Integer, pair_cutoff::Real,
                       lmax::AbstractVector{<:Integer}, isotropy::Bool = false)
        nbody >= 1 || throw(ArgumentError("nbody must be ≥ 1; got $nbody"))
        (pair_cutoff > 0 && !isnan(pair_cutoff)) ||
            throw(ArgumentError("pair_cutoff must be positive (or Inf for the full " *
                                "Wigner–Seitz cell); got $pair_cutoff"))
        isempty(lmax) && throw(ArgumentError("lmax must have one entry per species"))
        all(>=(0), lmax) ||
            throw(ArgumentError("lmax entries must be ≥ 0; got $(collect(lmax))"))
        return new(Int(nbody), Float64(pair_cutoff), collect(Int, lmax), isotropy)
    end
end

"""
    SCEBasis(crystal, spec; backend = NoSymmetry(), tol = 1e-5,
             images = MinimumImage())

Build the SCE basis for `crystal`: analyze symmetry, enumerate cluster orbits, and
construct the symmetry-adapted SALC basis. Pass `backend = SpglibBackend()` (with
`using Spglib`) for real space-group symmetry.

`images` selects how periodic images are admitted (see [`AbstractImageSelection`](@ref)):
[`MinimumImage`](@ref) (the default) keeps only the minimum-image, resolvable pairs of
a plain-PBC supercell; [`AllImages`](@ref) keeps every image within the cutoff (the
generalized-Bloch / spin-spiral seam, finite cutoff only). The `images` value is also
applied to the cluster-edge admissibility, so the neighbor list and the clusters stay
consistent.
"""
struct SCEBasis
    crystal::Crystal
    spacegroup::SpaceGroup
    salc_basis::SALCBasis
    spec::BasisSpec
end

function SCEBasis(crystal::Crystal, spec::BasisSpec;
                 backend::AbstractSymmetryBackend = NoSymmetry(), tol::Real = 1e-5,
                 images::AbstractImageSelection = MinimumImage())::SCEBasis
    sg = analyze_symmetry(backend, crystal; tol = tol)
    nl = build_neighbor_list(crystal, spec.pair_cutoff, images)
    clusters = build_clusters(crystal, nl, sg; nbody = spec.nbody, selection = images)
    salcs = build_salc_basis(crystal, sg, clusters;
                             lmax_by_species = spec.lmax, isotropy = spec.isotropy)
    return SCEBasis(crystal, sg, salcs, spec)
end

"""
    n_salcs(basis::SCEBasis) -> Int

The number of SALC basis functions in `basis` — equivalently, the number of
design-matrix columns / fitted coefficients.
"""
n_salcs(b::SCEBasis) = length(b.salc_basis)

"""
    salcs(basis::SCEBasis) -> Vector{SALC}

The ordered SALC basis functions of `basis` (column order of the design matrix).
"""
salcs(b::SCEBasis) = b.salc_basis.salcs

"""
    SCEDataset(basis, configs, energies)
    SCEDataset(basis, configs, energies, torques)

Pair an [`SCEBasis`](@ref) with training data: spin configurations `configs`
(each `3 × n_atoms`, unit columns) and their `energies`. Materializes the energy
design matrix `X_E[config, salc] = evaluate_salc(salc, config)`.

The four-argument form additionally takes per-configuration torques (each
`3 × n_atoms`, the DFT torque `τ_a = m_a × B_a = −e_a × ∂E/∂e_a` on every atom) and builds the
torque design matrix `X_T` for an energy+torque co-fit (see [`fit`](@ref)). Its
rows are flattened config-major, then atom-major, then `xyz`:
`row = 3·n_atoms·(config−1) + 3·(atom−1) + component`. Torque-free datasets leave
`X_T`/`y_T` empty.
"""
struct SCEDataset
    basis::SCEBasis
    configs::Vector{Matrix{Float64}}
    X_E::Matrix{Float64}
    y_E::Vector{Float64}
    X_T::Matrix{Float64}
    y_T::Vector{Float64}
end

"""
    has_torque(dataset) -> Bool

Whether `dataset` carries torque training data (so an energy+torque co-fit is
possible).
"""
has_torque(d::SCEDataset) = !isempty(d.y_T)

# Spin configurations must be `3 × n_atoms` with finite, unit-norm columns — the
# contract the harmonic kernels assume (they call `Zlm_unsafe`, skipping per-call
# checks). Enforce it once at the data boundary so malformed DFT input fails loudly
# here instead of silently biasing the design matrix / prediction. `label` carries the
# config index into the message; `atol` is the unit-norm tolerance.
function _validate_config(c::AbstractMatrix{<:Real}, nat::Int; atol::Real = 1e-6,
                          label::AbstractString = "spin config")
    size(c, 1) == 3 ||
        throw(ArgumentError("$label must have 3 rows (got $(size(c, 1)))"))
    size(c, 2) == nat ||
        throw(DimensionMismatch("$label has $(size(c, 2)) atoms, basis expects $nat"))
    @inbounds for a in axes(c, 2)
        u = SVector{3,Float64}(c[1, a], c[2, a], c[3, a])
        all(isfinite, u) ||
            throw(ArgumentError("$label column $a is not finite ($(Tuple(u)))"))
        abs(norm(u) - 1) <= atol ||
            throw(ArgumentError("$label column $a is not a unit vector (‖e‖ = $(norm(u)))"))
    end
    return nothing
end

function _validate_configs(basis::SCEBasis, cfgs::Vector{Matrix{Float64}}; atol::Real = 1e-6)
    nat = n_atoms(basis.crystal)
    for (i, c) in enumerate(cfgs)
        _validate_config(c, nat; atol = atol, label = "config $i")
    end
    return nothing
end

function SCEDataset(basis::SCEBasis, configs::AbstractVector, energies::AbstractVector;
                   atol::Real = 1e-6)::SCEDataset
    length(configs) == length(energies) ||
        throw(DimensionMismatch("got $(length(configs)) configs but $(length(energies)) energies"))
    cfgs = [Matrix{Float64}(c) for c in configs]
    _validate_configs(basis, cfgs; atol = atol)
    X = _design_energy(basis, cfgs)
    empty_T = Matrix{Float64}(undef, 0, size(X, 2))
    return SCEDataset(basis, cfgs, X, collect(Float64, energies), empty_T, Float64[])
end

function SCEDataset(basis::SCEBasis, configs::AbstractVector, energies::AbstractVector,
                   torques::AbstractVector; atol::Real = 1e-6)::SCEDataset
    length(configs) == length(energies) ||
        throw(DimensionMismatch("got $(length(configs)) configs but $(length(energies)) energies"))
    cfgs = [Matrix{Float64}(c) for c in configs]
    length(torques) == length(cfgs) ||
        throw(ArgumentError("got $(length(torques)) torque blocks for $(length(cfgs)) configs"))
    _validate_configs(basis, cfgs; atol = atol)
    X_E = _design_energy(basis, cfgs)
    X_T = _design_torque(basis, cfgs)
    y_T = _flatten_torques(torques, cfgs)
    return SCEDataset(basis, cfgs, X_E, collect(Float64, energies), X_T, y_T)
end
"""
    SCEPredictor

The **lightweight, persistable predictor**: the [`SCEBasis`](@ref) plus the fitted
reference energy `j0` and SALC coefficients `jphi` (one per basis function), with
`keys` recording each coefficient's [`SALCKey`](@ref). Unlike [`SCEFit`](@ref) it
carries no training data, design matrices, or diagnostics — just enough to evaluate
[`predict_energy`](@ref) / [`predict_torque`](@ref) and to round-trip through
`SCEFitting.save` / `SCEFitting.load`.

Obtain one from a fit with `SCEPredictor(fit)`. Within a session `jphi[k]` pairs with
`salcs(basis)[k]` positionally; on reload the coefficients are re-paired to a freshly
built basis **by `SALCKey`** (`keys`), so a persisted model stays valid across sessions.
"""
struct SCEPredictor
    basis::SCEBasis
    j0::Float64
    jphi::Vector{Float64}
    keys::Vector{SALCKey}
end

"""
    SCEPredictor(basis, j0, jphi) -> SCEPredictor

Assemble a predictor directly from a basis and coefficients — `jphi[k]` pairs with
`salcs(basis)[k]` positionally and the `keys` are filled in from the basis. The public
way to build a synthetic model (hand-set couplings for tests, demos, or model studies)
without touching the SALC-basis internals; a fitted model comes from `SCEPredictor(fit)`.
"""
function SCEPredictor(basis::SCEBasis, j0::Real, jphi::AbstractVector{<:Real})::SCEPredictor
    length(jphi) == n_salcs(basis) ||
        throw(DimensionMismatch("jphi has $(length(jphi)) coefficients for " *
                                "$(n_salcs(basis)) SALC basis functions"))
    return SCEPredictor(basis, Float64(j0), collect(Float64, jphi), basis.salc_basis.keys)
end

"""
    SCEFit

The **full result of [`fit`](@ref)**: the [`SCEDataset`](@ref) (with its design
matrices), the fitted `j0`/`jphi`, the `estimator`, the `torque_weight` used, and the
(energy) residuals. This is the heavyweight, data-bearing object you query for
diagnostics ([`r2_energy`](@ref), [`rmse_energy`](@ref), [`residuals_energy`](@ref), …).
For prediction and storage, convert it to the lightweight [`SCEPredictor`](@ref) with
`SCEPredictor(fit)`.
"""
struct SCEFit
    dataset::SCEDataset
    j0::Float64
    jphi::Vector{Float64}
    estimator::AbstractEstimator
    residuals::Vector{Float64}
    torque_weight::Float64
end
