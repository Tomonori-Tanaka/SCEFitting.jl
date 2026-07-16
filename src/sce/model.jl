"""
    BasisSpec([labels_or_crystal]; nbody, lmax, cutoff, lsum = nothing,
              isotropy = false)

An SCE basis **specification** — the knobs that define the basis, not an
interaction term itself. Stored in resolved, dense canonical form:

- `nbody::Int` — maximum body order.
- `lmax::Vector{Int}` — per-species maximum `l` (per site, every body order;
  `0` removes the species from the basis entirely).
- `lsum::Vector{Int}` — per-body-order cap on `Σl` over the cluster sites
  (index = body order; `LSUM_UNCAPPED` where uncapped).
- `cutoff::Vector{Matrix{Float64}}` — per-body-order symmetric species-pair
  cutoff radii in Å (`cutoff[N - 1][a, b]` for body order `N ≥ 2`); an
  `N`-body cluster is kept iff **every** internal edge is within its own
  species-pair radius. `Inf` = no cutoff (every resolvable pair — the whole
  Wigner–Seitz cell under [`MinimumImage`](@ref)); `0` excludes the pair.
- `isotropy::Bool` — keep only the scalar (`Lf = 0`) channel.
- `species_labels::Vector{String}` — the labels the spec was resolved against
  (empty when constructed index-keyed).

The keyword constructor also accepts ergonomic forms and resolves them here:

```julia
BasisSpec(crystal; nbody = 3, isotropy = true,          # or BasisSpec(labels; ...)
    lmax   = ["*" => 3, "B" => 0],                       # label-keyed, "*" fallback
    lsum   = [1 => 0, 2 => 4, 3 => 4],                   # body-keyed (or one Int)
    cutoff = [2 => Inf, 3 => ["Fe-*" => 6.0, "*-*" => 8.0]])
```

`lmax` may be an `Int` (all species), a per-species `Vector{Int}`, or
label-keyed pairs; `lsum` an `Int` (all body orders), body-keyed pairs, or
`nothing` (no cap); `cutoff` a scalar, a pair table (all body orders),
body-keyed pairs of either, or the canonical `Vector{Matrix}`. Pair keys are
unordered (`"Fe-Nd" ≡ "Nd-Fe"`) and resolve by specificity (concrete beats
`"A-*"` beats `"*-*"`); equal-specificity conflicts, unknown labels, uncovered
species/pairs, and body orders outside `nbody` are errors. Label-keyed forms
need the labels: pass them (or a `Crystal`) as the first argument.
"""
struct BasisSpec
    nbody::Int
    lmax::Vector{Int}
    lsum::Vector{Int}
    cutoff::Vector{Matrix{Float64}}
    isotropy::Bool
    species_labels::Vector{String}

    # Positional canonical-form constructor; all construction paths validate here.
    function BasisSpec(nbody::Int, lmax::Vector{Int}, lsum::Vector{Int},
                       cutoff::Vector{Matrix{Float64}}, isotropy::Bool,
                       species_labels::Vector{String})
        nbody >= 1 || throw(ArgumentError("nbody must be ≥ 1; got $nbody"))
        isempty(lmax) && throw(ArgumentError("lmax must have one entry per species"))
        all(>=(0), lmax) ||
            throw(ArgumentError("lmax entries must be ≥ 0; got $lmax"))
        length(lsum) == nbody ||
            throw(ArgumentError("lsum has $(length(lsum)) entries for nbody = $nbody"))
        all(>=(0), lsum) || throw(ArgumentError("lsum entries must be ≥ 0; got $lsum"))
        length(cutoff) == nbody - 1 ||
            throw(ArgumentError("cutoff has $(length(cutoff)) matrices for body " *
                                "orders 2:$nbody"))
        nkd = length(lmax)
        for (k, M) in enumerate(cutoff)
            size(M) == (nkd, nkd) ||
                throw(ArgumentError("cutoff body$(k + 1) matrix is $(size(M)), " *
                                    "expected ($nkd, $nkd)"))
            M == M' ||
                throw(ArgumentError("cutoff body$(k + 1) matrix is not symmetric"))
            all(v -> !isnan(v) && v >= 0, M) ||
                throw(ArgumentError("cutoff body$(k + 1) entries must be ≥ 0 or Inf"))
        end
        isempty(species_labels) || length(species_labels) == nkd ||
            throw(ArgumentError("$(length(species_labels)) species labels for " *
                                "$nkd lmax entries"))
        allunique(species_labels) ||
            throw(ArgumentError("species labels must be unique; got $species_labels"))
        return new(nbody, lmax, lsum, cutoff, isotropy, species_labels)
    end
end

function BasisSpec(labels::AbstractVector{<:AbstractString} = String[];
                   nbody::Integer, lmax, cutoff, lsum = nothing,
                   isotropy::Bool = false, pair_cutoff = nothing)
    pair_cutoff === nothing ||
        throw(ArgumentError("`pair_cutoff` was replaced by `cutoff` (a scalar is " *
                            "equivalent: cutoff = $pair_cutoff; see the BasisSpec " *
                            "docstring for per-body / per-pair forms)"))
    nbody >= 1 || throw(ArgumentError("nbody must be ≥ 1; got $nbody"))
    lv = collect(String, labels)
    nkd = !isempty(lv) ? length(lv) :
          lmax isa AbstractVector{<:Integer} ? length(lmax) :
          throw(ArgumentError("the species count is unknown: pass the labels " *
                              "(`BasisSpec(labels; ...)` / `BasisSpec(crystal; ...)`) " *
                              "or give `lmax` as a per-species Vector{Int}"))
    return BasisSpec(Int(nbody),
                     _resolve_species_table(lmax, nkd, lv, "lmax"),
                     _resolve_lsum(lsum, Int(nbody)),
                     _resolve_cutoff(cutoff, Int(nbody), nkd, lv),
                     isotropy, lv)
end

BasisSpec(crystal::Crystal; kwargs...) = BasisSpec(crystal.species_labels; kwargs...)

Base.:(==)(a::BasisSpec, b::BasisSpec) =
    a.nbody == b.nbody && a.lmax == b.lmax && a.lsum == b.lsum &&
    a.cutoff == b.cutoff && a.isotropy == b.isotropy &&
    a.species_labels == b.species_labels

function Base.show(io::IO, ::MIME"text/plain", sp::BasisSpec)
    nkd = length(sp.lmax)
    labels = isempty(sp.species_labels) ? ["#$k" for k = 1:nkd] : sp.species_labels
    println(io, "BasisSpec: nbody = ", sp.nbody, ", isotropy = ", sp.isotropy)
    println(io, "  lmax:  ", join(("$(labels[k]) = $(sp.lmax[k])" for k = 1:nkd), ", "))
    lsumstr(v) = v == LSUM_UNCAPPED ? "—" : string(v)
    println(io, "  lsum:  ",
            join(("body$N = $(lsumstr(sp.lsum[N]))" for N = 1:sp.nbody), ", "))
    for N = 2:sp.nbody
        M = sp.cutoff[N - 1]
        if all(==(M[1, 1]), M)
            println(io, "  cutoff body$N: ", M[1, 1], " Å (all pairs)")
        else
            println(io, "  cutoff body$N (Å):")
            for i = 1:nkd
                println(io, "    ", rpad(labels[i], 6), " ",
                        join((string(M[i, j]) for j = 1:nkd), "  "))
            end
        end
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
    length(spec.lmax) == length(crystal.species_labels) ||
        throw(ArgumentError("spec covers $(length(spec.lmax)) species, crystal has " *
                            "$(length(crystal.species_labels))"))
    isempty(spec.species_labels) || spec.species_labels == crystal.species_labels ||
        throw(ArgumentError("spec species labels $(spec.species_labels) do not match " *
                            "the crystal's $(crystal.species_labels)"))
    sg = analyze_symmetry(backend, crystal; tol = tol)
    # The neighbor list is built at the per-pair superset radius (element-wise max
    # over body orders); each body order's own radii then trim edges per cluster in
    # `candidate_clusters`.
    nl = build_neighbor_list(crystal, _superset_cutoff(spec), images)
    clusters = build_clusters(crystal, nl, sg; nbody = spec.nbody, selection = images,
                              cutoff = spec.cutoff)
    salcs = build_salc_basis(crystal, sg, clusters;
                             lmax_by_species = spec.lmax, lsum_by_body = spec.lsum,
                             isotropy = spec.isotropy)
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
