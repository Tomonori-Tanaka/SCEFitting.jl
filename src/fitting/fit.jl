# Fitting and prediction: assemble the regression problem, run the estimator
# (fit/refit), and evaluate the fitted surface (predict_energy/predict_torque).

# Assemble the regression problem the estimator actually sees: the column-centered
# energy-only design (`w == 0`), or the whitened, stacked energy+torque design
# (`w > 0`), together with the centering constants `xbar`/`ybar` (for the analytic `j0`)
# and the per-row resampling-unit `groups`. Shared by `fit` and `refit` so the two build
# byte-identical designs.
function _assemble_problem(dataset::SCEDataset, w::Float64)
    X_E = dataset.X_E
    y_E = dataset.y_E
    xbar = vec(mean(X_E; dims = 1))
    ybar = mean(y_E)
    if w > 0
        n_E = length(y_E)
        n_T = length(dataset.y_T)
        se = sqrt((1 - w) / n_E)
        sm = sqrt(w / n_T)
        X = vcat((X_E .- xbar') .* se, dataset.X_T .* sm)
        y = vcat((y_E .- ybar) .* se, dataset.y_T .* sm)
        # Resampling-unit labels for grouped cross-validation: a configuration's energy
        # row and all its (config-major) torque-component rows share its index, so a
        # CV-based estimator never splits one configuration across folds.
        block = div(n_T, n_E)                       # = 3·n_atoms torque rows per config
        groups = vcat(collect(1:n_E), repeat(1:n_E; inner = block))
    else
        # THE SAME 1/√n_E, so the objective is CONTINUOUS at w = 0. This branch used
        # to return the centered design unscaled, i.e. it minimized the energy SSE while
        # every other weight setting minimizes the energy MSE — which is what `fit`'s
        # docstring says the objective is. Nothing changes for OLS (scale-invariant), but
        # for every penalized estimator the Gram jumped by `n_E` across the boundary, so
        # `lambda` silently meant something `n_E` times different: measured on 60
        # configurations, `Ridge(lambda = 1.0)` gave `rmse_energy` 0.0038 at `w = 0` and
        # 0.086 at `w = 1e-12` — an infinitesimal weight was not an infinitesimal change.
        # `AdaptiveRidge`/`GroupAdaptiveRidge`'s `epsilon` floor and every λ a selection
        # path records read the same scale, so they move with it.
        # [Backported from SLCE.jl 572cbe0. BREAKING: λ for an energy-only penalized
        # fit is n_E times smaller than before.]
        se = sqrt(1 / length(y_E))
        X = (X_E .- xbar') .* se
        y = (y_E .- ybar) .* se
        groups = nothing                            # each row is its own configuration
    end
    return X, y, xbar, ybar, groups
end

"""
    fit(SCEFit, dataset, estimator; torque_weight = 0.0) -> SCEFit

Fit the SCE coefficients. The energy design matrix is column-centered, so the
reference energy `j0` is recovered analytically as `mean(y_E − X_E·jϕ)`
(independent of the estimator), and the centered problem is handed to
`solve_coefficients`.

With `torque_weight = w ∈ (0, 1]` and a torque-carrying `dataset` (see
[`SCEDataset`](@ref)), the fit minimizes the per-sample-normalized objective

    L = (1 − w)·MSE_energy + w·MSE_torque

by whitening: the centered energy block is row-scaled by `√((1−w)/n_E)` and the
torque block by `√(w/n_T)`, then stacked. `j0` does not enter the torque block, so
it is still recovered from the energy data alone. `w = 0` (the default) is a
pure-energy fit; requesting `w > 0` without torque data is an error.

A resampling estimator (a cross-validating [`ElasticNet`](@ref) / [`Lasso`](@ref))
receives per-row group labels so its folds are **grouped by configuration** — a
configuration's energy row and its torque-component rows are never split across the
train/holdout boundary, which would otherwise leak within-configuration structure
into the CV estimate and bias λ selection.

**Scale of the two terms.** `MSE_energy` is the mean squared error of the **total**
energy of a configuration, `MSE_torque` that of one torque component of one site.
On a cell of `N` atoms the energy residual carries the site residuals summed `N`
times, so at equal per-site accuracy `MSE_energy ≈ N²·MSE_torque`-ish and a
weight that *reads* torque-dominated can still be energy-dominated: on a
1296-site cell `w = 0.99` leaves `0.01·MSE_energy` above `0.99·MSE_torque` and
the torque error of the fit only starts to move near `w = 0.999–0.9999`. Use
[`torque_weight_per_site`](@ref) to state the weight on the per-site energy
scale instead, or scan `w` on a held-out set ([`cross_validate`](@ref)).
"""
function fit(::Type{SCEFit}, dataset::SCEDataset, estimator::AbstractEstimator;
             torque_weight::Real = 0.0)::SCEFit
    isempty(dataset.y_E) && throw(ArgumentError("dataset has no observations"))
    w = Float64(torque_weight)
    (0.0 <= w <= 1.0) || throw(ArgumentError("torque_weight must be in [0, 1]; got $w"))
    if w > 0 && !has_torque(dataset)
        throw(ArgumentError("torque_weight = $w but the dataset has no torque data; " *
                            "build it with SCEDataset(basis, configs, energies, torques)"))
    end
    _check_metric_provenance(estimator, :energy, dataset.basis.salc_basis.fingerprint, w)
    X, y, xbar, ybar, groups = _assemble_problem(dataset, w)
    jphi = solve_coefficients(estimator, X, y; groups = groups)
    j0 = ybar - dot(xbar, jphi)
    residuals = dataset.y_E .- (j0 .+ dataset.X_E * jphi)
    return SCEFit(dataset, j0, jphi, estimator, residuals, w, nothing)
end

"""
    torque_weight_per_site(w_site, n_atoms) -> Float64

The `torque_weight` that makes [`fit`](@ref)'s objective

    L = (1 − w)·MSE_energy + w·MSE_torque

equal, up to an overall factor, to the **per-site** objective

    L_site = (1 − w_site)·MSE_energy / n_atoms² + w_site·MSE_torque,

i.e. the weight the user means when the energy is thought of per site
(`E / n_atoms`) rather than per cell. Matching the two ratios gives

    w = w_site·n_atoms² / (w_site·n_atoms² + 1 − w_site),

with `w_site = 0 → 0` and `w_site = 1 → 1`. Every estimator that is invariant
to an overall scale of the design (OLS, Ridge at a rescaled `lambda`) returns
the same coefficients for `fit(…; torque_weight = torque_weight_per_site(w_site,
n))` as for the per-site-scaled problem; a penalized estimator's `lambda` is on
the scale of the objective `fit` actually assembles (see [`ElasticNet`](@ref)).
"""
function torque_weight_per_site(w_site::Real, n_atoms::Integer)::Float64
    (0 <= w_site <= 1) || throw(ArgumentError("w_site must be in [0, 1]; got $w_site"))
    n_atoms >= 1 || throw(ArgumentError("n_atoms must be ≥ 1; got $n_atoms"))
    ws = Float64(w_site)
    n2 = Float64(n_atoms)^2
    return ws * n2 / (ws * n2 + 1 - ws)
end

"""
    refit(f::SCEFit, estimator = OLS(); threshold = 0.0) -> SCEFit

Re-fit on the **support** of an existing fit `f` — the classic de-biasing step that
follows a sparse fit: select the basis functions a regularized estimator kept, then
re-solve on just those columns (by default with [`OLS`](@ref)) to remove the shrinkage
bias on the survivors. The dataset and `torque_weight` of `f` are reused, so the design is
assembled exactly as in [`fit`](@ref).

A column `j` is in the support when its **scaled-magnitude contribution** exceeds
`threshold`:

    |coef(f)[j]| · ‖X[:, j]‖ > threshold

where `X` is the assembled (centered / whitened) design. `threshold = 0` (the default)
keeps the columns of nonzero scaled magnitude — the nonzero support of `f`, minus any
column the centering annihilates to a zero column. Coefficients off the support are set to
zero; `j0` is recovered analytically as in `fit`. If the support is empty, an all-zero `jϕ`
is returned (with a warning) and `j0` falls back to `mean(y_E)`.

A [`PrecomputedPilot`](@ref) (or an [`AdaptiveLasso`](@ref) whose pilot is one) is
rejected: its fixed coefficient vector has the original column count, not the refit
support length, and is meaningless once a support has been chosen.
"""
function refit(f::SCEFit, estimator::AbstractEstimator = OLS();
               threshold::Real = 0.0)::SCEFit
    threshold >= 0 || throw(ArgumentError("refit threshold must be ≥ 0; got $threshold"))
    _reject_precomputed_pilot(estimator)
    # Same shape of up-front refusal: the group labels are indexed by the FULL basis
    # design, and a support has been chosen — letting it through dies deeper with a
    # DimensionMismatch blaming labels built on a different basis, which is not
    # what happened. [Backported from SLCE.jl 54457ca, review M3.]
    estimator isa GroupAdaptiveRidge && throw(ArgumentError(
        "refit does not accept a GroupAdaptiveRidge: its column_groups are indexed " *
        "by the full basis design, not the chosen support. De-bias with OLS() (the " *
        "default) or Ridge — the group structure already did its job selecting."))
    dataset = f.dataset
    w = f.torque_weight
    # The same door `fit` / `select_fit` / `cross_validate` carry: a metric built for
    # another channel, another basis, or another torque weight is invisible to every
    # numerical gate, and a support reduction only ever checks its LENGTH.
    _check_metric_provenance(estimator, :energy, dataset.basis.salc_basis.fingerprint, w)
    X, y, xbar, ybar, groups = _assemble_problem(dataset, w)
    jphi_in = f.jphi
    # Scaled-magnitude support on the assembled design: a column survives when its
    # contribution magnitude |jϕ_j|·‖X[:, j]‖ clears the threshold.
    support = findall(j -> abs(jphi_in[j]) * norm(@view X[:, j]) > threshold,
                      eachindex(jphi_in))
    jphi = zeros(Float64, length(jphi_in))
    if isempty(support)
        # Short-circuit before `solve_coefficients`: a GLMNet-backed estimator would error
        # on a zero-column design, and the all-zero `jϕ` is a well-defined (degenerate) fit
        # whose `j0` falls back to `mean(y_E)`.
        @warn "refit: empty support — every coefficient is at or below the " *
              "scaled-magnitude threshold. Returning an all-zero jϕ; " *
              "j0 falls back to mean(y_E)." threshold
    else
        # A penalty metric is per-column data, so it has to be cut down to the chosen
        # support with the design — otherwise the solve dies on a length check whose
        # message blames a basis mismatch that did not happen. `_reduce_to_active` is
        # a no-op for an estimator without one.
        active = falses(length(jphi_in))
        active[support] .= true
        jphi[support] .= solve_coefficients(_reduce_to_active(estimator, active),
                                            view(X, :, support), y; groups = groups)
    end
    j0 = ybar - dot(xbar, jphi)
    residuals = dataset.y_E .- (j0 .+ dataset.X_E * jphi)
    return SCEFit(dataset, j0, jphi, estimator, residuals, w, sort(support))
end

# A `PrecomputedPilot` carries a fixed, full-design coefficient vector. `refit` and
# `cross_validate` must both reject it (with rationales of their own), so the predicate
# lives here as the single definition both rejection sites share.
_carries_precomputed_pilot(e::AbstractEstimator)::Bool =
    e isa PrecomputedPilot || (e isa AdaptiveLasso && e.pilot isa PrecomputedPilot)

# `refit` re-solves on a column sub-matrix `X[:, support]`, so a `PrecomputedPilot` —
# whose fixed coefficient vector carries the original full column count — would throw a
# `DimensionMismatch` deep in `solve_coefficients`, and is meaningless once a support has
# been chosen. Reject it (and an `AdaptiveLasso` carrying one) upfront with a clear message.
function _reject_precomputed_pilot(e::AbstractEstimator)
    if e isa PrecomputedPilot
        throw(ArgumentError("refit does not accept a PrecomputedPilot: its fixed " *
            "coefficient vector has the original column count, not the refit support " *
            "length. Pass a fresh estimator such as OLS()."))
    elseif e isa AdaptiveLasso && e.pilot isa PrecomputedPilot
        throw(ArgumentError("refit does not accept an AdaptiveLasso whose pilot is a " *
            "PrecomputedPilot: the pilot's fixed coefficient vector has the original " *
            "column count, not the refit support length. Use a fresh pilot such as OLS()."))
    end
    return nothing
end

"""
    SCEPredictor(f::SCEFit) -> SCEPredictor

Extract the lightweight predictor from a fit. The fit's coefficients live on the
design columns ([`n_columns`](@ref)); here they are expanded to one coefficient per
SALC — the SALCs of a tied cross-orbit alias group each receive the tied column's
coefficient times their per-bond weight (`±1` for transported copies), and are
marked `split = :convention`.
"""
function SCEPredictor(f::SCEFit)::SCEPredictor
    basis = f.dataset.basis
    return SCEPredictor(basis, f.j0, _expand_coefficients(f.jphi, _column_ties(basis)),
                        basis.salc_basis.keys, _split_status(_column_ties(basis)))
end

# Read-out doors (`multipole_terms`, `bilinear_terms` / `to_sunny`) call this so a
# convention-split coefficient never leaves the package as a silent measurement.
function _warn_convention_split(model::SCEPredictor, what::AbstractString)
    n = count(j -> model.split[j] !== :free && model.jphi[j] != 0.0, eachindex(model.jphi))
    n == 0 && return nothing
    @warn "$what: $n coefficient(s) belong to tied cross-orbit alias groups and carry " *
          "a split between periodic images that the training cell did not determine " *
          "(`split = :convention`; `:legacy` for a pre-v7 file). Training-cell " *
          "energies, torques and every q = 0 quantity are unaffected; bond-resolved " *
          "J(R) tables, J(q) off Γ, and tiled-supercell energies of configurations " *
          "that are not periodic in the training cell depend on it. See " *
          "`alias_groups(model.basis)` and the `split` column of `coeftable(model)`." maxlog = 1 _id = Symbol("convention_split_", model.basis.salc_basis.fingerprint)
    return nothing
end

"""
    predict_energy(model, data) -> Float64 or Vector{Float64}

Predict the energy of a spin configuration (`3 × n_atoms` matrix) or a vector of
configurations. The vector form is evaluated in parallel over `Threads.nthreads()`
threads (set `julia -t` / `JULIA_NUM_THREADS`); the result is thread-count-independent.
"""
function predict_energy(model::SCEPredictor, config::AbstractMatrix{<:Real})::Float64
    _validate_config(config, n_atoms(model.basis.crystal))
    salcs = model.basis.salc_basis.salcs
    e = model.j0
    scratch = SALCScratch()                 # reused workspace (dnPl + harmonic tables)
    @inbounds for k in eachindex(model.jphi)
        e += model.jphi[k] * evaluate_salc(salcs[k], config, scratch)
    end
    return e
end
function predict_energy(model::SCEPredictor, configs::AbstractVector)::Vector{Float64}
    nat = n_atoms(model.basis.crystal)
    for (i, c) in enumerate(configs)                   # serial: clean errors, not wrapped
        _validate_config(c, nat; label = "config $i")
    end
    out = Vector{Float64}(undef, length(configs))
    Threads.@threads for i in eachindex(configs, out)  # independent slots
        out[i] = predict_energy(model, configs[i])
    end
    return out
end
predict_energy(f::SCEFit, data) = predict_energy(SCEPredictor(f), data)

"""
    predict_torque(model, config) -> Matrix{Float64}

Predict the per-atom torque `τ_a = −e_a × ∂E/∂e_a` of a spin configuration
(`3 × n_atoms`), returned as a `3 × n_atoms` matrix. A vector of configurations
returns a vector of such matrices. This is the Landau–Lifshitz / physical torque
`m_a × B_eff,a` (the negative rotation-gradient of the energy), the analytic
derivative of the same surface [`predict_energy`](@ref) evaluates, so the two are
consistent by construction (`τ = −e × ∇E`). The vector form is evaluated in parallel
over `Threads.nthreads()` threads; the result is thread-count-independent.
"""
function predict_torque(model::SCEPredictor, config::AbstractMatrix{<:Real})::Matrix{Float64}
    nat = n_atoms(model.basis.crystal)
    _validate_config(config, nat)
    salcs = model.basis.salc_basis.salcs
    G = zeros(Float64, 3, nat)
    scratch = SALCScratch()                 # reused workspace (dnPl + harmonic tables)
    @inbounds for k in eachindex(model.jphi)
        accumulate_grad!(G, salcs[k], config, model.jphi[k], scratch)
    end
    T = Matrix{Float64}(undef, 3, nat)
    @inbounds for a = 1:nat
        ea = SVector{3,Float64}(config[1, a], config[2, a], config[3, a])
        ga = SVector{3,Float64}(G[1, a], G[2, a], G[3, a])
        t = cross(ga, ea)                # τ = ∇E × e = −e × ∇E  (physical / LL torque)
        T[1, a] = t[1]
        T[2, a] = t[2]
        T[3, a] = t[3]
    end
    return T
end
function predict_torque(model::SCEPredictor, configs::AbstractVector)::Vector{Matrix{Float64}}
    nat = n_atoms(model.basis.crystal)
    for (i, c) in enumerate(configs)                   # serial: clean errors, not wrapped
        _validate_config(c, nat; label = "config $i")
    end
    out = Vector{Matrix{Float64}}(undef, length(configs))
    Threads.@threads for i in eachindex(configs, out)  # independent slots
        out[i] = predict_torque(model, configs[i])
    end
    return out
end
predict_torque(f::SCEFit, data) = predict_torque(SCEPredictor(f), data)
