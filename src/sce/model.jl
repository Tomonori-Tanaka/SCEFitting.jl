"""
    Interaction(; nbody, pair_cutoff, lmax, isotropy = false)

An SCE interaction **specification** — the knobs that define the basis, not an
interaction term itself: maximum body order, the 2-body cutoff (Å), per-species
maximum `l`, and whether to keep only the isotropic (`Lf = 0`) channel.

`pair_cutoff` may be `Inf`, meaning "every resolvable pair" — the whole
Wigner–Seitz cell of the (super)cell — under the default [`MinimumImage`](@ref)
selection (cf. Magesty's `-1` sentinel).
"""
struct Interaction
    nbody::Int
    pair_cutoff::Float64
    lmax::Vector{Int}
    isotropy::Bool
end
function Interaction(; nbody::Integer, pair_cutoff::Real, lmax::AbstractVector{<:Integer},
                     isotropy::Bool = false)
    nbody >= 1 || throw(ArgumentError("nbody must be ≥ 1; got $nbody"))
    (pair_cutoff > 0 && !isnan(pair_cutoff)) ||
        throw(ArgumentError("pair_cutoff must be positive (or Inf for the full " *
                            "Wigner–Seitz cell); got $pair_cutoff"))
    return Interaction(Int(nbody), Float64(pair_cutoff), collect(Int, lmax), isotropy)
end

"""
    SCEBasis(crystal, interaction; backend = NoSymmetry(), tol = 1e-5,
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
    interaction::Interaction
end

function SCEBasis(crystal::Crystal, interaction::Interaction;
                 backend::AbstractSymmetryBackend = NoSymmetry(), tol::Real = 1e-5,
                 images::AbstractImageSelection = MinimumImage())::SCEBasis
    sg = analyze_symmetry(backend, crystal; tol = tol)
    nl = build_neighbor_list(crystal, interaction.pair_cutoff, images)
    clusters = build_clusters(crystal, nl, sg; nbody = interaction.nbody, selection = images)
    salcs = build_salc_basis(crystal, sg, clusters;
                             lmax_by_species = interaction.lmax, isotropy = interaction.isotropy)
    return SCEBasis(crystal, sg, salcs, interaction)
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
`3 × n_atoms`, the DFT torque `τ_a = e_a × ∂E/∂e_a` on every atom) and builds the
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

function SCEDataset(basis::SCEBasis, configs::AbstractVector, energies::AbstractVector)::SCEDataset
    cfgs = [Matrix{Float64}(c) for c in configs]
    X = _design_energy(basis, cfgs)
    empty_T = Matrix{Float64}(undef, 0, size(X, 2))
    return SCEDataset(basis, cfgs, X, collect(Float64, energies), empty_T, Float64[])
end

function SCEDataset(basis::SCEBasis, configs::AbstractVector, energies::AbstractVector,
                   torques::AbstractVector)::SCEDataset
    cfgs = [Matrix{Float64}(c) for c in configs]
    length(torques) == length(cfgs) ||
        throw(ArgumentError("got $(length(torques)) torque blocks for $(length(cfgs)) configs"))
    X_E = _design_energy(basis, cfgs)
    X_T = _design_torque(basis, cfgs)
    y_T = _flatten_torques(torques, cfgs)
    return SCEDataset(basis, cfgs, X_E, collect(Float64, energies), X_T, y_T)
end

# Energy design matrix: X_E[config, salc] = Φ_salc(config).
function _design_energy(basis::SCEBasis, cfgs::Vector{Matrix{Float64}})::Matrix{Float64}
    salcs = basis.salc_basis.salcs
    n = length(cfgs)
    m = length(salcs)
    X = Matrix{Float64}(undef, n, m)
    # Columns are independent; each task owns whole columns, so the writes are
    # disjoint (no false sharing) and the result is identical at any thread count.
    Threads.@threads for j = 1:m
        @inbounds for i = 1:n        # column-major: stride-1 writes down each column
            X[i, j] = evaluate_salc(salcs[j], cfgs[i])
        end
    end
    return X
end

# Torque design matrix: for each SALC column, each config, each atom, the three
# components of τ_a = −e_a × ∂Φ/∂e_a (the physical / Landau–Lifshitz torque) stacked
# config-major / atom-major / xyz.
function _design_torque(basis::SCEBasis, cfgs::Vector{Matrix{Float64}})::Matrix{Float64}
    salcs = basis.salc_basis.salcs
    n = length(cfgs)
    m = length(salcs)
    nat = isempty(cfgs) ? 0 : size(cfgs[1], 2)
    block = 3 * nat
    X = Matrix{Float64}(undef, n * block, m)
    # Columns are independent; each task owns whole columns. `G` must be
    # task-local (a shared buffer would race across threads).
    Threads.@threads for j = 1:m     # column per SALC (column-major writes)
        G = Matrix{Float64}(undef, 3, nat)
        @inbounds for ci = 1:n
            c = cfgs[ci]
            fill!(G, 0.0)
            accumulate_grad!(G, salcs[j], c, 1.0)
            row_off = block * (ci - 1)
            for a = 1:nat
                ea = SVector{3,Float64}(c[1, a], c[2, a], c[3, a])
                ga = SVector{3,Float64}(G[1, a], G[2, a], G[3, a])
                t = cross(ga, ea)        # τ = ∇Φ × e = −e × ∇Φ  (physical / LL torque)
                rb = row_off + 3 * (a - 1)
                X[rb + 1, j] = t[1]
                X[rb + 2, j] = t[2]
                X[rb + 3, j] = t[3]
            end
        end
    end
    return X
end

# Flatten per-config torque targets in the same row order as `_design_torque`.
function _flatten_torques(torques::AbstractVector, cfgs::Vector{Matrix{Float64}})::Vector{Float64}
    nat = isempty(cfgs) ? 0 : size(cfgs[1], 2)
    block = 3 * nat
    y = Vector{Float64}(undef, length(cfgs) * block)
    @inbounds for ci in eachindex(torques)
        τ = torques[ci]
        size(τ) == (3, nat) ||
            throw(ArgumentError("torque block $ci must be 3 × $nat (got $(size(τ)))"))
        rb = block * (ci - 1)
        for a = 1:nat, d = 1:3
            y[rb + 3 * (a - 1) + d] = τ[d, a]
        end
    end
    return y
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
        X = X_E .- xbar'
        y = y_E .- ybar
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
    X, y, xbar, ybar, groups = _assemble_problem(dataset, w)
    jphi = solve_coefficients(estimator, X, y; groups = groups)
    j0 = ybar - dot(xbar, jphi)
    residuals = dataset.y_E .- (j0 .+ dataset.X_E * jphi)
    return SCEFit(dataset, j0, jphi, estimator, residuals, w)
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
    dataset = f.dataset
    w = f.torque_weight
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
        jphi[support] .= solve_coefficients(estimator, view(X, :, support), y;
                                            groups = groups)
    end
    j0 = ybar - dot(xbar, jphi)
    residuals = dataset.y_E .- (j0 .+ dataset.X_E * jphi)
    return SCEFit(dataset, j0, jphi, estimator, residuals, w)
end

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

Extract the lightweight predictor from a fit.
"""
SCEPredictor(f::SCEFit) = SCEPredictor(f.dataset.basis, f.j0, f.jphi, f.dataset.basis.salc_basis.keys)

"""
    predict_energy(model, data) -> Float64 or Vector{Float64}

Predict the energy of a spin configuration (`3 × n_atoms` matrix) or a vector of
configurations. The vector form is evaluated in parallel over `Threads.nthreads()`
threads (set `julia -t` / `JULIA_NUM_THREADS`); the result is thread-count-independent.
"""
function predict_energy(model::SCEPredictor, config::AbstractMatrix{<:Real})::Float64
    salcs = model.basis.salc_basis.salcs
    e = model.j0
    @inbounds for k in eachindex(model.jphi)
        e += model.jphi[k] * evaluate_salc(salcs[k], config)
    end
    return e
end
function predict_energy(model::SCEPredictor, configs::AbstractVector)::Vector{Float64}
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
    salcs = model.basis.salc_basis.salcs
    nat = size(config, 2)
    G = zeros(Float64, 3, nat)
    @inbounds for k in eachindex(model.jphi)
        accumulate_grad!(G, salcs[k], config, model.jphi[k])
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
    out = Vector{Matrix{Float64}}(undef, length(configs))
    Threads.@threads for i in eachindex(configs, out)  # independent slots
        out[i] = predict_torque(model, configs[i])
    end
    return out
end
predict_torque(f::SCEFit, data) = predict_torque(SCEPredictor(f), data)

"""
    coef(fit_or_model) -> Vector{Float64}

The fitted SALC coefficients `Jϕ`, one per design-matrix column (in [`SALCKey`](@ref)
order). The reference energy `j0` is separate; read it with [`intercept`](@ref).
"""
coef(f::SCEFit) = f.jphi
coef(m::SCEPredictor) = m.jphi

"""
    intercept(fit_or_model) -> Float64

The reference energy `j0` (the SCE intercept), recovered analytically in [`fit`](@ref).
"""
intercept(f::SCEFit) = f.j0
intercept(m::SCEPredictor) = m.j0

"""
    nobs(f::SCEFit) -> Int

The number of energy observations used in the fit.
"""
nobs(f::SCEFit) = length(f.dataset.y_E)

"""
    dof(f::SCEFit) -> Int

Degrees of freedom consumed by the fit: `length(coef(f)) + 1` — one per SCE coefficient
plus the intercept `j0`. This is the full parametric count, not an effective / selected
count, even for a sparse estimator.
"""
dof(f::SCEFit)::Int = length(f.jphi) + 1

"""
    residuals_energy(f::SCEFit) -> Vector{Float64}

The per-configuration energy residuals `y_E − (j0 + X_E·jϕ)` (eV), one per training
configuration — the same residuals the energy `R²` / RMSE are built from.
"""
residuals_energy(f::SCEFit)::Vector{Float64} = f.residuals

"""
    residuals_torque(f::SCEFit) -> Vector{Float64}

The torque residuals `y_T − X_T·jϕ` over the flattened (config-major, then atom-major,
then `xyz`) per-atom torque components (eV). Requires a torque-carrying dataset; the
torque prediction has no intercept (`j0` does not enter it).
"""
function residuals_torque(f::SCEFit)::Vector{Float64}
    has_torque(f.dataset) || throw(ArgumentError("dataset has no torque data"))
    return f.dataset.y_T .- f.dataset.X_T * f.jphi
end

"""
    rss_energy(f::SCEFit) -> Float64

Energy residual sum of squares `Σ (y_E − ŷ_E)²` (eV²).
"""
rss_energy(f::SCEFit)::Float64 = sum(abs2, residuals_energy(f))

"""
    rss_torque(f::SCEFit) -> Float64

Torque residual sum of squares `Σ (y_T − X_T·jϕ)²` (eV²). Requires a torque-carrying
dataset.
"""
rss_torque(f::SCEFit)::Float64 = sum(abs2, residuals_torque(f))

"""
    r2_energy(f::SCEFit) -> Float64

In-sample energy coefficient of determination `R²` (1 = the SALC span reproduces every
training energy). A degenerate target (all training energies equal ⇒ zero total variance)
returns `1.0` by convention.
"""
function r2_energy(f::SCEFit)::Float64
    y = f.dataset.y_E
    ss_tot = sum(abs2, y .- mean(y))
    return ss_tot == 0 ? 1.0 : 1 - rss_energy(f) / ss_tot
end

"""
    rmse_energy(f::SCEFit) -> Float64

In-sample energy root-mean-square error (same unit as the training energies, typically eV).
"""
rmse_energy(f::SCEFit)::Float64 = sqrt(rss_energy(f) / nobs(f))

"""
    r2_torque(f::SCEFit) -> Float64

In-sample torque `R²`, computed over the flattened per-atom torque components. Requires a
torque-carrying dataset. The torque prediction `τ = X_T·jϕ` has no intercept (`j0` does
not enter it), so the `R²` baseline is the uncentered total sum of squares `Σ τ²` (the
error of the zero predictor), not the mean-centered one.
"""
function r2_torque(f::SCEFit)::Float64
    ss_res = rss_torque(f)          # validates the dataset has torque (throws if not)
    ss_tot = sum(abs2, f.dataset.y_T)  # no-intercept model ⇒ uncentered baseline
    return ss_tot == 0 ? 1.0 : 1 - ss_res / ss_tot
end

"""
    rmse_torque(f::SCEFit) -> Float64

In-sample torque root-mean-square error over the flattened per-atom torque components
(eV). Requires a torque-carrying dataset.
"""
rmse_torque(f::SCEFit)::Float64 = sqrt(rss_torque(f) / length(f.dataset.y_T))
