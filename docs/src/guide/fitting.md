# Data and fitting

```@meta
CurrentModule = SCEFitting
```

With a basis in hand, the second half of the workflow pairs it with DFT data, fits the
coefficients, and reports how well the fit did. The objects are [`SCEDataset`](@ref),
[`fit`](@ref) with a pluggable estimator, and [`SCEFit`](@ref) / [`SCEPredictor`](@ref).

## Datasets

A [`SCEDataset`](@ref) materializes the energy design matrix `X_E[config, salc] =
Φ_salc(config)` from a list of spin configurations (each `3 × n_atoms`, unit columns) and
their energies:

```julia
dataset = SCEDataset(basis, configs, energies)
```

The four-argument form additionally takes per-configuration torques and builds the torque
design matrix for an energy + torque co-fit:

```julia
dataset = SCEDataset(basis, configs, energies, torques)   # torques: each 3 × n_atoms
```

You can also go straight from a DFT source (see [Persistence and I/O](io.md)):
`SCEDataset(basis, src)`.

## The fit and the analytic intercept

[`fit`](@ref) column-centers the design matrix, so the reference energy `j0` is recovered
analytically as `mean(y_E − X_E·jϕ)` — *independent of the estimator* — and the centered
problem is handed to [`solve_coefficients`](@ref). Every estimator therefore returns only
the slope coefficients and adds no intercept of its own.

```julia
f  = fit(SCEFit, dataset, OLS())
J  = coef(f)         # the fitted coefficients jϕ, in SALCKey order
j0 = intercept(f)    # the reference energy
```

## Energy + torque co-fit

The SCE's second observable is the per-atom torque
``\boldsymbol\tau_a = -\hat{\boldsymbol e}_a \times \partial E / \partial \hat{\boldsymbol e}_a``
(the Landau–Lifshitz / physical torque ``\boldsymbol m_a \times \boldsymbol B_{\mathrm{eff},a}``).
Because [`predict_torque`](@ref) is the analytic gradient of the same surface
[`predict_energy`](@ref) evaluates, the two are consistent by construction. A
`torque_weight ∈ (0, 1]` runs a co-fit that minimizes

```math
L = (1 - w)\,\mathrm{MSE}_{\text{energy}} + w\,\mathrm{MSE}_{\text{torque}},
```

implemented by whitening: the centered energy block is row-scaled by ``\sqrt{(1-w)/n_E}``
and the torque block by ``\sqrt{w/n_T}``, then stacked. `j0` does not enter the torque
block, so it stays an energy-only quantity.

```julia
f = fit(SCEFit, SCEDataset(basis, configs, energies, torques), OLS(); torque_weight = 0.5)
(r2_energy(f), r2_torque(f))
```

The [Getting started](../getting_started.md#Add-the-torque) page runs a full co-fit
end-to-end.

## Estimators

The estimator is the regression strategy, dispatched on [`AbstractEstimator`](@ref). Three
are in-tree (closed form, no dependencies):

| Estimator | Penalty | Notes |
|-----------|---------|-------|
| [`OLS`](@ref) | none | ordinary least squares (QR) |
| [`Ridge`](@ref) | ``\lambda\lVert\beta\rVert_2^2`` | L2, closed form |
| [`AdaptiveRidge`](@ref) | ``\lambda\sum_j w_j\beta_j^2`` | iterative reweighted ridge, an L0 approximation |

[`AdaptiveRidge`](@ref) (Frommlet & Nuel 2016) repeatedly refits a per-coefficient
weighted ridge with ``w_j = 1/(\beta_j^2 + \varepsilon)``, so large coefficients get a
light penalty and small ones a heavy penalty — iterating drives the small ones toward
zero. Each subproblem is the analytic weighted ridge, so it needs no extension:

```julia
fit(SCEFit, dataset, AdaptiveRidge(lambda = 1e-3))     # L0-like selection, closed form
```

The penalized-path estimators — the Lasso, the elastic net, and the adaptive Lasso — are
provided by a **GLMNet extension** that lights up under `using GLMNet`:

```julia
using GLMNet                               # activates the estimator extension

fit(SCEFit, dataset, Lasso())                          # CV-selected λ, sparse model
fit(SCEFit, dataset, Lasso(select = :lambda_1se))      # the parsimonious 1-SE model
fit(SCEFit, dataset, ElasticNet(alpha = 0.5))          # an L1/L2 mix
fit(SCEFit, dataset, Lasso(lambda = 1e-3))             # a fixed penalty (no CV)
fit(SCEFit, dataset, AdaptiveLasso())                  # data-driven reweighted Lasso
```

[`ElasticNet`](@ref) minimizes GLMNet's
``\tfrac{1}{2n}\lVert y - X\beta\rVert^2 + \lambda[(1-\alpha)/2\,\lVert\beta\rVert_2^2 +
\alpha\,\lVert\beta\rVert_1]`` on the centered design (so `j0` stays analytic), with
column standardization. With `lambda = nothing` the penalty is chosen by `nfolds`-fold
cross-validation (`select = :lambda_min` or `:lambda_1se`); a numeric `lambda` fits at
exactly that penalty. For a co-fit the CV folds are **grouped by configuration** so that a
configuration's energy and torque rows never split across folds.

[`AdaptiveLasso`](@ref) (Zou 2006) runs a `pilot` estimator first, then a weighted Lasso
with per-column penalty factors ``w_j = 1/\max(|\hat\beta_j^{\text{pilot}}|,
\varepsilon)^\gamma`` — penalizing columns the pilot found small and sparing those it
found large, which gives it oracle selection. `gamma = 0` reduces to a plain Lasso; the
pilot defaults to [`OLS`](@ref) but any estimator works, including a
[`PrecomputedPilot`](@ref) that reuses a prior fit's coefficients:

```julia
fit(SCEFit, dataset, AdaptiveLasso(gamma = 1.0))                       # OLS pilot (Zou 2006)
fit(SCEFit, dataset, AdaptiveLasso(pilot = Ridge(lambda = 1e-4)))      # for ill-conditioned designs
fit(SCEFit, dataset, AdaptiveLasso(pilot = PrecomputedPilot(coef(prior)), lambda = 1e-3))
```

It shares `lambda` / `standardize` / CV behavior with [`ElasticNet`](@ref) (`lambda =
nothing` selects λ by configuration-grouped CV with the adaptive weights held fixed).

Implementing your own estimator is one method:

```julia
struct MyEstimator <: AbstractEstimator end
SCEFitting.solve_coefficients(::MyEstimator, X, y; groups = nothing) = X \ y  # centered (X, y)
```

## Refitting on a selected support

After a sparse fit (`Lasso` / `AdaptiveLasso` / `AdaptiveRidge`), the surviving
coefficients are shrunk toward zero by the penalty. [`refit`](@ref) removes that bias: it
keeps the support of an existing fit and re-solves on just those columns — by default with
[`OLS`](@ref), the textbook de-biasing step.

```julia
fsparse = fit(SCEFit, dataset, Lasso())     # selects a support (some jϕ exactly zero)
fdebias = refit(fsparse)                     # OLS on that support — unshrunk survivors
```

A column survives when its scaled-magnitude contribution `|coef(f)[j]|·‖X[:, j]‖` exceeds
`threshold` (default `0`, i.e. exactly the nonzero support); pass a positive `threshold` to
prune further. `refit` reuses the dataset and `torque_weight` of the input fit, so the
co-fit whitening is identical.

## Diagnostics

A fitted [`SCEFit`](@ref) answers the usual questions:

```julia
r2_energy(f);  rmse_energy(f)        # in-sample energy R² / RMSE
r2_torque(f);  rmse_torque(f)        # torque equivalents (need a co-fit dataset)
nobs(f)                              # number of energy observations
dof(f)                               # degrees of freedom: length(coef(f)) + 1
rss_energy(f);  rss_torque(f)        # residual sums of squares
residuals_energy(f);  residuals_torque(f)   # the raw residual vectors
```

The energy and torque blocks are reported separately throughout (the rebuild does not fold
them into one combined residual): `residuals_energy(f)` is `y_E − (j0 + X_E·jϕ)` and
`residuals_torque(f)` is `y_T − X_T·jϕ` over the flattened torque components.

To inspect the coefficients as a table, use [`coeftable`](@ref) — a Tables.jl source with
one row per SALC. See [Persistence and I/O](io.md#Tabular-coefficients).

Next: [Persistence and I/O](io.md).
