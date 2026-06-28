# Data and fitting

```@meta
CurrentModule = MagestyRebuild
```

With a basis in hand, the second half of the workflow pairs it with DFT data, fits the
coefficients, and reports how well the fit did. The objects are [`SCEDataset`](@ref),
[`fit`](@ref) with a pluggable estimator, and [`SCEFit`](@ref) / [`SCEModel`](@ref).

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

The estimator is the regression strategy, dispatched on [`AbstractEstimator`](@ref). Two
are in-tree (closed form, no dependencies):

| Estimator | Penalty | Notes |
|-----------|---------|-------|
| [`OLS`](@ref) | none | ordinary least squares (QR) |
| [`Ridge`](@ref) | ``\lambda\lVert\beta\rVert_2^2`` | L2, closed form |

Two more — the Lasso and the elastic net, with cross-validated regularization paths — are
provided by a **GLMNet extension** that lights up under `using GLMNet`:

```julia
using GLMNet                               # activates the estimator extension

fit(SCEFit, dataset, Lasso())                          # CV-selected λ, sparse model
fit(SCEFit, dataset, Lasso(select = :lambda_1se))      # the parsimonious 1-SE model
fit(SCEFit, dataset, ElasticNet(alpha = 0.5))          # an L1/L2 mix
fit(SCEFit, dataset, Lasso(lambda = 1e-3))             # a fixed penalty (no CV)
```

[`ElasticNet`](@ref) minimizes GLMNet's
``\tfrac{1}{2n}\lVert y - X\beta\rVert^2 + \lambda[(1-\alpha)/2\,\lVert\beta\rVert_2^2 +
\alpha\,\lVert\beta\rVert_1]`` on the centered design (so `j0` stays analytic), with
column standardization. With `lambda = nothing` the penalty is chosen by `nfolds`-fold
cross-validation (`select = :lambda_min` or `:lambda_1se`); a numeric `lambda` fits at
exactly that penalty. For a co-fit the CV folds are **grouped by configuration** so that a
configuration's energy and torque rows never split across folds.

Implementing your own estimator is one method:

```julia
struct MyEstimator <: AbstractEstimator end
MagestyRebuild.solve_coefficients(::MyEstimator, X, y; groups = nothing) = X \ y  # centered (X, y)
```

## Diagnostics

A fitted [`SCEFit`](@ref) answers the usual questions:

```julia
r2_energy(f);  rmse_energy(f)        # in-sample energy R² / RMSE
r2_torque(f);  rmse_torque(f)        # torque equivalents (need a co-fit dataset)
nobs(f)                              # number of energy observations
```

To inspect the coefficients as a table, use [`coeftable`](@ref) — a Tables.jl source with
one row per SALC. See [Persistence and I/O](io.md#Tabular-coefficients).

Next: [Persistence and I/O](io.md).
