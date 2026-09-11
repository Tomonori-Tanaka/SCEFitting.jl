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
`SCEDataset(basis, src)`. On that path every atom the SALC basis references must
carry a nonzero magnetic moment in every configuration — a quenched moment would
enter the fit through a placeholder direction and silently bias it, so it is an
error. Species that are genuinely non-magnetic belong outside the basis
(`lmax = 0`), and their moments are then never consulted.

### Slicing and concatenation

Datasets slice by configuration and concatenate without recomputing design
matrices, which makes train/test splits, filtering, and incremental data addition
cheap:

```julia
train, test = dataset[1:80], dataset[81:end]   # ranges, index vectors, Bool masks, :
sce_fit = fit(SCEFit, train, OLS())
rmse_holdout = sqrt(sum(abs2, predict_energy(sce_fit, test.configs) .- test.y_E) / length(test))

more = SCEDataset(basis2, new_configs, new_energies)
dataset = vcat(dataset, more)     # basis2 must be the same basis (fingerprint-checked)
```

`vcat` accepts parts built on a persisted-and-reloaded basis (the check is the
SALC-basis fingerprint, not object identity); mixing torque-bearing and
energy-only datasets is an error.

## The fit and the analytic intercept

[`fit`](@ref) column-centers the design matrix, so the reference energy `j0` is recovered
analytically as `mean(y_E − X_E·jϕ)` — *independent of the estimator* — and the centered
problem is handed to [`solve_coefficients`](@ref). Every estimator therefore returns only
the slope coefficients and adds no intercept of its own.

```julia
sce_fit = fit(SCEFit, dataset, OLS())
J       = coef(sce_fit)         # the fitted coefficients jϕ, one per design column
                                # (= per SALC in SALCKey order, unless cross-orbit alias
                                # groups are tied — then SCEPredictor(sce_fit) expands them)
j0      = intercept(sce_fit)    # the reference energy
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
sce_fit = fit(SCEFit, SCEDataset(basis, configs, energies, torques), OLS(); torque_weight = 0.5)
(r2_energy(sce_fit), r2_torque(sce_fit))
```

The [Getting started](../getting_started.md#Add-the-torque) page runs a full co-fit
end-to-end.

## Estimators

The estimator is the regression strategy, dispatched on [`AbstractEstimator`](@ref). Four
are in-tree (closed form, no dependencies):

| Estimator | Penalty | Notes |
|-----------|---------|-------|
| [`OLS`](@ref) | none | ordinary least squares (QR) |
| [`Ridge`](@ref) | ``\lambda\sum_j m_j\beta_j^2`` | L2, closed form |
| [`AdaptiveRidge`](@ref) | ``\lambda\sum_j D_j\beta_j^2`` | iterative reweighted ridge, an L0 approximation |
| [`GroupAdaptiveRidge`](@ref) | as above, ``w`` shared per group | group-L0, fixed cost weights |

[`AdaptiveRidge`](@ref) (Frommlet & Nuel 2016) repeatedly refits a per-coefficient
weighted ridge with ``D_j = m_j/(m_j\beta_j^2 + \varepsilon)``, so large coefficients
get a light penalty and small ones a heavy penalty — iterating drives the small ones
toward
zero. Each subproblem is the analytic weighted ridge, so it needs no extension:

```julia
fit(SCEFit, dataset, AdaptiveRidge(lambda = 1e-3))     # L0-like selection, closed form
```

[`GroupAdaptiveRidge`](@ref) is its group extension: all columns of a group share one
weight ``w_j = v_g/(\sum_{k\in g} m_k\beta_k^2 + p_g\varepsilon)``, so whole groups —
not individual columns — are driven to zero, and the fixed multiplier ``v_g`` prices
each group (at convergence a surviving group pays exactly ``\lambda v_g``). This is the
estimator behind the Monte-Carlo-cost-aware selection workflow below.

## The penalty metric

``m_j`` above is the **penalty metric**: a per-column scale, `nothing` (uniform) or a
vector from [`penalty_metric`](@ref). It is on by default whenever an estimator is
built from a basis.

``\lambda\lVert\beta\rVert^2`` is not invariant under rescaling a design column, and
SALC column norms are set by basis conventions — an orbit's member count, and the
ordering multiplicity the member fold absorbs — not by physics. A column with a larger
norm carries a smaller coefficient at the same physical effect, so it is shrunk *less*:
the unweighted penalty quietly prefers large orbits and high body order. The metric

```math
m_j(w) = (1-w)\,\mathrm{Var}[\Phi_j]
       + w\,\frac{1}{3 n_\text{atoms}}\,
         \mathbb{E}\Bigl[\sum_a \lVert (\partial\Phi_j/\partial e_a)\times e_a\rVert^2\Bigr]
```

— the reference norm of the column as the estimator sees it, over uniform-random spin
configurations — removes that, leaving only the prior you state deliberately through
`theta`.

```julia
estimator = GroupAdaptiveRidge(basis; lambda = 1e-5, theta = 1.0)   # metric attached
m   = penalty_metric(basis; torque_weight = 0.3)              # or build it yourself
est_w = GroupAdaptiveRidge(basis; lambda = 1e-5, torque_weight = 0.3)
fit(SCEFit, dataset, est_w; torque_weight = 0.3)              # weights must agree
est_plain = GroupAdaptiveRidge(salc_groups(basis), weights; lambda = 1e-5)  # uniform

# a λ sweep: build the metric ONCE and move it along the path
fits = [fit(SCEFit, dataset, SCEFitting.with_lambda(estimator, l)) for l in lambdas]
```

Use `SCEFitting.with_lambda` rather than rebuilding the estimator by hand. A hand
rebuild (`GroupAdaptiveRidge(estimator.column_groups, estimator.group_weights; lambda = l)`)
silently drops the metric, and a dropped metric is indistinguishable from a deliberate
uniform one — no door will complain. Rebuilding through the basis-aware constructor
instead re-runs `penalty_metric` at every point, which is the expensive half: the
metric is a Monte-Carlo average — ~3.3 % relative standard error per column at the
default `metric_nconfig = 2048` (5.3 % on the worst column, `1/√nconfig` from there),
costing 0.56 s on a 54-atom, 21-column basis and rising linearly in both. The λ path
itself is unaffected — the metric is one extra multiply per column per iteration.

Three things follow from that definition and are worth knowing:

- **`torque_weight` is part of the metric.** The assembled design mixes the energy and
  torque blocks by ``w``, so the column scales move with it. Build the metric at the
  weight you will fit at; [`fit`](@ref), [`select_fit`](@ref) and
  [`cross_validate`](@ref) refuse a mismatch, along with a metric built on a different
  basis. (A metric you assembled by hand carries no provenance and is never second-
  guessed.)
- **It is a property of the basis, not of the training data.** That is deliberate: λ
  becomes comparable between cells, cross-validation needs no per-fold recomputation to
  stay leak-free, and the prior sits on the function space rather than on how strongly
  one training set happened to excite each column. The price is the converse — a column
  the reference ensemble excites weakly but your data drives hard is effectively
  under-penalized — worth keeping in mind when your configurations are far from
  uniform (near-collinear low-temperature states, say) rather than spread over the
  sphere the way the reference ensemble is.
- **The support rule is a separate scale, and already invariant.** [`refit`](@ref) and
  [`select_support`](@ref) threshold ``|\beta_j|\cdot\lVert X[:,j]\rVert``, which does
  not move under a column rescaling at all, so the metric does not touch it. The two
  knobs are independent by construction rather than by tuning.

An entry of exactly `0` marks a column **unpenalized** — that is how the moment
channel's ``\mu_0`` intercepts are kept out of the penalty (see the
[moment guide](moment.md)), for every estimator rather than only the group form.

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

A column survives when its scaled-magnitude contribution `|coef(sce_fit)[j]|·‖X[:, j]‖` exceeds
`threshold` (default `0`, i.e. exactly the nonzero support); pass a positive `threshold` to
prune further. `refit` reuses the dataset and `torque_weight` of the input fit, so the
co-fit whitening is identical.

## Cost-weighted group selection

A Monte-Carlo sweep over a fitted SCE pays per **contraction entry**, and an entry
vanishes only when *every* SALC of its `(body, orbit, l-multiset)` group has a zero
coefficient — so the quantity to minimize alongside the fit error is the summed cost of
the surviving groups, not the coefficient count. The workflow prices each group up
front and lets the penalty act at exactly that granularity:

```julia
estimator = GroupAdaptiveRidge(basis; lambda = 1.0, theta = 1.0)  # cost-proportional
path      = select_fit(dataset, estimator; lambdas = 10.0 .^ range(2, -8; length = 25))
fbest     = refit(path.fit; threshold = path.threshold)  # de-bias the alive support
```

(The reweighted ridge crushes dead groups to tiny — not exactly zero — values, so
"alive" is decided by a relative floor on the scaled magnitudes; `path.threshold` is
the effective absolute value at the selected λ, and passing it to [`refit`](@ref)
realizes exactly the support the path reported.)

The convenience constructor bundles `SCEFitting.salc_groups` (the column → group
labels) with `SCEFitting.cost_weights`, which sets the fixed per-group weights to
``v_g = \sqrt{p_g}\,(c_g/\bar c)^\theta`` — ``c_g`` being the group's a-priori
Monte-Carlo cost (its distinct contraction entries, `SCEFitting.group_costs`). Two
independent knobs shape the cost–accuracy trade:

- **`theta ∈ [0, 1]`** tilts the *penalty*: `theta = 0` ignores cost (plain group
  selection), `theta = 1` makes an expensive group earn its keep with a
  correspondingly larger error reduction. Different `theta` change the *order* in
  which groups die along the λ path.
- **`delta ≥ 0`** tilts the *selection*: [`select_fit`](@ref) scores every λ (GCV by
  default, or configuration-grouped CV with `criterion = :cv`) and picks the
  **cheapest** λ whose score is within `(1 + delta)` of the path minimum — the
  cost-aware generalization of the conventional `:lambda_1se` rule.

The returned [`SelectionPath`](@ref) is a Tables.jl source with one row per λ
(`lambda`, `score`, `edof`, `n_alive`, `cost`, `selected`); sweeping `theta` and taking
the lower envelope of the per-θ paths traces the full (cost, error) Pareto front. For
a torque co-fit prefer `criterion = :cv` — see the caveat in [`gcv`](@ref).

On real data the group-magnitude spectrum is usually **continuous** — there is no
clean alive/dead gap for the λ path to expose, and most of the cost–error trade
lives in the support threshold itself. [`select_support`](@ref) is the second knob:
it sweeps the alive threshold at a fixed fit, de-biases with [`refit`](@ref) at each
point (one cheap OLS per point — no re-solving the penalty), scores each refit on an
evaluation dataset, and applies the same Pareto rule:

```julia
train, held = dataset[1:80], dataset[81:100]        # dataset slicing
sce_fit = fit(SCEFit, train, GroupAdaptiveRidge(basis; lambda = 1e-5, theta = 1.0))
front   = select_support(sce_fit; npoints = 25, evalset = held, delta = 0.05)
front.fit                                           # the selected de-biased refit
```

Pass a held-out `evalset` for an honest error axis (the default is the in-sample
training set). On a production Nd₂Fe₁₄B model this front offered, e.g., 38 % of the
Monte-Carlo cost at a held-out torque RMSE *better* than the full model, and 3 % of
the cost at +19 % — trades the λ path alone cannot see.

## Cross-validation

[`cross_validate`](@ref) is the generic, honest assessment behind all of the above:
configuration-grouped K-fold CV of any `fit` call, refitting each fold from scratch
(centering and torque whitening stay inside the training fold — nothing leaks) and
scoring the held-out configurations in prediction space:

```julia
cv_result = cross_validate(dataset, GroupAdaptiveRidge(basis; lambda = 1e-5);
                    torque_weight = 1.0, nfolds = 5)
cv_result.pooled_rmse_energy, cv_result.pooled_rmse_torque   # both error axes, out-of-fold
cv_result.score                                       # per-fold (1−w)·MSE_E + w·MSE_T
```

Both RMSEs are reported whenever the dataset carries torque data, **independent of
`torque_weight`** — an energy-only fit still gets its torque error measured, which
is exactly what comparing `torque_weight` settings needs. The `pooled_*` fields
aggregate the out-of-fold residuals (every configuration is held out exactly once);
the per-fold columns give the spread. Use it where a single train/holdout split is
too noisy — e.g. to back a [`select_support`](@ref) point with a K-fold error bar,
or to rank estimators on an equal footing. It differs from
`select_fit(criterion = :cv)`, which whitens globally and only *ranks* a λ path;
`cross_validate` is the generalization-error estimate.

## Diagnostics

A fitted [`SCEFit`](@ref) answers the usual questions:

```julia
r2_energy(sce_fit);  rmse_energy(sce_fit)        # in-sample energy R² / RMSE
r2_torque(sce_fit);  rmse_torque(sce_fit)        # torque equivalents (need a co-fit dataset)
nobs(sce_fit)                              # number of energy observations
dof(sce_fit)                               # degrees of freedom: length(coef(sce_fit)) + 1
rss_energy(sce_fit);  rss_torque(sce_fit)        # residual sums of squares
residuals_energy(sce_fit);  residuals_torque(sce_fit)   # the raw residual vectors
effective_dof(sce_fit)                     # hat-matrix trace + 1 (linear estimators)
gcv(sce_fit)                               # generalized cross-validation score
```

For a *linear* estimator (`islinear`: `OLS` / `Ridge` / `AdaptiveRidge` /
`GroupAdaptiveRidge`) two closed-form model-selection diagnostics come for free:
[`effective_dof`](@ref) is the trace of the hat matrix plus the intercept — the
*effective* parameter count a penalized fit actually spends, as opposed to
[`dof`](@ref)'s raw count — and [`gcv`](@ref) is the generalized cross-validation
score `n·RSS/(n − df)²` built on it, the fast λ-selection criterion used by
[`select_fit`](@ref).

The energy and torque blocks are reported separately throughout (the rebuild does not fold
them into one combined residual): `residuals_energy(sce_fit)` is `y_E − (j0 + X_E·jϕ)` and
`residuals_torque(sce_fit)` is `y_T − X_T·jϕ` over the flattened torque components.

The generic names — [`predict`](@ref), [`residuals`](@ref), [`r2`](@ref), plus `coef` /
`fit` / `nobs` / `dof` / `coeftable` / `islinear` — **extend StatsAPI** (imported, not
shadowed) and default to the energy block, so they compose with `using StatsBase` /
`using GLM` instead of clashing.

To inspect the coefficients as a table, use [`coeftable`](@ref) — a Tables.jl source with
one row per SALC. See [Persistence and I/O](io.md#Tabular-coefficients).

Next: [Persistence and I/O](io.md).
