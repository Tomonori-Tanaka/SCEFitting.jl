"""
    AbstractEstimator

Regression strategy. Implement one by subtyping this and defining

    solve_coefficients(est, X, y) -> Vector{Float64}

where `(X, y)` is already **column-centered** (so the solver adds no intercept;
the reference energy `j0` is recovered analytically downstream). In-tree
estimators: [`OLS`](@ref), [`Ridge`](@ref). GLMNet-backed estimators (Lasso /
ElasticNet) are provided by a package extension.
"""
abstract type AbstractEstimator end

"""
    MetricProvenance

Where a penalty metric came from: the `channel` it was built for (`:energy` or
`:moment`), the `torque_weight` its reference norms were taken at (`0.0` on the moment
channel, which has no torque block), the reference-ensemble size `nconfig` and its
`seed`, and the `fingerprint` of the SALC basis it was built on. Carried on the
estimator so a fit can refuse a metric built for a different problem. The metric enters
every penalized coefficient, and a mismatch is invisible to the scale-invariance gates
— those hold for any `m ∝ c²`, right or wrong — so it has to be checked from outside.
"""
const MetricProvenance = @NamedTuple{channel::Symbol, torque_weight::Float64,
                                     nconfig::Int, seed::Int, fingerprint::UInt64}

# Metric accessors that work for EVERY estimator, so the fitting doors can check
# provenance without dispatching on which estimators happen to carry a metric.
_estimator_metric(::AbstractEstimator)::Union{Nothing,Vector{Float64}} = nothing
_estimator_provenance(::AbstractEstimator)::Union{Nothing,MetricProvenance} = nothing

# A metric built for another problem is invisible to every numerical gate — scale
# invariance holds for any `m ∝ c²`, right or wrong — so the doors check the
# provenance instead. An estimator carrying a metric with NO provenance is a
# hand-built one and passes: the caller owns it.
function _check_metric_provenance(est::AbstractEstimator, channel::Symbol,
                                  fingerprint::UInt64, torque_weight::Float64)
    pv = _estimator_provenance(est)
    pv === nothing && return nothing
    pv.channel === channel || throw(ArgumentError(
        "the estimator's penalty metric was built for the $(pv.channel) channel, but " *
        "this is the $channel channel. The two channels have different design rows " *
        "and different reference norms; build the metric from THIS basis."))
    pv.fingerprint == fingerprint || throw(ArgumentError(
        "the estimator's penalty metric was built on a different basis (SALC " *
        "fingerprint mismatch). Column scales are basis-specific; rebuild it with " *
        "`penalty_metric` on this basis."))
    pv.torque_weight == torque_weight || throw(ArgumentError(
        "the estimator's penalty metric was built at torque_weight = " *
        "$(pv.torque_weight), but this fit uses $torque_weight. The assembled design " *
        "mixes the energy and torque blocks by that weight, so the column scales " *
        "move with it; rebuild the metric at the weight you are fitting at."))
    return nothing
end

# --- penalty metric: shared validation and use ------------------------------------

# Validate a penalty metric at construction. `nothing` means "uniform" and is kept as a
# distinct value from an all-ones vector so `show`, the provenance check, and the
# byte-equality gates can tell "no metric" from "a metric that happens to be flat".
# Length is NOT checked here — a constructor does not see the design matrix, so that
# check lives at the `solve_coefficients` door.
function _validated_metric(metric)::Union{Nothing,Vector{Float64}}
    metric === nothing && return nothing
    m = Vector{Float64}(metric)
    isempty(m) &&
        throw(ArgumentError("metric must be nonempty (one entry per design column)"))
    all(x -> isfinite(x) && x >= 0, m) ||
        throw(ArgumentError("metric entries must be finite and ≥ 0"))
    any(>(0), m) || throw(ArgumentError(
        "metric is zero on every column, so no column would be penalized — that is " *
        "`OLS()`, not a penalized fit. Use `OLS()` if that is the intent."))
    return m
end

# The per-column penalty scale the solvers see: the estimator's metric, or the uniform
# metric when none was supplied. `1.0 * x === x` for every Float64, so a `nothing`
# metric reproduces the unweighted penalty bit for bit and the solvers below need no
# separate no-metric branch.
function _metric_vector(metric::Union{Nothing,Vector{Float64}}, p::Int,
                        what::AbstractString)::Vector{Float64}
    metric === nothing && return ones(Float64, p)
    length(metric) == p || throw(DimensionMismatch(
        "$what metric length $(length(metric)) does not match design-matrix column " *
        "count $p; the metric was likely built on a different basis."))
    return metric
end

# Unpenalized columns (`mⱼ == 0`) leave the normal equations singular unless that block
# is itself of full column rank: `v'(X'X + λD)v = ‖Xv‖² + λΣ Dⱼvⱼ²` vanishes only for
# `v_pen = 0` and `X_free v_free = 0`, so positive definiteness is exactly a rank
# condition on the unpenalized block. A `Symmetric` solve on a singular matrix does not
# reliably throw — it can return numerical garbage instead (measured at ‖β‖ ~ 1e16
# against OLS's 0.48 on a design with one duplicated column) — so refuse by name first.
#
# The test runs on the free block's Gram, whose eigenvalues are the SQUARED singular
# values of `X_free`, and cuts at `√eps` relative: the normal equations square the
# condition number, so `κ(X_free) ≳ 1e8` is the point past which the unpenalized
# coefficients stop being meaningful whether or not the factorization reports failure.
# A plain `cholesky(...; check = false)` is not enough — on an exactly duplicated
# column the trailing pivot lands a rounding step above zero and succeeds. The cut is
# deliberately stricter than the rank test `_edof_free` applies to the same block, so a
# solve that is accepted here is never rejected by the diagnostic afterwards.
function _check_free_block(XtX::AbstractMatrix, metric::Vector{Float64},
                           what::AbstractString)
    free = findall(iszero, metric)
    isempty(free) && return nothing
    ev = eigvals(Symmetric(Matrix{Float64}(XtX[free, free])))
    if !(ev[end] > 0 && ev[1] > sqrt(eps(Float64)) * ev[end])
        throw(ArgumentError(
            "$what: the unpenalized columns $free are linearly dependent (or nearly " *
            "so) on this design — the block's condition number is " *
            "$(ev[1] > 0 ? sqrt(ev[end] / ev[1]) : Inf), against a limit of 1e8 — so " *
            "the penalized normal equations are singular and those coefficients are " *
            "not identified. Penalize the columns, or drop the dependent ones."))
    end
    return nothing
end

# Relative ∞-norm change of an IRLS iterate, measured in metric coordinates `√mⱼ·βⱼ`
# and over the PENALIZED columns only. `√mⱼ·βⱼ` is the invariant of a column rescaling
# (`Φⱼ → cⱼΦⱼ` sends `βⱼ → βⱼ/cⱼ` and `mⱼ → cⱼ²mⱼ`), so the stopping rule moves with the
# penalty instead of against it; restricting to penalized columns keeps one large
# unpenalized coefficient from turning the relative tolerance into an absolute one far
# coarser than the coefficients it is meant to converge. With a uniform metric and no
# unpenalized column this is exactly the plain relative ∞-norm rule.
function _irls_rel_change(beta_new::Vector{Float64}, beta::Vector{Float64},
                          metric::Vector{Float64})::Float64
    num = 0.0
    den = 0.0
    @inbounds for j in eachindex(beta_new)
        metric[j] > 0.0 || continue
        s = sqrt(metric[j])
        d = s * abs(beta_new[j] - beta[j])
        d > num && (num = d)
        a = s * abs(beta_new[j])
        a > den && (den = a)
    end
    return num / max(den, eps(Float64))
end

"""
    OLS()

Ordinary least squares, solved by pivoted QR (`X \\ y`) — **not** the normal
equations, which would square the condition number. On a rank-deficient design this
returns the minimum-norm solution rather than an arbitrary one, which is why every
other closed-form estimator routes here at zero penalty.
"""
struct OLS <: AbstractEstimator end

"""
    Ridge(; lambda = 0.0, metric = nothing, metric_provenance = nothing)

L2-penalized least squares with penalty `lambda ≥ 0`. `lambda = 0` routes to
[`OLS`](@ref)'s QR path, exactly (not within a tolerance).

For `lambda > 0` this solves the normal equations `(XᵀX + λ·Diagonal(m))β = Xᵀy`, which
is safe because the penalty makes the matrix positive definite, but which does square
the condition number of `X`. On an ill-conditioned design prefer a larger `lambda`, or
`OLS` if the intent was no penalty at all.

`metric` is the per-column penalty scale `m` (`nothing` = uniform). It makes the
penalty invariant under rescaling a design column — `‖β‖²` on its own is not, and SALC
column norms are set by basis conventions rather than by physics — and an entry of
exactly `0` marks a column as **unpenalized**. Use
[`penalty_metric`](@ref) to build one from a basis; `metric_provenance` records what it
was built from so a fit can refuse a metric made for a different problem. The
unpenalized block must have full column rank or the solve is refused by name.
"""
struct Ridge <: AbstractEstimator
    lambda::Float64
    metric::Union{Nothing,Vector{Float64}}
    metric_provenance::Union{Nothing,MetricProvenance}

    function Ridge(lambda::Real, metric, metric_provenance)
        (lambda >= 0 && isfinite(lambda)) ||
            throw(ArgumentError("lambda must be finite and ≥ 0; got $lambda"))
        return new(Float64(lambda), _validated_metric(metric), metric_provenance)
    end
end
Ridge(lambda::Real) = Ridge(lambda, nothing, nothing)
Ridge(; lambda::Real = 0.0, metric = nothing, metric_provenance = nothing) =
    Ridge(lambda, metric, metric_provenance)

"""
    ElasticNet(; alpha = 1.0, lambda = nothing, standardize = true, nfolds = 10,
               select = :lambda_min, seed = 1, nlambda = 100)

Elastic-net (mixed L1/L2) regression, backed by GLMNet. The method lights up only
when GLMNet is loaded (`using GLMNet` activates the extension); the type itself
lives in the core package so it can be named and dispatched on without the heavy
dependency.

`alpha ∈ [0, 1]` is the elastic-net mixing — `alpha = 1` is the Lasso (pure L1),
`alpha = 0` is ridge-like (pure L2). GLMNet minimizes

    (1 / 2n)·‖y − Xβ‖² + λ·[ (1 − α)/2·‖β‖₂² + α·‖β‖₁ ]

over the **column-centered** `(X, y)` it is handed (see [`AbstractEstimator`](@ref));
no intercept is fitted, so `j0` is still recovered analytically downstream. Note the
`1/2n` and the L1 term: this `lambda` is on a different scale than [`Ridge`](@ref)'s
plain `λ‖β‖₂²`, and is **not** directly comparable.

`standardize = true` (the conventional Lasso default) scales each column to unit
variance before fitting, so with it the penalty effectively acts as `λ·std(colⱼ)` per
column — the L1 term treats differently-scaled SALC columns evenhandedly. Coefficients
are always returned on the original scale.

Penalty strength:
- `lambda = nothing` (default): chosen by `nfolds`-fold cross-validation over
  GLMNet's automatic λ path (the fold count is capped at `n ÷ 3` for small samples).
  `select` picks which λ — `:lambda_min` (lowest CV error) or `:lambda_1se` (the
  sparsest model within one standard error of the minimum, the conventional
  parsimonious choice). `seed` fixes the fold assignment so the fit is reproducible
  within a Julia session/version. `nlambda` sizes the λ grid.
- `lambda = ⟨number⟩`: fit at exactly that penalty, skipping CV.

For an energy+torque co-fit (see [`fit`](@ref)) the CV folds are **grouped by
configuration** — a configuration's energy row and all its torque-component rows go
into the same fold — so resampling never splits one configuration. Because the design
is the whitened, stacked energy+torque matrix, the per-column standardization (hence
the effective penalty on each SALC) depends on `torque_weight`; `lambda` is therefore
not directly comparable between an energy-only fit and a co-fit.
"""
struct ElasticNet <: AbstractEstimator
    alpha::Float64
    lambda::Union{Nothing,Float64}
    standardize::Bool
    nfolds::Int
    select::Symbol
    seed::Int
    nlambda::Int
end
function ElasticNet(; alpha::Real = 1.0, lambda::Union{Nothing,Real} = nothing,
                    standardize::Bool = true, nfolds::Integer = 10,
                    select::Symbol = :lambda_min, seed::Integer = 1,
                    nlambda::Integer = 100)
    (0 <= alpha <= 1) || throw(ArgumentError("alpha must be in [0, 1]; got $alpha"))
    (lambda === nothing || (lambda >= 0 && isfinite(lambda))) ||
        throw(ArgumentError("lambda must be ≥ 0 (or nothing for CV); got $lambda"))
    select in (:lambda_min, :lambda_1se) ||
        throw(ArgumentError("select must be :lambda_min or :lambda_1se; got :$select"))
    nfolds >= 2 || throw(ArgumentError("nfolds must be ≥ 2; got $nfolds"))
    nlambda >= 1 || throw(ArgumentError("nlambda must be ≥ 1; got $nlambda"))
    return ElasticNet(Float64(alpha), lambda === nothing ? nothing : Float64(lambda),
                      standardize, Int(nfolds), select, Int(seed), Int(nlambda))
end

"""
    Lasso(; kwargs...)

The Lasso (L1-penalized) estimator: [`ElasticNet`](@ref) with `alpha = 1`. Takes
the same keyword arguments except `alpha`.
"""
Lasso(; kwargs...) = ElasticNet(; alpha = 1.0, kwargs...)

"""
    AdaptiveLasso(; pilot = OLS(), lambda = nothing, gamma = 1.0,
                  epsilon = eps(Float64), standardize = true, nfolds = 10,
                  select = :lambda_min, seed = 1, nlambda = 100)

One-shot Adaptive Lasso (Zou 2006), backed by GLMNet — the solve lights up only
when GLMNet is loaded (`using GLMNet` activates the extension); the type lives in
the core package so it can be named and dispatched on without the heavy dependency.

The fit is two-stage. First `pilot` (any [`AbstractEstimator`](@ref)) is run on the
same centered `(X, y)` to get a pilot coefficient vector `β̂`. Then a **weighted**
Lasso is solved with per-column penalty factors

    wⱼ = 1 / max(|β̂ⱼ|, epsilon)^gamma

so columns the pilot found large are penalized little and columns it found small are
penalized heavily — the data-driven reweighting that gives the Adaptive Lasso its
oracle selection property. `gamma = 0` reduces to a plain Lasso (uniform weights);
`gamma = 1` (the default) is the Zou 2006 / ALAMODE recipe. `epsilon > 0` floors the
pilot magnitude so a numerically-zero pilot coefficient gives a large-but-finite
penalty rather than `Inf`. GLMNet internally rescales the weights to sum to the
column count, so `lambda` is **not** directly comparable to a plain [`Lasso`](@ref)'s
once `gamma > 0`. The pilot is fit on the un-standardized design and `standardize =
true` then penalizes in the standardized coordinate, so the effective per-column weight
carries an extra `std(colⱼ)` factor on top of `1/|β̂ⱼ|^γ` — the ALAMODE convention,
consistent with [`ElasticNet`](@ref) / [`Lasso`](@ref) here (set `standardize = false`
for the textbook Zou form).

`pilot = OLS()` matches Zou 2006. For a rank-deficient or near-collinear design (e.g.
`n_salcs ≥ nconfig` at `torque_weight = 0`) prefer `pilot = Ridge(lambda = small)`: the
OLS minimum-norm solution puts ~1e-10 noise in the null-space directions, which is
not floored by `epsilon = eps(Float64)` and miscalibrates those weights (the pilot
solve now also warns about the deficiency). Any pilot is
allowed, including a cross-validating [`ElasticNet`](@ref) or a [`PrecomputedPilot`](@ref)
reusing a prior fit's coefficients; the pilot's own solve receives the co-fit `groups`.

Penalty strength and the remaining keywords (`lambda`, `standardize`, `nfolds`,
`select`, `seed`, `nlambda`) behave exactly as in [`ElasticNet`](@ref): `lambda =
nothing` (default) selects λ by configuration-grouped cross-validation with the
adaptive weights held fixed across the path; a number fits at that fixed penalty.
"""
struct AdaptiveLasso <: AbstractEstimator
    pilot::AbstractEstimator
    lambda::Union{Nothing,Float64}
    gamma::Float64
    epsilon::Float64
    standardize::Bool
    nfolds::Int
    select::Symbol
    seed::Int
    nlambda::Int
end
function AdaptiveLasso(; pilot::AbstractEstimator = OLS(),
                       lambda::Union{Nothing,Real} = nothing, gamma::Real = 1.0,
                       epsilon::Real = eps(Float64), standardize::Bool = true,
                       nfolds::Integer = 10, select::Symbol = :lambda_min,
                       seed::Integer = 1, nlambda::Integer = 100)
    (lambda === nothing || (lambda >= 0 && isfinite(lambda))) ||
        throw(ArgumentError("lambda must be ≥ 0 (or nothing for CV); got $lambda"))
    gamma >= 0 || throw(ArgumentError("gamma must be ≥ 0; got $gamma"))
    epsilon > 0 || throw(ArgumentError("epsilon must be > 0; got $epsilon"))
    select in (:lambda_min, :lambda_1se) ||
        throw(ArgumentError("select must be :lambda_min or :lambda_1se; got :$select"))
    nfolds >= 2 || throw(ArgumentError("nfolds must be ≥ 2; got $nfolds"))
    nlambda >= 1 || throw(ArgumentError("nlambda must be ≥ 1; got $nlambda"))
    return AdaptiveLasso(pilot, lambda === nothing ? nothing : Float64(lambda),
                         Float64(gamma), Float64(epsilon), standardize, Int(nfolds),
                         select, Int(seed), Int(nlambda))
end

"""
    PrecomputedPilot(beta)

Estimator adapter that returns a fixed coefficient vector from `solve_coefficients`,
ignoring `(X, y)` except for a length check against `size(X, 2)`. Intended as an
[`AdaptiveLasso`](@ref) `pilot`: it lets the adaptive reweighting reuse coefficients
from a previous fit (`coef(f)`) instead of running a fresh pilot regression. The
vector is copied at construction, so later mutation of the caller's storage does not
leak in.
"""
struct PrecomputedPilot <: AbstractEstimator
    beta::Vector{Float64}
    PrecomputedPilot(b::AbstractVector{<:Real}) = new(Vector{Float64}(b))
end

"""
    AdaptiveRidge(; lambda, epsilon = 1e-8, max_iter = 50, tol = 1e-6,
                  metric = nothing, metric_provenance = nothing)

Iterative Adaptive Ridge (Frommlet & Nuel 2016). Approximates an L0-penalized fit by
repeatedly refitting a per-coefficient weighted ridge problem

    min_b ‖y − Xb‖² + lambda·Σⱼ wⱼ·bⱼ²

and updating the weights `wⱼ = mⱼ / (mⱼ·bⱼ² + epsilon)` between iterations, where `m`
is the per-column penalty scale ([`Ridge`](@ref)'s `metric`; `nothing` = uniform). The
metric sits in the **denominator** so that `Σ mⱼβⱼ²` — a rescaling invariant — is what
the weight map sees; outside it the estimator would not be scale invariant. Iteration zero
is a plain ridge solve (uniform weights); each subsequent step rebuilds the weights
from the current coefficients, so large coefficients get a light penalty and small
ones a heavy penalty — iterating drives the small ones toward zero (an L0
approximation). Each subproblem is the analytic weighted ridge
`b = (X'X + lambda·Diagonal(w)) \\ (X'y)`, the same closed-form family as
[`Ridge`](@ref); no GLMNet is involved, so this estimator works in the core package
without the extension.

The iteration stops when the relative infinity-norm change in the coefficient vector
drops below `tol`, or after `max_iter` reweighting steps. `lambda = 0` reduces to
[`OLS`](@ref). The penalty acts directly on the coefficients (no `standardize`
keyword): the per-cluster `(4π)^(N/2)` basis scale already places them on the
conventional spin-model scale, so a single `epsilon` is a roughly uniform magnitude
floor across clusters.
"""
struct AdaptiveRidge <: AbstractEstimator
    lambda::Float64
    epsilon::Float64
    max_iter::Int
    tol::Float64
    metric::Union{Nothing,Vector{Float64}}
    metric_provenance::Union{Nothing,MetricProvenance}
end
function AdaptiveRidge(; lambda::Real, epsilon::Real = 1e-8, max_iter::Integer = 50,
                       tol::Real = 1e-6, metric = nothing, metric_provenance = nothing)
    (lambda >= 0 && isfinite(lambda)) ||
        throw(ArgumentError("lambda must be finite and ≥ 0; got $lambda"))
    epsilon > 0 || throw(ArgumentError("epsilon must be > 0; got $epsilon"))
    max_iter >= 1 || throw(ArgumentError("max_iter must be ≥ 1; got $max_iter"))
    tol > 0 || throw(ArgumentError("tol must be > 0; got $tol"))
    return AdaptiveRidge(Float64(lambda), Float64(epsilon), Int(max_iter), Float64(tol),
                         _validated_metric(metric), metric_provenance)
end

"""
    GroupAdaptiveRidge(column_groups, group_weights; lambda, epsilon = 1e-8,
                       max_iter = 50, tol = 1e-6, metric = nothing,
                       metric_provenance = nothing)

Group extension of [`AdaptiveRidge`](@ref): iterative reweighted ridge approximating the
weighted **group-L0** penalty `lambda·Σ_g v_g·1{β_g ≠ 0}`. Each reweighting step solves
the analytic weighted ridge

    min_b ‖y − Xb‖² + lambda·Σⱼ Dⱼ·bⱼ²,
    Dⱼ = mⱼ·wⱼ,   wⱼ = v_g / (Σ_{k∈g} m_k·β_k² + p_g·epsilon)

where `g = column_groups[j]`, `β_g` is the current coefficient sub-vector of group `g`,
`p_g` its column count, `v_g = group_weights[g]` a **fixed** positive per-group weight,
and `m` the per-column penalty scale ([`Ridge`](@ref)'s `metric`; `nothing` = uniform).
The metric sits in the denominator, so the group norm the weight map sees is the
rescaling invariant `Σ m_k β_k²` and the fixed point below is preserved; a column with
`mⱼ = 0` is unpenalized and contributes to neither.

At a converged fixed point the penalty contribution of a surviving group tends to
`lambda·v_g` (a constant per alive group), so `v_g` is exactly the group-L0 weight —
e.g. a Monte-Carlo cost weight from `SCEFitting.cost_weights`. This holds with or
without a metric, which is the reason it sits in the denominator: outside it the
converged contribution would be `lambda·v_g·⟨m⟩_g` and `v_g` would stop being the
group-L0 weight. The denominator equals `p_g·(mean_j m_j·βⱼ² + epsilon)`, so `epsilon`
is a magnitude floor on the **metric-weighted** coefficients, independent of the group
size (the same calibration as [`AdaptiveRidge`](@ref), whose update this reproduces
exactly for singleton groups with unit weights).

`column_groups` labels the design-matrix **columns** with contiguous group ids `1:G`
(every label present; for an SCE fit use `SCEFitting.salc_groups`); it is unrelated to
the per-**row** `groups` keyword of [`solve_coefficients`](@ref), which this estimator
ignores. `group_weights` has length `G`. Both vectors are copied at construction.
`lambda = 0` reduces to [`OLS`](@ref). Like `AdaptiveRidge`, the converged fit is a
linear smoother in the fixed-weight sense ([`islinear`](@ref) is `true`), so
[`gcv`](@ref) applies. Iterating to convergence drives whole groups toward zero; follow
with [`refit`](@ref) to select the support and de-bias the survivors, or drive the λ
path with [`select_fit`](@ref).
"""
struct GroupAdaptiveRidge <: AbstractEstimator
    lambda::Float64
    column_groups::Vector{Int}
    group_weights::Vector{Float64}
    epsilon::Float64
    max_iter::Int
    tol::Float64
    group_sizes::Vector{Int}     # p_g, derived and cached by the inner constructor
    metric::Union{Nothing,Vector{Float64}}
    metric_provenance::Union{Nothing,MetricProvenance}

    # `metric` / `metric_provenance` are POSITIONAL and carry no default on purpose:
    # every internal site that rebuilds an estimator from an existing one (the λ-path
    # driver's cold re-solve, the active-column reduction) must pass them through, and
    # a site that forgets has to fail with a `MethodError` rather than silently return
    # an unweighted fit whose λ means something else.
    function GroupAdaptiveRidge(lambda::Real, column_groups::AbstractVector{<:Integer},
                                group_weights::AbstractVector{<:Real}, epsilon::Real,
                                max_iter::Integer, tol::Real, metric,
                                metric_provenance)
        (lambda >= 0 && isfinite(lambda)) ||
            throw(ArgumentError("lambda must be finite and ≥ 0; got $lambda"))
        epsilon > 0 || throw(ArgumentError("epsilon must be > 0; got $epsilon"))
        max_iter >= 1 || throw(ArgumentError("max_iter must be ≥ 1; got $max_iter"))
        tol > 0 || throw(ArgumentError("tol must be > 0; got $tol"))
        isempty(column_groups) &&
            throw(ArgumentError("column_groups must be nonempty (one label per column)"))
        minimum(column_groups) >= 1 ||
            throw(ArgumentError("column_groups labels must be ≥ 1; got " *
                                "$(minimum(column_groups))"))
        G = Int(maximum(column_groups))
        sizes = zeros(Int, G)
        for g in column_groups
            sizes[g] += 1
        end
        all(>(0), sizes) ||
            throw(ArgumentError("column_groups labels must cover 1:$G with no gaps; " *
                                "empty group(s): $(findall(==(0), sizes))"))
        length(group_weights) == G ||
            throw(ArgumentError("group_weights length $(length(group_weights)) ≠ " *
                                "number of groups $G"))
        all(v -> isfinite(v) && v > 0, group_weights) ||
            throw(ArgumentError("group_weights must be finite and > 0"))
        m = _validated_metric(metric)
        m === nothing || length(m) == length(column_groups) || throw(DimensionMismatch(
            "metric length $(length(m)) ≠ column_groups length " *
            "$(length(column_groups))"))
        return new(Float64(lambda), Vector{Int}(column_groups),
                   Vector{Float64}(group_weights), Float64(epsilon), Int(max_iter),
                   Float64(tol), sizes, m, metric_provenance)
    end
end
function GroupAdaptiveRidge(column_groups::AbstractVector{<:Integer},
                            group_weights::AbstractVector{<:Real};
                            lambda::Real, epsilon::Real = 1e-8, max_iter::Integer = 50,
                            tol::Real = 1e-6, metric = nothing,
                            metric_provenance = nothing)
    return GroupAdaptiveRidge(lambda, column_groups, group_weights, epsilon, max_iter,
                              tol, metric, metric_provenance)
end

# Compact display: the default struct printer would dump the whole `beta` vector
# (hundreds–thousands of entries) and recurse into the nested pilot.
Base.show(io::IO, p::PrecomputedPilot) =
    print(io, "PrecomputedPilot(", length(p.beta), " coefficients)")
Base.show(io::IO, e::AdaptiveLasso) =
    print(io, "AdaptiveLasso(pilot=", e.pilot, ", lambda=", e.lambda, ", gamma=",
          e.gamma, ", epsilon=", e.epsilon, ", standardize=", e.standardize, ")")
Base.show(io::IO, e::GroupAdaptiveRidge) =
    print(io, "GroupAdaptiveRidge(", length(e.column_groups), " columns in ",
          length(e.group_weights), " groups, lambda=", e.lambda,
          ", epsilon=", e.epsilon, ", ", _metric_summary(e.metric), ")")
Base.show(io::IO, e::Ridge) =
    print(io, "Ridge(lambda=", e.lambda, ", ", _metric_summary(e.metric), ")")
Base.show(io::IO, e::AdaptiveRidge) =
    print(io, "AdaptiveRidge(lambda=", e.lambda, ", epsilon=", e.epsilon, ", ",
          _metric_summary(e.metric), ")")

# One-line penalty-metric summary for `show`: how many columns it covers and how many
# of them it leaves unpenalized (the number a reader needs to sanity-check an
# intercept exemption without printing the whole vector).
_metric_summary(m::Nothing) = "metric=uniform"
_metric_summary(m::Vector{Float64}) =
    string("metric=", length(m), " columns, ", count(iszero, m), " unpenalized")

"""
    islinear(est) -> Bool

Whether `est` is a linear (closed-form) estimator — one whose fitted values are
`ŷ = H·y` with a hat matrix `H` independent of `y`. Gates closed-form diagnostics
(e.g. [`gcv`](@ref) / [`effective_dof`](@ref)). [`OLS`](@ref), [`Ridge`](@ref),
[`AdaptiveRidge`](@ref), and [`GroupAdaptiveRidge`](@ref) (the latter two linear in
the converged-weight sense, the standard adaptive-ridge approximation) return `true`;
the penalty-path / pilot estimators ([`ElasticNet`](@ref) / [`Lasso`](@ref) /
[`AdaptiveLasso`](@ref)) are not linear. Extends `StatsAPI.islinear`.
"""
islinear(::AbstractEstimator) = false
islinear(::OLS) = true
islinear(::Ridge) = true
islinear(::AdaptiveRidge) = true
islinear(::GroupAdaptiveRidge) = true

"""
    solve_coefficients(est, X, y; groups = nothing) -> Vector{Float64}

Solve for the SALC coefficients from the (already centered) design matrix `X` and
target `y`. See [`AbstractEstimator`](@ref) for the `(X, y)` contract.

`groups` is an optional per-row label vector (one entry per row of `X`) marking rows
that derive from the same physical sample — in an energy+torque co-fit, a
configuration's energy row and its torque-component rows share a label. Estimators
that resample rows (cross-validating [`ElasticNet`](@ref) / [`Lasso`](@ref)) keep
same-label rows together so the resampling does not leak within-sample structure;
estimators with a closed-form fit ([`OLS`](@ref) / [`Ridge`](@ref)) ignore it.
"""
function solve_coefficients(est::AbstractEstimator, X::AbstractMatrix, y::AbstractVector;
                            groups = nothing)
    error("solve_coefficients has no method for $(typeof(est)); load the backend " *
          "package (e.g. `using GLMNet` for Lasso/ElasticNet).")
end

# Relative conditioning threshold on the pivoted-QR diagonal (|R_ii| / |R|_max).
# This is a one-sided CONDITIONING gate, not an exact rank witness: since
# σ_min ≤ min|R_ii| and max|R_ii| ≤ σ_max, the diagonal ratio is ≥ 1/κ₂(X), so
# the warning CANNOT fire while κ₂(X) < 1e10 — no false alarms on well-conditioned
# designs — and on an exact dependence (a duplicated column; the MnTe bulk-SOC
# design whose silent max|coef| ~1e7 motivated this, measured ratio 1e-12–1e-14)
# it fires reliably. In the band 1e-15 ≲ 1/κ ≲ 1e-10 a design may be full rank
# and merely severely ill-conditioned (graded columns can reach it), and the QR
# solve itself truncates only at its own ~min(m,n)·eps tolerance — which is why
# the message says "or severely ill-conditioned" and does not claim an exact
# rank. Kahan-type pivoted-QR failures only produce false NEGATIVES (near-
# deficiency without a small diagonal), i.e. the pre-gate behavior. Not a
# refusal: OLS still returns the QR solution (the documented min-norm behavior
# on exact deficiency), but silently returning non-unique coefficients is the
# bug this warns about.
const _OLS_RANK_RTOL = 1e-10

function solve_coefficients(::OLS, X::AbstractMatrix, y::AbstractVector;
                            groups = nothing)::Vector{Float64}
    # Explicit pivoted QR — the exact factorization `X \ y` uses for a rectangular
    # dense matrix, reused here for the conditioning check so the solve is not paid
    # twice (verified bitwise-identical to `X \ y`, tall and wide, full-rank and
    # deficient). For a SQUARE X, `\` would route through LU: on a nonsingular
    # square design the QR path differs by ~1e-15 relative (oracle tolerances
    # absorb it), and on a singular one LU returned ~1e15 garbage where QR gives
    # the min-norm solution.
    F = qr(X, ColumnNorm())
    d = min(size(X)...)
    if d > 0
        dmax = maximum(i -> abs(F.R[i, i]), 1:d)
        r = dmax == 0.0 ? 0 : count(i -> abs(F.R[i, i]) > _OLS_RANK_RTOL * dmax, 1:d)
        if r < size(X, 2)
            # maxlog bounds the spam from resampling sweeps (CV folds /
            # select_support points re-solve the same deficient design many
            # times); the interactive logger honors it, @test_logs captures all.
            @warn "OLS design matrix is rank deficient or severely ill-" *
                  "conditioned: only $r of $(size(X, 2)) columns are independent " *
                  "at a 1e-10 diagonal threshold (κ ≳ 1e10). The least-squares " *
                  "solution is non-unique or unstable, so the returned " *
                  "coefficients (and anything read off them — coeftable, " *
                  "bilinear_terms, Sunny export) are one arbitrary " *
                  "representative; predicted energies are still well-defined. " *
                  "Causes: redundant basis functions the SALC reduction could " *
                  "not certify (see the build warning), supercell aliasing " *
                  "across orbits, or too few / degenerate training " *
                  "configurations. Consider Ridge(lambda > 0), more data, or a " *
                  "smaller basis." maxlog = 4
        end
    end
    return F \ y   # QR-based least squares (more robust than the normal equations)
end

function solve_coefficients(est::Ridge, X::AbstractMatrix, y::AbstractVector;
                            groups = nothing)::Vector{Float64}
    # Exact zero only: at lambda = 0 the penalty is gone and `XtX` may be singular, so
    # route to the QR (min-norm) OLS path — the same guard `AdaptiveRidge` and
    # `GroupAdaptiveRidge` carry, and for the same reason. Do NOT widen this to a
    # near-zero tolerance. Without it a rank-deficient design does not throw: the
    # symmetric factorization of a singular matrix returns numerical garbage, measured
    # at ‖β‖ ~ 1e16 against OLS's 0.48 on a design with one duplicated column — silently,
    # from an estimator whose docstring promises OLS at lambda = 0.
    # [Backported from SLCE.jl 896180e.]
    est.lambda == 0.0 && return solve_coefficients(OLS(), X, y)
    m = _metric_vector(est.metric, size(X, 2), "Ridge")
    XtX = X' * X
    _check_free_block(XtX, m, "Ridge")
    return Symmetric(XtX + est.lambda * Diagonal(m)) \ (X' * y)
end

function solve_coefficients(est::AdaptiveRidge, X::AbstractMatrix, y::AbstractVector;
                            groups = nothing)::Vector{Float64}
    # Exact zero only: at lambda = 0 the penalty diagonal is gone and `XtX` may be
    # singular, so route to the QR (min-norm) OLS path. Any lambda > 0 makes
    # `XtX + lambda·Diagonal(w)` SPD below, so the symmetric solve is safe — do NOT
    # widen this to a near-zero tolerance.
    est.lambda == 0.0 && return solve_coefficients(OLS(), X, y)
    XtX = X' * X
    Xty = X' * y
    p = size(XtX, 2)
    m = _metric_vector(est.metric, p, "AdaptiveRidge")
    _check_free_block(XtX, m, "AdaptiveRidge")
    # Iteration 0: the fixed-weight ridge `Dⱼ = mⱼ` (with a uniform metric, numerically
    # Ridge(lambda)). Only the penalty diagonal changes between iterations; the Gram
    # matrix `XtX` is formed once. The cold start is built from the metric for the same
    # reason the weight map carries it: a metric-blind iteration 0 would make the
    # converged point of this non-convex surrogate depend on the column scaling even
    # though its objective does not.
    beta = Symmetric(XtX + est.lambda * Diagonal(m)) \ Xty
    D = Vector{Float64}(undef, p)
    for _ = 1:est.max_iter
        @. D = m / (m * beta^2 + est.epsilon)
        beta_new = Symmetric(XtX + est.lambda * Diagonal(D)) \ Xty
        rel = _irls_rel_change(beta_new, beta, m)
        beta = beta_new
        rel < est.tol && break
    end
    return beta
end

# The group-adaptive-ridge weight map — the ONLY definition of the group-form update
# `wⱼ = v_g / (Σ_{k∈g} m_k·β_k² + p_g·ε)`. What it fills in and returns is the PENALTY
# DIAGONAL `Dⱼ = mⱼ·wⱼ`, not `wⱼ`: the metric belongs to the penalty, and every
# consumer (the solvers, and `_penalty_diagonal` behind the GCV / effective-dof
# diagnostics in fitting/selection.jl) wants the diagonal. Change the formula here and
# nowhere else. `normsq` is a reusable length-`G` accumulator; fills `D` in place.
function _gar_weights!(D::Vector{Float64}, beta::Vector{Float64},
                       column_groups::Vector{Int}, group_weights::Vector{Float64},
                       group_sizes::Vector{Int}, metric::Vector{Float64},
                       epsilon::Float64,
                       normsq::Vector{Float64})::Vector{Float64}
    fill!(normsq, 0.0)
    @inbounds for j in eachindex(beta)
        normsq[column_groups[j]] += metric[j] * abs2(beta[j])
    end
    @inbounds for j in eachindex(beta)
        g = column_groups[j]
        D[j] = metric[j] * group_weights[g] / (normsq[g] + group_sizes[g] * epsilon)
    end
    return D
end

# The reweighted-ridge iteration on a precomputed Gram system, shared by
# `solve_coefficients(::GroupAdaptiveRidge, …)` (`beta0 = nothing`, cold start) and the
# `select_fit` λ-path driver (`beta0` = the previous λ's solution, warm start; the
# caller forms `XtX`/`Xty` once for the whole path). Requires `lambda > 0` — the
# callers route the exact `lambda == 0` case to OLS.
function _solve_gar(XtX::Matrix{Float64}, Xty::Vector{Float64}, lambda::Float64,
                    column_groups::Vector{Int}, group_weights::Vector{Float64},
                    group_sizes::Vector{Int}, metric::Vector{Float64},
                    epsilon::Float64, max_iter::Int, tol::Float64;
                    beta0::Union{Nothing,Vector{Float64}} = nothing)::Vector{Float64}
    p = length(Xty)
    D = Vector{Float64}(undef, p)
    normsq = Vector{Float64}(undef, length(group_weights))
    _check_free_block(XtX, metric, "GroupAdaptiveRidge")
    # Iteration 0 (cold start): the fixed-weight ridge `Dⱼ = mⱼ·v_g` — with a uniform
    # metric and unit weights this is numerically Ridge(lambda), matching the
    # AdaptiveRidge initialization. A warm start replaces it with the caller's `beta0`
    # and enters the loop directly. With lambda > 0 and every penalized Dⱼ > 0
    # (epsilon > 0), the shifted matrix is SPD on the penalized block throughout — do
    # NOT widen the callers' `lambda == 0` routing to a tolerance.
    beta = if beta0 === nothing
        @inbounds for j = 1:p
            D[j] = metric[j] * group_weights[column_groups[j]]
        end
        Symmetric(XtX + lambda * Diagonal(D)) \ Xty
    else
        copy(beta0)
    end
    for _ = 1:max_iter
        _gar_weights!(D, beta, column_groups, group_weights, group_sizes, metric,
                      epsilon, normsq)
        beta_new = Symmetric(XtX + lambda * Diagonal(D)) \ Xty
        rel = _irls_rel_change(beta_new, beta, metric)
        beta = beta_new
        rel < tol && break
    end
    return beta
end

function solve_coefficients(est::GroupAdaptiveRidge, X::AbstractMatrix, y::AbstractVector;
                            groups = nothing)::Vector{Float64}
    length(est.column_groups) == size(X, 2) || throw(DimensionMismatch(
        "GroupAdaptiveRidge column_groups length $(length(est.column_groups)) does not " *
        "match design-matrix column count $(size(X, 2)); the labels were likely built " *
        "on a different SCEBasis."))
    # Exact zero only (see the AdaptiveRidge method above): at lambda = 0 the penalty
    # diagonal is gone and the Gram matrix may be singular — route to QR (min-norm) OLS.
    est.lambda == 0.0 && return solve_coefficients(OLS(), X, y)
    m = _metric_vector(est.metric, size(X, 2), "GroupAdaptiveRidge")
    return _solve_gar(Matrix{Float64}(X' * X), Vector{Float64}(X' * y), est.lambda,
                      est.column_groups, est.group_weights, est.group_sizes, m,
                      est.epsilon, est.max_iter, est.tol)
end

function solve_coefficients(est::PrecomputedPilot, X::AbstractMatrix, y::AbstractVector;
                            groups = nothing)::Vector{Float64}
    length(est.beta) == size(X, 2) || throw(DimensionMismatch(
        "PrecomputedPilot coefficient length $(length(est.beta)) does not match " *
        "design-matrix column count $(size(X, 2)); the pilot was likely fit on a " *
        "different SCEBasis."))
    return copy(est.beta)
end

_estimator_metric(e::Ridge) = e.metric
_estimator_metric(e::AdaptiveRidge) = e.metric
_estimator_metric(e::GroupAdaptiveRidge) = e.metric
_estimator_provenance(e::Ridge) = e.metric_provenance
_estimator_provenance(e::AdaptiveRidge) = e.metric_provenance
_estimator_provenance(e::GroupAdaptiveRidge) = e.metric_provenance

# --- column selection --------------------------------------------------------------

# A column-structured estimator must follow whatever column selection its caller made
# before the solve — the pointed fit's vanishing-column freeze, and `refit`'s chosen
# support — so its per-column metadata has to shrink with the design.
#
# For GroupAdaptiveRidge the reduction preserves every group NORM
# exactly (frozen coefficients are exact zeros), but the group SIZE p_g drops by
# the frozen count, so the weight w_g = v_g/(‖β_g‖² + p_g·ε) moves at O(ε) —
# material only for a group already at the ε floor, and defensible there: ε is
# documented as a per-coefficient floor and a frozen column carries no
# coefficient. Note the deliberate divergence from the energy side, where
# ASR-frozen columns STAY in column_groups and keep their p_g contribution
# (upstream SLCE.jl; this package's energy side has no frozen columns at all).
# Groups emptied by the reduction are relabeled away (the estimator's
# every-label-present contract). Every other estimator passes through unchanged.
#
# A penalty metric makes `Ridge` and `AdaptiveRidge` column-structured too, so they get
# their own methods: a metric is per-column data and has to be cut with the freeze or
# the solve hits a `DimensionMismatch`.
_reduce_to_active(estimator::AbstractEstimator, ::BitVector) = estimator

# The active-column slice of an estimator's penalty metric, or `nothing` when it has
# none. Refuses a metric built on a different column count by name.
function _reduce_metric(metric::Union{Nothing,Vector{Float64}}, active::BitVector,
                        what::AbstractString)::Union{Nothing,Vector{Float64}}
    metric === nothing && return nothing
    length(metric) == length(active) || throw(DimensionMismatch(
        "$what metric length $(length(metric)) does not match the pointed design " *
        "column count $(length(active)); build it on THIS basis (penalty_metric(mb))"))
    return metric[active]
end

function _reduce_to_active(estimator::Ridge, active::BitVector)::Ridge
    (estimator.metric === nothing || all(active)) && return estimator
    return Ridge(estimator.lambda, _reduce_metric(estimator.metric, active, "Ridge"),
                 estimator.metric_provenance)
end

function _reduce_to_active(estimator::AdaptiveRidge, active::BitVector)::AdaptiveRidge
    (estimator.metric === nothing || all(active)) && return estimator
    return AdaptiveRidge(estimator.lambda, estimator.epsilon, estimator.max_iter,
                         estimator.tol,
                         _reduce_metric(estimator.metric, active, "AdaptiveRidge"),
                         estimator.metric_provenance)
end

function _reduce_to_active(estimator::GroupAdaptiveRidge,
                           active::BitVector)::GroupAdaptiveRidge
    all(active) && return estimator
    length(estimator.column_groups) == length(active) || throw(DimensionMismatch(
        "GroupAdaptiveRidge column_groups length $(length(estimator.column_groups)) " *
        "does not match the pointed design column count $(length(active)); build " *
        "the labels on THIS basis (salc_groups(mb))"))
    sub = estimator.column_groups[findall(active)]
    remap = Dict{Int,Int}()
    labels = [get!(remap, g, length(remap) + 1) for g in sub]
    old = Vector{Int}(undef, length(remap))
    for (g, n) in remap
        old[n] = g
    end
    return GroupAdaptiveRidge(labels, estimator.group_weights[old];
                              lambda = estimator.lambda, epsilon = estimator.epsilon,
                              max_iter = estimator.max_iter, tol = estimator.tol,
                              metric = _reduce_metric(estimator.metric, active,
                                                      "GroupAdaptiveRidge"),
                              metric_provenance = estimator.metric_provenance)
end
