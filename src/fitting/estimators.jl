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
    OLS()

Ordinary least squares via the normal equations.
"""
struct OLS <: AbstractEstimator end

"""
    Ridge(; lambda = 0.0)

L2-penalized least squares with penalty `lambda ≥ 0` (`lambda = 0` is OLS).
"""
struct Ridge <: AbstractEstimator
    lambda::Float64

    function Ridge(lambda::Real)
        (lambda >= 0 && isfinite(lambda)) ||
            throw(ArgumentError("lambda must be finite and ≥ 0; got $lambda"))
        return new(Float64(lambda))
    end
end
Ridge(; lambda::Real = 0.0) = Ridge(lambda)

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
not floored by `epsilon = eps(Float64)` and miscalibrates those weights. Any pilot is
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
    AdaptiveRidge(; lambda, epsilon = 1e-8, max_iter = 50, tol = 1e-6)

Iterative Adaptive Ridge (Frommlet & Nuel 2016). Approximates an L0-penalized fit by
repeatedly refitting a per-coefficient weighted ridge problem

    min_b ‖y − Xb‖² + lambda·Σⱼ wⱼ·bⱼ²

and updating the weights `wⱼ = 1 / (bⱼ² + epsilon)` between iterations. Iteration zero
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
end
function AdaptiveRidge(; lambda::Real, epsilon::Real = 1e-8, max_iter::Integer = 50,
                       tol::Real = 1e-6)
    (lambda >= 0 && isfinite(lambda)) ||
        throw(ArgumentError("lambda must be finite and ≥ 0; got $lambda"))
    epsilon > 0 || throw(ArgumentError("epsilon must be > 0; got $epsilon"))
    max_iter >= 1 || throw(ArgumentError("max_iter must be ≥ 1; got $max_iter"))
    tol > 0 || throw(ArgumentError("tol must be > 0; got $tol"))
    return AdaptiveRidge(Float64(lambda), Float64(epsilon), Int(max_iter), Float64(tol))
end

# Compact display: the default struct printer would dump the whole `beta` vector
# (hundreds–thousands of entries) and recurse into the nested pilot.
Base.show(io::IO, p::PrecomputedPilot) =
    print(io, "PrecomputedPilot(", length(p.beta), " coefficients)")
Base.show(io::IO, e::AdaptiveLasso) =
    print(io, "AdaptiveLasso(pilot=", e.pilot, ", lambda=", e.lambda, ", gamma=",
          e.gamma, ", epsilon=", e.epsilon, ", standardize=", e.standardize, ")")

"""
    islinear(est) -> Bool

Whether `est` is a linear (closed-form) estimator — one whose fitted values are
`ŷ = H·y` with a hat matrix `H` independent of `y`. Gates closed-form diagnostics.
[`OLS`](@ref), [`Ridge`](@ref), and [`AdaptiveRidge`](@ref) (linear in the
converged-weight sense, the standard adaptive-ridge approximation) return `true`;
the penalty-path / pilot estimators ([`ElasticNet`](@ref) / [`Lasso`](@ref) /
[`AdaptiveLasso`](@ref)) are not linear. Extends `StatsAPI.islinear`.
"""
islinear(::AbstractEstimator) = false
islinear(::OLS) = true
islinear(::Ridge) = true
islinear(::AdaptiveRidge) = true

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

function solve_coefficients(::OLS, X::AbstractMatrix, y::AbstractVector;
                            groups = nothing)::Vector{Float64}
    return X \ y   # QR-based least squares (more robust than the normal equations)
end

function solve_coefficients(est::Ridge, X::AbstractMatrix, y::AbstractVector;
                            groups = nothing)::Vector{Float64}
    n = size(X, 2)
    return Symmetric(X' * X + est.lambda * I(n)) \ (X' * y)
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
    # Iteration 0: uniform-weight ridge (numerically Ridge(lambda)). Only the penalty
    # diagonal changes between iterations; the Gram matrix `XtX` is formed once. With
    # lambda > 0 and every wⱼ > 0 (epsilon > 0), the shifted matrix is SPD throughout.
    beta = Symmetric(XtX + est.lambda * I(p)) \ Xty
    w = Vector{Float64}(undef, p)
    for _ = 1:est.max_iter
        @. w = 1.0 / (beta^2 + est.epsilon)
        beta_new = Symmetric(XtX + est.lambda * Diagonal(w)) \ Xty
        # Relative ∞-norm change; the eps floor only guards the all-zero β case.
        rel = mapreduce((a, b) -> abs(a - b), max, beta_new, beta) /
              max(maximum(abs, beta_new), eps(Float64))
        beta = beta_new
        rel < est.tol && break
    end
    return beta
end

function solve_coefficients(est::PrecomputedPilot, X::AbstractMatrix, y::AbstractVector;
                            groups = nothing)::Vector{Float64}
    length(est.beta) == size(X, 2) || throw(DimensionMismatch(
        "PrecomputedPilot coefficient length $(length(est.beta)) does not match " *
        "design-matrix column count $(size(X, 2)); the pilot was likely fit on a " *
        "different SCEBasis."))
    return copy(est.beta)
end
