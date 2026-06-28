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

L2-penalized least squares with penalty `lambda`.
"""
struct Ridge <: AbstractEstimator
    lambda::Float64
end
Ridge(; lambda::Real = 0.0) = Ridge(Float64(lambda))

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
    islinear(est) -> Bool

Whether `est` is a linear (closed-form) estimator. Gates closed-form diagnostics.
Penalty-path estimators ([`ElasticNet`](@ref) / [`Lasso`](@ref)) are not linear.
"""
islinear(::AbstractEstimator) = false
islinear(::OLS) = true
islinear(::Ridge) = true

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
