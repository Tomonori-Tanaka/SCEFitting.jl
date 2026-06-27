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
    islinear(est) -> Bool

Whether `est` is a linear (closed-form) estimator. Gates closed-form diagnostics.
"""
islinear(::AbstractEstimator) = false
islinear(::OLS) = true
islinear(::Ridge) = true

"""
    solve_coefficients(est, X, y) -> Vector{Float64}

Solve for the SALC coefficients from the (already centered) design matrix `X` and
target `y`. See [`AbstractEstimator`](@ref) for the `(X, y)` contract.
"""
function solve_coefficients(est::AbstractEstimator, X::AbstractMatrix, y::AbstractVector)
    error("solve_coefficients has no method for $(typeof(est)); load the backend " *
          "package (e.g. `using GLMNet` for Lasso/ElasticNet).")
end

function solve_coefficients(::OLS, X::AbstractMatrix, y::AbstractVector)::Vector{Float64}
    return X \ y   # QR-based least squares (more robust than the normal equations)
end

function solve_coefficients(est::Ridge, X::AbstractMatrix, y::AbstractVector)::Vector{Float64}
    n = size(X, 2)
    return Symmetric(X' * X + est.lambda * I(n)) \ (X' * y)
end
