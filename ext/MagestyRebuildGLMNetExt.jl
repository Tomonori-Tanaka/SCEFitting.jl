module MagestyRebuildGLMNetExt

# Provides `solve_coefficients(::ElasticNet, X, y)` when GLMNet is loaded. The
# `ElasticNet` / `Lasso` types live in the core package (so they can be named and
# dispatched on without this dependency); only the actual GLMNet solve lights up
# here, mirroring the Spglib and Sunny extensions.
#
# Contract (see `AbstractEstimator`): `X` is already column-centered and `y`
# centered, so we fit with `intercept = false` — GLMNet adds no intercept and `j0`
# is recovered analytically downstream in `fit`. GLMNet minimizes
#   (1/2n)·‖y − Xβ‖² + λ·[(1−α)/2·‖β‖₂² + α·‖β‖₁]
# and returns β on the original (un-standardized) scale.

using MagestyRebuild: ElasticNet
import MagestyRebuild: solve_coefficients
using GLMNet: glmnet, glmnetcv

function solve_coefficients(est::ElasticNet, X::AbstractMatrix, y::AbstractVector;
                            groups = nothing)::Vector{Float64}
    Xf = Matrix{Float64}(X)
    yf = Vector{Float64}(y)
    if est.lambda === nothing
        return _solve_cv(est, Xf, yf, groups)
    else
        # A single user-specified penalty: GLMNet wants a λ vector and computes the
        # path warm-started, so a one-element path is a valid (if less warm-started)
        # request for exactly this λ.
        path = glmnet(Xf, yf; alpha = est.alpha, lambda = [est.lambda],
                      standardize = est.standardize, intercept = false)
        return collect(Float64, path.betas[:, 1])
    end
end

# Cross-validated penalty selection over GLMNet's automatic λ path. Folds are over
# resampling units — whole `groups` if given (so a co-fit configuration's energy and
# torque rows stay together), otherwise individual rows.
function _solve_cv(est::ElasticNet, Xf::Matrix{Float64}, yf::Vector{Float64},
                   groups)::Vector{Float64}
    n = length(yf)
    units = groups === nothing ? collect(1:n) : collect(groups)
    nunits = length(unique(units))
    # GLMNet's own default cap; each fold needs a few resampling units to fit.
    nf = min(est.nfolds, div(nunits, 3))
    nf >= 2 || error("ElasticNet cross-validation needs at least 6 resampling units " *
                     "for ≥ 2 folds; got $nunits. Pass an explicit `lambda` to skip CV.")
    cv = glmnetcv(Xf, yf; alpha = est.alpha, folds = _make_folds(units, nf, est.seed),
                  standardize = est.standardize, intercept = false, nlambda = est.nlambda)
    idx = _select_lambda(cv, est.select)
    return collect(Float64, cv.path.betas[:, idx])
end

# Deterministic, balanced, seed-controlled fold assignment (no RNG dependency): rank
# the distinct resampling units by a seeded hash, deal ranks round-robin into `nf`
# folds, then map each row to its unit's fold. Rows sharing a unit label never split
# across the train/holdout boundary. Same seed ⇒ identical folds ⇒ reproducible CV
# (within a Julia session/version, since `hash` is version-dependent).
function _make_folds(units::AbstractVector, nf::Int, seed::Int)::Vector{Int}
    uniq = unique(units)
    order = sortperm([hash((seed, u)) for u in uniq])
    foldof = Dict{eltype(uniq),Int}()
    @inbounds for rank in eachindex(order)
        foldof[uniq[order[rank]]] = mod1(rank, nf)
    end
    return [foldof[u] for u in units]
end

# Index into the descending CV λ path. `:lambda_min` is the lowest-CV-error λ;
# `:lambda_1se` is the most-regularized λ (largest λ ⇒ smallest index ⇒ sparsest)
# whose mean CV loss is still within one standard error of that minimum.
function _select_lambda(cv, select::Symbol)::Int
    imin = argmin(cv.meanloss)
    select === :lambda_min && return imin
    thresh = cv.meanloss[imin] + cv.stdloss[imin]
    for i = 1:imin
        cv.meanloss[i] <= thresh && return i
    end
    return imin
end

end # module MagestyRebuildGLMNetExt
