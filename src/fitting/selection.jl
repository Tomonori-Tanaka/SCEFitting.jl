# Model selection for the fit-accuracy-vs-Monte-Carlo-cost trade-off: group labels and
# a-priori MC costs over an `SCEBasis` (the fixed weights of `GroupAdaptiveRidge`), GCV /
# effective degrees of freedom for linear estimators, and the λ-path driver with the
# cost-aware Pareto selection rule. Lives after `fit.jl`/`diagnostics.jl` in the include
# order because it needs `SCEBasis`, `SCEDataset`, `SCEFit`, and `_assemble_problem`.

"""
    salc_groups(basis::SCEBasis) -> Vector{Int}

Per-design-matrix-column group labels (contiguous `1:G`, one label per SALC in
`SALCKey` order): columns grouped by `(key.body, key.orbit_id, key.ls)`. This is the
granularity at which Monte-Carlo contraction entries vanish — all `Lf` / `block`
channels of one cluster orbit and `l`-multiset share their entry support, so an entry
disappears only when **every** coefficient of the group is zero. Feed the labels to
[`GroupAdaptiveRidge`](@ref) (or use the `GroupAdaptiveRidge(basis; ...)` convenience
constructor, which calls this for you).
"""
function salc_groups(basis::SCEBasis)::Vector{Int}
    ks = basis.salc_basis.keys
    labels = Vector{Int}(undef, length(ks))
    g = 0
    for j in eachindex(ks)
        # keys are sorted by (body, orbit_id, ls, Lf, block), so equal
        # (body, orbit_id, ls) runs are contiguous — label at the change points
        if j == 1 || (ks[j].body, ks[j].orbit_id, ks[j].ls) !=
                     (ks[j-1].body, ks[j-1].orbit_id, ks[j-1].ls)
            g += 1
        end
        labels[j] = g
    end
    return labels
end

# One Monte-Carlo contraction entry: (member sites, l-assignment, tensor index). Two
# SALCs touching an identical entry key have their coefficients folded into a single
# contraction weight downstream, so the entry's cost is paid once per group.
const _EntryKey = Tuple{Vector{Int},Vector{SVector{3,Int}},Vector{Int},Vector{Int}}

# Function barrier over the rank-erased `folded::Array{Float64}` (runtime rank = body
# order): pushes one `_EntryKey` per nonzero tensor element.
function _push_entries!(set::Set{_EntryKey}, atoms::Vector{Int},
                        shifts::Vector{SVector{3,Int}}, ls::Vector{Int},
                        folded::Array{Float64})::Nothing
    for idx in CartesianIndices(folded)
        # exact != 0.0 is deliberate coupling: it mirrors SCEMonteCarlo's own
        # exact-zero entry skip, so the count matches what a sweep would touch —
        # do not soften to an isapprox
        if folded[idx] != 0.0
            push!(set, (atoms, shifts, ls, collect(Tuple(idx))))
        end
    end
    return nothing
end

# Contiguity check shared by `group_costs` (the `GroupAdaptiveRidge` inner constructor
# performs the same validation on its own copy): labels cover 1:G with no gaps.
function _validate_labels(labels::AbstractVector{<:Integer}, n::Int, what::String)
    length(labels) == n ||
        throw(ArgumentError("$what length $(length(labels)) ≠ number of SALCs $n"))
    isempty(labels) && return 0
    minimum(labels) >= 1 ||
        throw(ArgumentError("$what labels must be ≥ 1; got $(minimum(labels))"))
    G = Int(maximum(labels))
    counts = zeros(Int, G)
    for g in labels
        counts[g] += 1
    end
    all(>(0), counts) ||
        throw(ArgumentError("$what labels must cover 1:$G with no gaps; " *
                            "empty group(s): $(findall(==(0), counts))"))
    return G
end

"""
    group_costs(basis::SCEBasis,
                labels::AbstractVector{<:Integer} = salc_groups(basis)) -> Vector{Int}

Per-group Monte-Carlo cost: the number of **distinct contraction entries** — keys
`(member sites, l-assignment, nonzero tensor index)` — in the union over the group's
SALCs (canonical members, schema v4). This is an **a-priori proxy (a lower bound)**
for the entries a `SCEMonteCarlo` sweep realizes: an entry vanishes only when the
whole group is zero, and SALC channels whose tensors are proportional fold into one
realized entry, while non-proportional channels of one group are realized separately
downstream — the union undercounts those. The relative ordering it induces is what
the selection needs. Costs are additive across the [`salc_groups`](@ref) partition
(distinct `(body, orbit_id, ls)` groups never share an entry key); `labels` may also
be any coarser contiguous `1:G` partition of the columns.
"""
function group_costs(basis::SCEBasis,
                     labels::AbstractVector{<:Integer} = salc_groups(basis))::Vector{Int}
    sl = basis.salc_basis.salcs
    G = _validate_labels(labels, length(sl), "group_costs")
    sets = [Set{_EntryKey}() for _ = 1:G]
    for j in eachindex(sl)
        set = sets[labels[j]]
        for m in sl[j].members, t in m.terms
            _push_entries!(set, m.atoms, m.shifts, t.ls, t.folded)
        end
    end
    return length.(sets)
end

"""
    cost_weights(basis::SCEBasis; theta::Real = 1.0)
        -> (; labels::Vector{Int}, weights::Vector{Float64})

Fixed [`GroupAdaptiveRidge`](@ref) weights over the [`salc_groups`](@ref) partition:

    v_g = √p_g · (c_g / c̄)^theta

with `p_g` the group's column count (the Yuan–Lin group-size factor), `c_g =`
[`group_costs`](@ref) and `c̄` their mean. `theta ∈ [0, 1]` sets the cost-vs-accuracy
tilt of the penalty: `theta = 0` is cost-blind group selection, `theta = 1` penalizes
each group in proportion to its Monte-Carlo cost — an expensive group must then earn
its keep with a correspondingly larger error reduction. Sweeping `theta` changes the
*order* in which groups are eliminated along a λ path, so the lower envelope over
several `theta` values traces the (cost, error) Pareto front (see [`select_fit`](@ref)).
"""
function cost_weights(basis::SCEBasis; theta::Real = 1.0)
    (0 <= theta <= 1) || throw(ArgumentError("theta must be in [0, 1]; got $theta"))
    labels = salc_groups(basis)
    c = group_costs(basis, labels)
    isempty(c) && return (; labels, weights = Float64[])
    all(>(0), c) && all(isfinite, c) ||
        error("internal: a canonical SALC group has no nonzero tensor entry")
    G = length(c)
    p = zeros(Int, G)
    for g in labels
        p[g] += 1
    end
    cbar = mean(c)
    t = Float64(theta)
    weights = [sqrt(p[g]) * (c[g] / cbar)^t for g = 1:G]
    return (; labels, weights)
end

"""
    GroupAdaptiveRidge(basis::SCEBasis; lambda, theta = 1.0, epsilon = 1e-8,
                       max_iter = 50, tol = 1e-6)

Cost-weighted group estimator for `basis`: [`salc_groups`](@ref) column labels with the
fixed [`cost_weights`](@ref)`(basis; theta)` weights. See the primary
[`GroupAdaptiveRidge`](@ref) constructor for the estimator itself.
"""
function GroupAdaptiveRidge(basis::SCEBasis; lambda::Real, theta::Real = 1.0,
                            epsilon::Real = 1e-8, max_iter::Integer = 50,
                            tol::Real = 1e-6)
    lw = cost_weights(basis; theta = theta)
    return GroupAdaptiveRidge(lw.labels, lw.weights; lambda = lambda, epsilon = epsilon,
                              max_iter = max_iter, tol = tol)
end

# --- GCV / effective degrees of freedom -------------------------------------------

# The converged quadratic-penalty diagonal `(lambda, w)` of a linear estimator: the
# fitted values are `ŷ = X(X'X + λ·Diagonal(w))⁻¹X'y` with `w` frozen at the fitted
# coefficients. `w === nothing` ⇔ unpenalized (OLS, or λ = 0), where the effective dof
# is the design rank. The adaptive members recompute their converged diagonal from the
# fitted `beta` through the SAME weight formulas the solvers iterate (`AdaptiveRidge`'s
# `1/(β² + ε)`, `_gar_weights!` for the group form) — the formula lives in one place.
_penalty_diagonal(::OLS, beta::Vector{Float64}) = (0.0, nothing)
function _penalty_diagonal(est::Ridge, beta::Vector{Float64})
    est.lambda == 0.0 && return (0.0, nothing)
    return (est.lambda, ones(Float64, length(beta)))
end
function _penalty_diagonal(est::AdaptiveRidge, beta::Vector{Float64})
    est.lambda == 0.0 && return (0.0, nothing)
    return (est.lambda, @.(1.0 / (beta^2 + est.epsilon)))
end
function _penalty_diagonal(est::GroupAdaptiveRidge, beta::Vector{Float64})
    est.lambda == 0.0 && return (0.0, nothing)
    length(beta) == length(est.column_groups) || throw(DimensionMismatch(
        "coefficient length $(length(beta)) ≠ column_groups length " *
        "$(length(est.column_groups))"))
    w = Vector{Float64}(undef, length(beta))
    normsq = Vector{Float64}(undef, length(est.group_weights))
    _gar_weights!(w, beta, est.column_groups, est.group_weights, est.group_sizes,
                  est.epsilon, normsq)
    return (est.lambda, w)
end
_penalty_diagonal(est::AbstractEstimator, beta::Vector{Float64}) =
    throw(ArgumentError("gcv/effective_dof require a linear estimator " *
                        "(`islinear`); got $(typeof(est))"))

# Numerical-rank dof of the unpenalized smoother (LinearAlgebra.rank-style relative
# tolerance on the singular values).
function _rank_df(X::Matrix{Float64})::Float64
    s = svdvals(X)
    (isempty(s) || s[1] == 0.0) && return 0.0
    tolr = maximum(size(X)) * eps(Float64) * s[1]
    return Float64(count(>(tolr), s))
end

# Penalized effective dof `tr(X(X'X + λ·Diagonal(w))⁻¹X') = Σᵢ sᵢ/(sᵢ + λ)` over the
# eigenvalues `s` of the weighted Gram `X̃'X̃` (`p ≤ n`) or its dual `X̃X̃'` (`n < p`),
# `X̃ = X·D^{-1/2}` — always an eigenproblem on the smaller side, never an `n × p` SVD.
# A λ-path caller passes its cached `XtX` so the `p ≤ n` branch touches only `p × p`
# data per λ. Tiny negative eigenvalues from roundoff are clamped out.
function _edof(X::Matrix{Float64}, lambda::Float64, w::Vector{Float64};
               XtX::Union{Nothing,Matrix{Float64}} = nothing)::Float64
    n, p = size(X)
    M = if p <= n
        if XtX === nothing
            Xt = X ./ sqrt.(w)'
            Symmetric(Xt' * Xt)
        else
            isq = 1.0 ./ sqrt.(w)
            Symmetric(isq .* XtX .* isq')          # D^{-1/2}·X'X·D^{-1/2}
        end
    else
        Xt = X ./ sqrt.(w)'
        Symmetric(Xt * Xt')
    end
    df = 0.0
    for s in eigvals(M)
        s > 0.0 || continue
        df += s / (s + lambda)
    end
    return df
end

"""
    effective_dof(f::SCEFit) -> Float64

Effective degrees of freedom of a linear-estimator fit: `tr(H) + 1`, where `H =
X(X'X + λ·Diagonal(w))⁻¹X'` is the hat matrix of the assembled (centered / whitened)
problem with the penalty diagonal frozen at the fitted coefficients, and the `+1`
counts the analytic intercept `j0`. For an unpenalized fit ([`OLS`](@ref), or
`lambda = 0`) this is the design rank `+1`. Distinct from [`dof`](@ref), the raw
parametric count. Linear estimators only ([`islinear`](@ref)); the adaptive members
([`AdaptiveRidge`](@ref) / [`GroupAdaptiveRidge`](@ref)) are handled in the standard
converged-weight sense.
"""
function effective_dof(f::SCEFit)::Float64
    islinear(f.estimator) || throw(ArgumentError(
        "effective_dof requires a linear estimator (`islinear`); " *
        "got $(typeof(f.estimator))"))
    X, _, _, _, _ = _assemble_problem(f.dataset, f.torque_weight)
    lambda, w = _penalty_diagonal(f.estimator, f.jphi)
    df = w === nothing ? _rank_df(X) : _edof(X, lambda, w)
    return df + 1.0
end

# The GCV score (and the effective dof it used, intercept included) on an already-
# assembled problem; shared by `gcv(::SCEFit)` and the `select_fit` λ-path driver
# (which passes its cached `XtX`). `w === nothing` ⇔ unpenalized. The score is `Inf`
# when `df` approaches the row count `n` (the near-interpolating regime, where the GCV
# denominator loses meaning); the ≥ 1 slack keeps the score from exploding on rounding
# when `df ≈ n`.
function _gcv_score(X::Matrix{Float64}, y::Vector{Float64}, beta::Vector{Float64},
                    lambda::Float64, w::Union{Nothing,Vector{Float64}};
                    XtX::Union{Nothing,Matrix{Float64}} = nothing,
                    )::Tuple{Float64,Float64}
    n = size(X, 1)
    df = (w === nothing ? _rank_df(X) : _edof(X, lambda, w; XtX = XtX)) + 1.0
    n - df < max(1.0, 1e-8 * n) && return (Inf, df)
    rss = sum(abs2, y .- X * beta)
    return (n * rss / (n - df)^2, df)
end

"""
    gcv(f::SCEFit) -> Float64

Generalized cross-validation score of a linear-estimator fit:

    GCV = n·RSS / (n − df)²

over the `n` rows of the assembled (centered / whitened) problem, with `df =`
[`effective_dof`](@ref). Returns `Inf` in the near-interpolating regime `df → n`.
Linear estimators only ([`islinear`](@ref)).

!!! warning "Torque co-fits"
    With `torque_weight > 0` the energy row and the torque-component rows of one
    configuration are correlated, but GCV treats all rows as exchangeable — the score
    is then optimistic (the same leak configuration-grouped CV folds avoid). For a
    co-fit, prefer the grouped cross-validation criterion
    ([`select_fit`](@ref)`(...; criterion = :cv)`); use this GCV as a fast reference.
"""
function gcv(f::SCEFit)::Float64
    islinear(f.estimator) || throw(ArgumentError(
        "gcv requires a linear estimator (`islinear`); got $(typeof(f.estimator))"))
    X, y, _, _, _ = _assemble_problem(f.dataset, f.torque_weight)
    lambda, w = _penalty_diagonal(f.estimator, f.jphi)
    return first(_gcv_score(X, y, f.jphi, lambda, w))
end

# --- λ-path driver with the cost-aware Pareto selection rule ----------------------

# Deterministic, balanced, seed-controlled fold assignment (no RNG dependency): rank
# the distinct resampling units by a seeded hash, deal ranks round-robin into `nf`
# folds, then map each row to its unit's fold. Rows sharing a unit label never split
# across the train/holdout boundary. A core port of the GLMNet extension's
# `_make_folds` (the extension cannot be referenced from here); same seed ⇒ identical
# folds within a Julia session/version (`hash` is version-dependent).
function _grouped_folds(units::AbstractVector, nf::Int, seed::Int)::Vector{Int}
    uniq = unique(units)
    order = sortperm([hash((seed, u)) for u in uniq])
    foldof = Dict{eltype(uniq),Int}()
    @inbounds for rank in eachindex(order)
        foldof[uniq[order[rank]]] = mod1(rank, nf)
    end
    return [foldof[u] for u in units]
end

# The default relative alive floor: a group counts as alive at a given λ when one of
# its columns' scaled magnitude |βⱼ|·‖X[:,j]‖ exceeds this fraction of the largest
# scaled magnitude at that λ. The adaptive-ridge alive/dead gap is many orders of
# magnitude, so the exact ratio is uncritical — any value inside the gap partitions
# identically.
const _ALIVE_RTOL = 1e-6

# The cost-aware selection rule: among the path points whose score is within
# `(1 + delta)` of the finite minimum, pick the one of smallest predicted cost; ties
# go to the earlier index (the larger λ, i.e. the more-regularized fit). `Inf` scores
# (near-interpolating fits) are never eligible.
function _select_pareto(scores::Vector{Float64}, costs::Vector{Float64},
                        delta::Float64)::Int
    delta >= 0 || throw(ArgumentError("delta must be ≥ 0; got $delta"))
    length(scores) == length(costs) || throw(DimensionMismatch(
        "scores length $(length(scores)) ≠ costs length $(length(costs))"))
    emin = Inf
    for s in scores
        isfinite(s) && s < emin && (emin = s)
    end
    isfinite(emin) || throw(ArgumentError(
        "no finite score on the λ path (every fit is near-interpolating); " *
        "extend `lambdas` toward larger values"))
    best = 0
    for i in eachindex(scores)
        scores[i] <= (1 + delta) * emin || continue
        if best == 0 || costs[i] < costs[best]
            best = i
        end
    end
    return best
end

"""
    SelectionPath

Result of [`select_fit`](@ref): the descending λ path with the per-λ selection score
(GCV, or grouped-CV mean squared error), effective dof (`NaN` under `criterion = :cv`,
where it is not computed), alive-group count, and predicted Monte-Carlo cost
`Σ_{g alive} c_g`; plus the selection tolerance `delta`, the effective absolute alive
`threshold` at the selected λ, the `selected` index, and the selected `fit` (re-solved
cold at the selected λ, so `fit(SCEFit, dataset, estimator)` reproduces it — and its
row of the table is re-derived from that cold solve, so `fit` / `threshold` /
`n_alive[selected]` / `cost[selected]` are mutually consistent). De-bias with
`refit(path.fit; threshold = path.threshold)`, which reproduces exactly the reported
alive support. A Tables.jl source with one row per λ (columns `lambda`, `score`,
`edof`, `n_alive`, `cost`, `selected`).
"""
struct SelectionPath
    lambda::Vector{Float64}       # descending
    score::Vector{Float64}
    criterion::Symbol             # :gcv | :cv
    edof::Vector{Float64}         # NaN per entry when criterion == :cv
    n_alive::Vector{Int}
    cost::Vector{Float64}
    delta::Float64
    threshold::Float64            # effective absolute alive threshold at `selected`
    selected::Int
    fit::SCEFit
end

Tables.istable(::Type{SelectionPath}) = true
Tables.columnaccess(::Type{SelectionPath}) = true
Tables.columns(p::SelectionPath) =
    (; lambda = p.lambda, score = p.score, edof = p.edof, n_alive = p.n_alive,
       cost = p.cost, selected = [i == p.selected for i in eachindex(p.lambda)])

function Base.show(io::IO, ::MIME"text/plain", p::SelectionPath)
    print(io, "SelectionPath (criterion = :", p.criterion, ", delta = ", p.delta,
          "; ", length(p.lambda), " λ):")
    for i in eachindex(p.lambda)
        print(io, "\n  λ = ", round(p.lambda[i]; sigdigits = 4),
              "  score = ", round(p.score[i]; sigdigits = 5),
              "  n_alive = ", p.n_alive[i],
              "  cost = ", round(p.cost[i]; sigdigits = 5))
        i == p.selected && print(io, "   ← selected")
    end
end

"""
    select_fit(dataset::SCEDataset, est::GroupAdaptiveRidge;
               lambdas, torque_weight = 0.0, criterion = :gcv, delta = 0.05,
               costs = nothing, threshold = nothing, nfolds = 5, seed = 1)
        -> SelectionPath

Fit `est` along the descending λ path `lambdas` (each solve warm-started from the
previous λ's coefficients), score every fit, and select the **cheapest** λ whose score
is within `(1 + delta)` of the path minimum — the cost-aware generalization of the
conventional `:lambda_1se` rule. The returned [`SelectionPath`](@ref) carries the full
per-λ table and the selected fit (re-solved cold, so it is reproducible by a plain
[`fit`](@ref) call); follow with [`refit`](@ref) to de-bias the surviving groups.

`est` supplies the column groups, fixed weights, and IRLS controls; **its own `lambda`
is ignored** (the path is `lambdas`). Scoring:

- `criterion = :gcv` (default) — the [`gcv`](@ref) score from the closed-form hat
  matrix; fast, no refitting. See the co-fit caveat in the [`gcv`](@ref) docstring:
  with `torque_weight > 0` prefer `:cv`.
- `criterion = :cv` — `nfolds`-fold configuration-grouped cross-validation (folds
  never split a configuration's energy/torque rows; deterministic seeded fold
  assignment). The centering/whitening constants stay global — the score ranks λ, it
  is not an unbiased error estimate.

A group is **alive** at a given λ when any of its columns clears the
scaled-magnitude rule `|jϕⱼ|·‖X[:, j]‖ > threshold` on the assembled design (the same
support rule as [`refit`](@ref)). `GroupAdaptiveRidge` crushes dead groups to
tiny-but-**nonzero** values (the `epsilon` floor), so the default
`threshold = nothing` uses a per-λ **relative** floor, `1e-6` of that λ's largest
scaled magnitude — the alive/dead gap spans many orders of magnitude, so any ratio
inside the gap gives the same partition; an absolute number reproduces `refit`'s rule
verbatim (and `threshold = 0` counts every group alive — the cost column is then
flat). The effective absolute threshold at the selected λ is returned as
`path.threshold`; de-bias with `refit(path.fit; threshold = path.threshold)` to
realize exactly the reported support. The predicted Monte-Carlo cost of a fit is
`Σ_{g alive} c_g` with `c_g` from `costs` (default: `SCEFitting.group_costs` of the
dataset's basis under `est`'s column partition). `delta` sets the accuracy tolerance
of the cost–error trade; sweep the `theta` of `SCEFitting.cost_weights` to tilt the
penalty itself and trace a Pareto front over both knobs.
"""
function select_fit(dataset::SCEDataset, est::GroupAdaptiveRidge;
                    lambdas::AbstractVector{<:Real}, torque_weight::Real = 0.0,
                    criterion::Symbol = :gcv, delta::Real = 0.05,
                    costs::Union{Nothing,AbstractVector{<:Real}} = nothing,
                    threshold::Union{Nothing,Real} = nothing, nfolds::Integer = 5,
                    seed::Integer = 1)::SelectionPath
    isempty(dataset.y_E) && throw(ArgumentError("dataset has no observations"))
    w = Float64(torque_weight)
    (0.0 <= w <= 1.0) || throw(ArgumentError("torque_weight must be in [0, 1]; got $w"))
    if w > 0 && !has_torque(dataset)
        throw(ArgumentError("torque_weight = $w but the dataset has no torque data"))
    end
    isempty(lambdas) && throw(ArgumentError("lambdas must be nonempty"))
    all(l -> isfinite(l) && l >= 0, lambdas) ||
        throw(ArgumentError("lambdas must be finite and ≥ 0"))
    criterion in (:gcv, :cv) ||
        throw(ArgumentError("criterion must be :gcv or :cv; got :$criterion"))
    delta >= 0 || throw(ArgumentError("delta must be ≥ 0; got $delta"))
    threshold === nothing || threshold >= 0 ||
        throw(ArgumentError("threshold must be ≥ 0 (or nothing for the relative " *
                            "default); got $threshold"))
    nfolds >= 2 || throw(ArgumentError("nfolds must be ≥ 2; got $nfolds"))
    G = length(est.group_weights)
    length(est.column_groups) == n_salcs(dataset.basis) || throw(DimensionMismatch(
        "estimator column_groups length $(length(est.column_groups)) ≠ basis column " *
        "count $(n_salcs(dataset.basis))"))
    cg = costs === nothing ?
         Float64.(group_costs(dataset.basis, est.column_groups)) : Float64.(costs)
    length(cg) == G ||
        throw(ArgumentError("costs length $(length(cg)) ≠ number of groups $G"))

    lams = sort!(unique(Float64.(lambdas)); rev = true)
    nl = length(lams)
    X, y, _, _, rowgroups = _assemble_problem(dataset, w)
    n = size(X, 1)
    XtX = Matrix{Float64}(X' * X)
    Xty = Vector{Float64}(X' * y)
    colnorms = [norm(view(X, :, j)) for j = 1:size(X, 2)]

    # Warm-started descending path: each IRLS is seeded with the previous (more
    # regularized, already group-sparse) λ's solution.
    betas = Vector{Vector{Float64}}(undef, nl)
    prev = nothing
    for i = 1:nl
        b = lams[i] == 0.0 ? (X \ y) :
            _solve_gar(XtX, Xty, lams[i], est.column_groups, est.group_weights,
                       est.group_sizes, est.epsilon, est.max_iter, est.tol;
                       beta0 = prev)
        betas[i] = b
        prev = b
    end

    # Alive groups (the refit support rule, per group) and the predicted MC cost. The
    # adaptive-ridge iteration never produces exact zeros (the ε floor), so the
    # `threshold === nothing` default applies a per-λ *relative* floor on the scaled
    # magnitudes — the alive/dead gap spans many orders, so any ratio inside it gives
    # the same partition. Returns (alive count, predicted cost, absolute threshold).
    thr_abs = threshold === nothing ? nothing : Float64(threshold)
    function alive_stats!(alive::BitVector, b::Vector{Float64})
        smax = 0.0
        for j in eachindex(b)
            m = abs(b[j]) * colnorms[j]
            m > smax && (smax = m)
        end
        t = thr_abs === nothing ? _ALIVE_RTOL * smax : thr_abs
        fill!(alive, false)
        for j in eachindex(b)
            abs(b[j]) * colnorms[j] > t && (alive[est.column_groups[j]] = true)
        end
        c = sum(cg[g] for g = 1:G if alive[g]; init = 0.0)
        return count(alive), c, t
    end
    n_alive = Vector{Int}(undef, nl)
    cost = Vector{Float64}(undef, nl)
    alive = falses(G)
    for i = 1:nl
        n_alive[i], cost[i], _ = alive_stats!(alive, betas[i])
    end

    edof = fill(NaN, nl)
    score = Vector{Float64}(undef, nl)
    if criterion === :gcv
        wv = Vector{Float64}(undef, length(Xty))
        normsq = Vector{Float64}(undef, G)
        for i = 1:nl
            if lams[i] == 0.0
                score[i], edof[i] = _gcv_score(X, y, betas[i], 0.0, nothing)
            else
                _gar_weights!(wv, betas[i], est.column_groups, est.group_weights,
                              est.group_sizes, est.epsilon, normsq)
                score[i], edof[i] = _gcv_score(X, y, betas[i], lams[i], wv; XtX = XtX)
            end
        end
    else
        units = rowgroups === nothing ? collect(1:n) : rowgroups
        nunits = length(unique(units))
        nf = min(Int(nfolds), div(nunits, 3))
        nf >= 2 || throw(ArgumentError(
            "cross-validation needs at least 6 resampling units for ≥ 2 folds; got " *
            "$nunits. Use criterion = :gcv or pass more data."))
        nf < Int(nfolds) &&
            @warn "select_fit: reducing CV folds so every fold keeps ≥ 3 resampling " *
                  "units" requested = Int(nfolds) effective = nf units = nunits
        folds = _grouped_folds(units, nf, seed)
        sse = zeros(Float64, nl)
        for k = 1:nf
            ho = findall(==(k), folds)
            Xho = X[ho, :]
            yho = y[ho]
            # training Gram by downdating the cached full Gram — no per-fold X pass
            XtX_tr = XtX .- Xho' * Xho
            Xty_tr = Xty .- Xho' * yho
            prevf = nothing
            for i = 1:nl
                bf = if lams[i] == 0.0
                    tr_rows = findall(!=(k), folds)
                    X[tr_rows, :] \ y[tr_rows]
                else
                    _solve_gar(XtX_tr, Xty_tr, lams[i], est.column_groups,
                               est.group_weights, est.group_sizes, est.epsilon,
                               est.max_iter, est.tol; beta0 = prevf)
                end
                prevf = bf
                sse[i] += sum(abs2, yho .- Xho * bf)
            end
        end
        score .= sse ./ n
    end

    sel = _select_pareto(score, cost, Float64(delta))
    est_sel = GroupAdaptiveRidge(lams[sel], est.column_groups, est.group_weights,
                                 est.epsilon, est.max_iter, est.tol)
    fsel = fit(SCEFit, dataset, est_sel; torque_weight = w)
    # Re-derive the selected row from the cold re-solve, so `fit` / `threshold` /
    # `n_alive[selected]` / `cost[selected]` are mutually consistent (warm and cold
    # agree only within the IRLS tol — a knife-edge coefficient could differ).
    n_alive[sel], cost[sel], t_sel = alive_stats!(alive, fsel.jphi)
    return SelectionPath(lams, score, criterion, edof, n_alive, cost, Float64(delta),
                         t_sel, sel, fsel)
end
