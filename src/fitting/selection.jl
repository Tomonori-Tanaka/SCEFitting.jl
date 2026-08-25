# Model selection for the fit-accuracy-vs-Monte-Carlo-cost trade-off: group labels and
# a-priori MC costs over an `SCEBasis` (the fixed weights of `GroupAdaptiveRidge`), GCV /
# effective degrees of freedom for linear estimators, and the λ-path driver with the
# cost-aware Pareto selection rule. Lives after `fit.jl`/`diagnostics.jl` in the include
# order because it needs `SCEBasis`, `SCEDataset`, `SCEFit`, and `_assemble_problem`.

"""
    salc_groups(basis::SCEBasis) -> Vector{Int}

Per-design-matrix-column group labels (contiguous `1:G`, one label per SALC in
`SALCKey` order): columns grouped by `(key.body, key.orbit_id, key.decors)`. This is
the granularity at which Monte-Carlo contraction entries vanish — all `L_S` / `Lf` /
`block` channels of one cluster orbit and decoration multiset share their entry
support, so an entry
disappears only when **every** coefficient of the group is zero. Feed the labels to
[`GroupAdaptiveRidge`](@ref) (or use the `GroupAdaptiveRidge(basis; ...)` convenience
constructor, which calls this for you).
"""
function salc_groups(basis::SCEBasis)::Vector{Int}
    ks = basis.salc_basis.keys
    labels = Vector{Int}(undef, length(ks))
    g = 0
    for j in eachindex(ks)
        # keys are sorted by (body, orbit_id, decors, L_S, Lf, block), so equal
        # (body, orbit_id, decors) runs are contiguous — label at the change points
        if j == 1 || (ks[j].body, ks[j].orbit_id, ks[j].decors) !=
                     (ks[j-1].body, ks[j-1].orbit_id, ks[j-1].decors)
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

Per-group Monte-Carlo **sweep** cost: over the union of the group's **distinct
contraction entries** — keys `(member sites, l-assignment, nonzero tensor index)`
(canonical members, schema v4) — the summed **member-site count** of each entry.
A Metropolis sweep evaluates each entry once per member site position (every site
a member touches carries the entry in its site program), so an N-body entry costs
N site-program slots per sweep, not 1: pricing by the bare entry count is the size
of the *energy* program (walked once per run) and mis-ranked a 3-body group at 2/3
of its real sweep cost relative to a 2-body group at equal entry count.
[Backported from SLCE.jl a596ea3.] This is an **a-priori proxy (a lower bound)**
for what a sweep realizes: an entry vanishes only when the whole group is zero.
The relative ordering it induces is what the selection needs. Costs are additive
across the [`salc_groups`](@ref) partition (distinct `(body, orbit_id, decors)`
groups never share an entry key); `labels` may also be any coarser contiguous `1:G`
partition of the columns.
"""
function group_costs(basis::SCEBasis,
                     labels::AbstractVector{<:Integer} = salc_groups(basis))::Vector{Int}
    sl = basis.salc_basis.salcs
    all(is_pure_spin(k) for k in basis.salc_basis.keys) || throw(ArgumentError(
        "group_costs: the basis carries displacement-decorated sectors; the " *
        "Monte-Carlo entry-key model is pure-spin"))
    G = _validate_labels(labels, length(sl), "group_costs")
    sets = [Set{_EntryKey}() for _ = 1:G]
    for j in eachindex(sl)
        set = sets[labels[j]]
        for m in sl[j].members, t in m.terms
            # Pure-spin entry key (the MC-contract migration moves this to
            # slots); identical values to the v4 per-site ls on today's bases.
            _push_entries!(set, m.atoms, m.shifts, _term_spin_ls(t), t.folded)
        end
    end
    # One site-program slot per member site of each distinct entry (k[3] is the
    # entry's SPIN-axis rank list; on a pure-spin basis every axis is a distinct
    # site, so its length is the member's site count = body order).
    return [sum(k -> length(k[3]), s; init = 0) for s in sets]
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
                       max_iter = 50, tol = 1e-6, torque_weight = 0.0,
                       metric = :basis, metric_nconfig = 2048, metric_seed = 1)

Cost-weighted group estimator for `basis`: [`salc_groups`](@ref) column labels with the
fixed [`cost_weights`](@ref)`(basis; theta)` weights, and by default the basis-intrinsic
[`penalty_metric`](@ref)`(basis; torque_weight, ...)`. See the primary
[`GroupAdaptiveRidge`](@ref) constructor for the estimator itself, and
[`penalty_metric`](@ref) for what the metric does and why it is on by default.

`torque_weight` must be the weight the fit will run at — the metric is a property of
the assembled design, which mixes the two blocks by it, and the fitting doors refuse a
mismatch. Pass `metric = nothing` for the unweighted penalty, or a vector of your own.
"""
function GroupAdaptiveRidge(basis::SCEBasis; lambda::Real, theta::Real = 1.0,
                            epsilon::Real = 1e-8, max_iter::Integer = 50,
                            tol::Real = 1e-6, torque_weight::Real = 0.0,
                            metric = :basis, metric_nconfig::Integer = 2048,
                            metric_seed::Integer = 1)
    lw = cost_weights(basis; theta = theta)
    m, pv = _basis_metric(basis, metric, torque_weight, metric_nconfig, metric_seed)
    return GroupAdaptiveRidge(lw.labels, lw.weights; lambda = lambda, epsilon = epsilon,
                              max_iter = max_iter, tol = tol, metric = m,
                              metric_provenance = pv)
end

"""
    Ridge(basis::SCEBasis; lambda, torque_weight = 0.0, metric = :basis,
          metric_nconfig = 2048, metric_seed = 1)

Ridge for `basis`, carrying the basis-intrinsic
[`penalty_metric`](@ref)`(basis; torque_weight, ...)` so that λ means the same thing
across the three penalized estimators.

`torque_weight` is the only metric keyword without a `metric_` prefix, deliberately:
it must equal the `torque_weight` the fit runs at, because the assembled design mixes
the energy and torque blocks by it and the column scales move with it. The fitting
doors refuse a mismatch. `metric = nothing` gives the unweighted penalty and a vector
is taken as given; anything else is refused by name. Reuse one metric across a λ sweep
with `SCEFitting.with_lambda` rather than rebuilding it per point.
"""
function Ridge(basis::SCEBasis; lambda::Real, torque_weight::Real = 0.0,
               metric = :basis, metric_nconfig::Integer = 2048,
               metric_seed::Integer = 1)
    m, pv = _basis_metric(basis, metric, torque_weight, metric_nconfig, metric_seed)
    return Ridge(lambda, m, pv)
end

"""
    AdaptiveRidge(basis::SCEBasis; lambda, epsilon = 1e-8, max_iter = 50, tol = 1e-6,
                  torque_weight = 0.0, metric = :basis, metric_nconfig = 2048,
                  metric_seed = 1)

The per-coefficient adaptive ridge for `basis`, carrying the basis-intrinsic
[`penalty_metric`](@ref)`(basis; torque_weight, ...)`. Keyword semantics as in
[`Ridge`](@ref)`(basis; ...)`.
"""
function AdaptiveRidge(basis::SCEBasis; lambda::Real, epsilon::Real = 1e-8,
                       max_iter::Integer = 50, tol::Real = 1e-6,
                       torque_weight::Real = 0.0, metric = :basis,
                       metric_nconfig::Integer = 2048, metric_seed::Integer = 1)
    m, pv = _basis_metric(basis, metric, torque_weight, metric_nconfig, metric_seed)
    return AdaptiveRidge(; lambda = lambda, epsilon = epsilon, max_iter = max_iter,
                         tol = tol, metric = m, metric_provenance = pv)
end

# Resolve the `metric` keyword of a basis-aware constructor: `:basis` builds the
# reference metric and stamps its provenance, anything else is taken as given (and
# carries no provenance — the caller owns it).
function _basis_metric(basis::SCEBasis, metric, torque_weight::Real, nconfig::Integer,
                       seed::Integer)
    metric === :basis || return (_checked_metric_keyword(metric), nothing)
    m = penalty_metric(basis; torque_weight = torque_weight, nconfig = nconfig,
                       seed = seed)
    pv = MetricProvenance(:energy, torque_weight, false, nconfig, seed,
                          basis.salc_basis.fingerprint)
    return (m, pv)
end

"""
    penalty_metric(basis::SCEBasis; torque_weight = 0.0, nconfig = 2048, seed = 1)
        -> Vector{Float64}

The per-column penalty scale of `basis`: one entry per SALC column, in design-column
(`SALCKey`) order, for the `metric` field of [`Ridge`](@ref) / [`AdaptiveRidge`](@ref) /
[`GroupAdaptiveRidge`](@ref).

`λ·Σⱼβⱼ²` is not invariant under rescaling a design column, and SALC column norms are
set by basis conventions — an orbit's member count, and the ordering multiplicity the
member fold absorbs — rather than by physics. A larger column norm means a smaller
fitted coefficient at the same physical effect, hence *less* shrinkage: the plain
penalty carries an accidental prior in favour of large orbits and high body order.
Weighting the penalty by `mⱼ` removes it, leaving whatever prior the caller states
deliberately (the `theta` of [`cost_weights`](@ref)).

The scale is the reference norm of the column **as the estimator sees it** — the
assembled, centered / whitened design of `_assemble_problem` at this `torque_weight`:

    mⱼ(w) = (1 − w)·Var[Φⱼ] + w·(1 / 3n_atoms)·E[ Σ_a ‖(∂Φⱼ/∂e_a) × e_a‖² ]

with both moments taken over `nconfig` independent uniform-random spin configurations.
The torque term's `1/(3·n_atoms)` is the per-row average the assembly's `√(w/n_T)`
already applies (`n_T = n_E·3·n_atoms`); dropping it would misscale the two blocks
against each other by the atom count.

The metric is a property of the **basis**, not of the training data — deliberately, so
that λ can be compared between cells, so that cross-validation needs no per-fold
recomputation to stay leak-free, and so the prior sits on the function space rather
than on how strongly a particular training set happened to excite each column. The
price is that a column the reference ensemble excites weakly but the training data
drives hard is effectively under-penalized; `nothing` (uniform) remains available, and
the ratio `mⱼ / Var_train[Φⱼ]` is worth looking at when the training set is far from
uniform (near-collinear low-temperature configurations, say).

`nconfig` / `seed` control the reference ensemble. The generator is specified inside
this package rather than taken from `Random`, so the sequence does not move between
Julia versions; the estimate converges as `1/√nconfig` to a closed-form expectation.

The default is sized from that convergence, not guessed. Measured on bcc Fe 2×2×2
(`lmax = 2`, 2- and 3-body columns) as the spread of `mⱼ` over eight independent
seeds: the relative standard error is **3.3 % median / 5.3 % worst column at
`nconfig = 2048`**, falling to 1.6 % / 2.8 % at 8192. Body order barely moves it
(2-body and 3-body columns agree within the spread), so `1/√nconfig` from these
numbers sizes any basis. That residual is a seed-dependent wobble on the *prior*, an order of magnitude
smaller than the systematic factor the metric removes — orbit size times the ordering
multiplicity of the member fold, which spans decades — but it is not zero, so build the
metric ONCE and reuse it across a λ sweep (`SCEFitting.with_lambda`) rather than
rebuilding per point with a different seed.

See also [`Ridge`](@ref)`(basis; ...)` and [`GroupAdaptiveRidge`](@ref)`(basis; ...)`,
which attach the metric and its provenance for you.
"""
function penalty_metric(basis::SCEBasis; torque_weight::Real = 0.0,
                        nconfig::Integer = 2048, seed::Integer = 1)::Vector{Float64}
    w = Float64(torque_weight)
    (isfinite(w) && 0 <= w <= 1) ||
        throw(ArgumentError("torque_weight must be in [0, 1]; got $torque_weight"))
    nat = n_atoms(basis.crystal)
    cfgs = _reference_configs(nat, Int(nconfig), Int(seed))
    sal = basis.salc_basis.salcs
    p = length(sal)
    K = length(cfgs)
    m = Vector{Float64}(undef, p)
    # PRECONDITION on the torque block: the assembly whitens it by `√(w/n_T)` with
    # `n_T = length(dataset.y_T)`, and this restates that as `n_T = n_E·3·n_atoms` —
    # true when every training configuration carries torques, which is the only shape
    # `SCEDataset` builds today. On a dataset with torques for only SOME configurations
    # the two differ, and the energy/torque mixture the metric encodes stops matching
    # the assembled column norms. `MetricProvenance` records `torque_weight` but not
    # coverage, so the fitting doors cannot see it: this is the price of the metric
    # being a property of the BASIS rather than of a dataset.
    # Columns are independent and each task owns one, so the result is identical at
    # any thread count. The torque block is accumulated, never materialized: the full
    # `K·3·n_atoms × p` design would be hundreds of MB for a supercell basis.
    # `:greedy` because column cost grows steeply with body order and the columns are
    # in sorted-key order, so the expensive ones are contiguous at the end — the
    # default (contiguous per-thread chunks) would put them all in one task.
    #
    # The three `w` regimes are separate loops rather than one loop with a test: at
    # `w = 1` the energy term is multiplied by zero, and evaluating it anyway costs
    # roughly the whole `w = 0` column.
    Threads.@threads :greedy for j = 1:p
        scratch = SALCScratch()
        s1 = 0.0
        s2 = 0.0
        st = 0.0
        if w < 1.0
            @inbounds for c in cfgs
                phi = evaluate_salc(sal[j], c, scratch)
                s1 += phi
                s2 += phi * phi
            end
        end
        if w > 0.0
            G = Matrix{Float64}(undef, 3, nat)
            @inbounds for c in cfgs
                fill!(G, 0.0)
                accumulate_grad!(G, sal[j], c, 1.0, scratch)
                for a = 1:nat
                    ea = SVector{3,Float64}(c[1, a], c[2, a], c[3, a])
                    ga = SVector{3,Float64}(G[1, a], G[2, a], G[3, a])
                    st += sum(abs2, cross(ga, ea))
                end
            end
        end
        # The textbook one-pass variance. It can go slightly negative on a column whose
        # variance is genuinely zero, which the clamp absorbs; the reference ensemble is
        # centered enough (`E[Φ] = 0` for every all-`l ≥ 1` label) that the
        # cancellation this form is known for does not bite here.
        varj = max(0.0, s2 / K - (s1 / K)^2)
        m[j] = (1 - w) * varj + w * (st / K) / (3 * nat)
    end
    _refuse_zero_metric(m, "penalty_metric(::SCEBasis)")
    return m
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
    return (est.lambda, _metric_vector(est.metric, length(beta), "Ridge"))
end
function _penalty_diagonal(est::AdaptiveRidge, beta::Vector{Float64})
    est.lambda == 0.0 && return (0.0, nothing)
    m = _metric_vector(est.metric, length(beta), "AdaptiveRidge")
    return (est.lambda, @.(m / (m * beta^2 + est.epsilon)))
end
function _penalty_diagonal(est::GroupAdaptiveRidge, beta::Vector{Float64})
    est.lambda == 0.0 && return (0.0, nothing)
    length(beta) == length(est.column_groups) || throw(DimensionMismatch(
        "coefficient length $(length(beta)) ≠ column_groups length " *
        "$(length(est.column_groups))"))
    m = _metric_vector(est.metric, length(beta), "GroupAdaptiveRidge")
    D = Vector{Float64}(undef, length(beta))
    normsq = Vector{Float64}(undef, length(est.group_weights))
    _gar_weights!(D, beta, est.column_groups, est.group_weights, est.group_sizes, m,
                  est.epsilon, normsq)
    return (est.lambda, D)
end
_penalty_diagonal(est::AbstractEstimator, beta::Vector{Float64}) =
    throw(ArgumentError("gcv/effective_dof require a linear estimator " *
                        "(`islinear`); got $(typeof(est))"))

# Numerical-rank dof of the unpenalized smoother. The tolerance follows
# `LinearAlgebra.rank`'s default (`minimum(size(X))·eps·s₁`) — `maximum(size(X))`
# inflated the cut by n/p on tall designs and could under-report the rank of a
# borderline-conditioned design. [Backported from SLCE.jl 896180e.]
function _rank_df(X::Matrix{Float64})::Float64
    s = svdvals(X)
    (isempty(s) || s[1] == 0.0) && return 0.0
    tolr = minimum(size(X)) * eps(Float64) * s[1]
    return Float64(count(>(tolr), s))
end

# Penalized effective dof `tr(X(X'X + λ·Diagonal(w))⁻¹X') = Σᵢ sᵢ/(sᵢ + λ)` over the
# eigenvalues `s` of the weighted Gram `X̃'X̃` (`p ≤ n`) or its dual `X̃X̃'` (`n < p`),
# `X̃ = X·D^{-1/2}` — always an eigenproblem on the smaller side, never an `n × p` SVD.
# A λ-path caller passes its cached `XtX` so the `p ≤ n` branch touches only `p × p`
# data per λ. Tiny negative eigenvalues from roundoff are clamped out.
#
# `w[j] == 0` marks column `j` as UNPENALIZED. `X̃` is then undefined — and the cached
# `XtX` form divides by zero without so much as an `Inf` in the trace — so the split is
# taken FIRST, before the `XtX` keyword is consulted.
function _edof(X::Matrix{Float64}, lambda::Float64, w::Vector{Float64};
               XtX::Union{Nothing,Matrix{Float64}} = nothing,
               columns::Union{Nothing,Vector{Int}} = nothing)::Float64
    any(iszero, w) && return _edof_free(X, lambda, w; columns = columns)
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

# `_edof` when part of the penalty diagonal is exactly zero. Split the design as
# `X = [X_F X_P]` (unpenalized / penalized), `D = diag(0, W)`, `A = X'X + λD`:
#
#   tr(H) = tr(A⁻¹X'X) = p − λ·tr((A⁻¹)_PP W),   (A⁻¹)_PP = (X_P'M X_P + λW)⁻¹
#
# with `M = I − X_F(X_F'X_F)⁻¹X_F'` the projector off the unpenalized columns, so
#
#   tr(H) = p_F + Σᵢ sᵢ/(sᵢ + λ),   sᵢ = eig(W^{-1/2}(X_P'M X_P)W^{-1/2}),
#
# recovering the penalized form at `p_F = 0` and giving `df → p_F` as `λ → ∞` (an
# unpenalized column always costs its full degree of freedom). Preconditions: `λ > 0`,
# and `X_F` of full column rank — `v'Av = ‖Xv‖² + λΣ_P w_j v_j²` vanishes only for
# `v_P = 0` and `X_F v_F = 0`, so `A ≻ 0` is exactly a rank condition on `X_F` and the
# penalized block is safe for any `W ≻ 0`. A rank-deficient `X_F` is refused by name:
# the model is not identified and any finite dof reported for it would be meaningless.
#
# `M` is never formed `n × n`. With `Q` the thin-QR basis of `X_F`, the matrix
# `B = X̃_P − Q(Q'X̃_P)` has the same nonzero singular values as `M X_P W^{-1/2}`, and
# the eigenproblem is taken on the smaller of `B'B` (`p_P × p_P`) and `BB'` (`n × n`).
# A cached `XtX` is deliberately unused: this branch runs once per diagnostic, not once
# per point of a λ path.
function _edof_free(X::Matrix{Float64}, lambda::Float64, w::Vector{Float64};
                    columns::Union{Nothing,Vector{Int}} = nothing)::Float64
    free = findall(iszero, w)
    pen = findall(!iszero, w)
    XF = X[:, free]
    # ONE factorization: the SVD supplies both the rank test and the orthonormal basis
    # of the projector, and doing those with two different factorizations would leave
    # the test and the projection able to disagree. The tolerance is `_rank_df`'s.
    F = svd(XF)
    tolr = isempty(F.S) ? 0.0 : minimum(size(XF)) * eps(Float64) * F.S[1]
    p_F = count(>(tolr), F.S)
    if p_F != length(free)
        named = columns === nothing ? free : columns[free]
        throw(ArgumentError(
            "effective dof with unpenalized columns: the unpenalized block must have " *
            "full column rank, but columns $named have numerical rank $p_F < " *
            "$(length(free)). The penalized least-squares problem is then singular " *
            "and its hat matrix undefined — drop the dependent columns, or penalize " *
            "them, rather than reporting a finite dof for a model that is not " *
            "identified"))
    end
    isempty(pen) && return Float64(p_F)
    Xp = X[:, pen] ./ sqrt.(w[pen])'
    Q = view(F.U, :, 1:p_F)
    B = Xp .- Q * (Q' * Xp)
    G = length(pen) <= size(X, 1) ? Symmetric(B' * B) : Symmetric(B * B')
    df = Float64(p_F)
    for s in eigvals(G)
        s > 0.0 || continue
        df += s / (s + lambda)
    end
    return df
end

# The informative row count: at `torque_weight == 1` the energy rows enter the
# assembled design with weight exactly zero, so they carry no information and must
# not be counted by GCV. The zero-weight condition is the SAME expression
# `_assemble_problem` scales with (`√((1 − w)/n_E)`).
_gcv_neff(f::SCEFit)::Int =
    (f.torque_weight == 1.0 ? 0 : size(f.dataset.X_E, 1)) +
    (f.torque_weight > 0 ? length(f.dataset.y_T) : 0)

# `effective_dof`/`gcv` reconstruct the FULL design from `dataset`; a `refit`
# solved on a column support, so the reconstruction describes a model the refit
# did not return — measured upstream: df 40.0 where ≈ 5 was honest, and `gcv` a
# silent `Inf` through the `n − df` guard. Refuse by name rather than answer
# wrong. [Backported from SLCE.jl 54457ca, review M3.]
function _refuse_refit_diagnostic(f::SCEFit, what::AbstractString)
    f.support === nothing || throw(ArgumentError(
        "$what on a `refit` result: the diagnostic reconstructs the full " *
        "$(length(f.jphi))-column design, not the $(length(f.support))-column " *
        "support the refit solved on, so its value would describe a model that " *
        "was not returned. Score the REGULARIZED fit before de-biasing (that is " *
        "what selects λ), and validate the refit by holdout or `cross_validate`."))
    return nothing
end

"""
    effective_dof(f::SCEFit) -> Float64

Effective degrees of freedom of a linear-estimator fit: `tr(H) + 1`, where `H =
X(X'X + λ·Diagonal(w))⁻¹X'` is the hat matrix of the assembled (centered / whitened)
problem with the penalty diagonal frozen at the fitted coefficients, and the `+1`
counts the analytic intercept `j0` (charged only while the energy block carries
weight — at `torque_weight == 1` the `j0` estimate rides rows GCV does not count).
For an unpenalized fit ([`OLS`](@ref), or `lambda = 0`) this is the design rank
`+1`. Distinct from [`dof`](@ref), the raw parametric count. Linear estimators
only ([`islinear`](@ref)).

!!! warning "Adaptive estimators: a lower bound"
    For the adaptive members ([`AdaptiveRidge`](@ref) /
    [`GroupAdaptiveRidge`](@ref)) the penalty diagonal is **frozen** at the
    fitted coefficients (the standard converged-weight treatment). That ignores
    the selection the data-dependent weights perform, so the value is a **lower
    bound** on the true `tr(∂ŷ/∂y)` — measured upstream up to 14 % low at small
    λ against a numerical trace — and a [`gcv`](@ref) built on it is
    correspondingly **optimistic**, most where the weights select hardest.
    Grouped cross-validation ([`select_fit`](@ref)`(...; criterion = :cv)`) is
    the honest criterion.

A [`refit`](@ref) result is **refused**: the diagnostics reconstruct the full
design, not the support the refit solved on, so the value would describe a model
that was not returned. Score the regularized fit before de-biasing.
"""
function effective_dof(f::SCEFit)::Float64
    islinear(f.estimator) || throw(ArgumentError(
        "effective_dof requires a linear estimator (`islinear`); " *
        "got $(typeof(f.estimator))"))
    _refuse_refit_diagnostic(f, "effective_dof")
    X, _, _, _, _ = _assemble_problem(f.dataset, f.torque_weight)
    lambda, w = _penalty_diagonal(f.estimator, f.jphi)
    df = w === nothing ? _rank_df(X) : _edof(X, lambda, w)
    # The intercept is charged only when the energy block carries weight (same
    # zero test as `_gcv_neff` / the assembly's scale expression).
    return df + (f.torque_weight == 1.0 ? 0.0 : 1.0)
end

# The GCV score (and the effective dof it used, intercept included) on an already-
# assembled problem; shared by `gcv(::SCEFit)` and the `select_fit` λ-path driver
# (which passes its cached `XtX`). `w === nothing` ⇔ unpenalized. `n_eff` is the
# informative row count (`size(X, 1)` minus zero-weight rows). The score is `Inf`
# when `df` approaches `n_eff` (the near-interpolating regime, where the GCV
# denominator loses meaning); the ≥ 1 slack keeps the score from exploding on
# rounding when `df ≈ n_eff`.
function _gcv_score(X::Matrix{Float64}, y::Vector{Float64}, beta::Vector{Float64},
                    lambda::Float64, w::Union{Nothing,Vector{Float64}};
                    XtX::Union{Nothing,Matrix{Float64}} = nothing,
                    n_eff::Int = size(X, 1),
                    intercept::Float64 = 1.0,
                    )::Tuple{Float64,Float64}
    n = n_eff
    # `intercept` is 0.0 when the energy block carries zero weight (`w == 1`):
    # j0 is then estimated from rows `n_eff` does not count, so charging it to the
    # counted rows inflated the score by ((n−df)/(n−df−1))² — ~3 % at n_eff = 72
    # and unbounded as df → n_eff. [Backported from SLCE.jl 54457ca, review M4.]
    df = (w === nothing ? _rank_df(X) : _edof(X, lambda, w; XtX = XtX)) + intercept
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

See the [`effective_dof`](@ref) warnings, which apply here unchanged: for the
adaptive estimators the frozen-weight `df` makes the score **optimistic**, and a
[`refit`](@ref) result is **refused** (the reconstruction is not the problem the
refit solved).
"""
function gcv(f::SCEFit)::Float64
    islinear(f.estimator) || throw(ArgumentError(
        "gcv requires a linear estimator (`islinear`); got $(typeof(f.estimator))"))
    _refuse_refit_diagnostic(f, "gcv")
    X, y, _, _, _ = _assemble_problem(f.dataset, f.torque_weight)
    lambda, w = _penalty_diagonal(f.estimator, f.jphi)
    return first(_gcv_score(X, y, f.jphi, lambda, w;
                            n_eff = _gcv_neff(f),
                            intercept = f.torque_weight == 1.0 ? 0.0 : 1.0))
end

# --- λ-path driver with the cost-aware Pareto selection rule ----------------------

# Deterministic, balanced, seed-controlled fold assignment (no RNG dependency): rank
# the distinct resampling units by a seeded hash, deal ranks round-robin into `nf`
# folds, then map each row to its unit's fold. Rows sharing a unit label never split
# across the train/holdout boundary. A core port of the GLMNet extension's
# `_make_folds` (the extension cannot be referenced from here); same seed ⇒ identical
# folds within a Julia session/version (`hash` is version-dependent).
# The fold count every grouped cross-validation in the package settles on: at most the
# requested `nfolds`, and never so many that a fold would hold fewer than three
# resampling units. Fewer than two folds is refused; a reduction warns, because a
# silently different fold count makes two runs incomparable. The `unit` noun is the
# caller's, since "configuration" and "resampling unit" mean the same thing here only
# for an energy-only dataset.
function _cv_fold_count(nunits::Int, nfolds::Integer, caller::AbstractString,
                        unit::AbstractString, extra::AbstractString = "")::Int
    nf = min(Int(nfolds), div(nunits, 3))
    nf >= 2 || throw(ArgumentError(
        "cross-validation needs at least 6 $unit for ≥ 2 folds; got $nunits." * extra))
    nf < Int(nfolds) &&
        @warn "$caller: reducing CV folds so every fold keeps ≥ 3 $unit" requested =
              Int(nfolds) effective = nf units = nunits
    return nf
end

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
  matrix; fast, no refitting. Two caveats from the [`gcv`](@ref) docstring apply
  here with force: with `torque_weight > 0` prefer `:cv` (row-correlation leak),
  and the frozen-weight df of the adaptive penalty makes the score **optimistic
  exactly where it matters for this driver** — most at small λ, where the weights
  select hardest — so `:gcv` leans toward keeping extra groups alive and the
  Pareto rule then buys real Monte-Carlo cost against a biased score. `:cv` is the
  honest criterion; use `:gcv` for fast scans.
- `criterion = :cv` — `nfolds`-fold configuration-grouped cross-validation (folds
  never split a configuration's energy/torque rows; deterministic seeded fold
  assignment). The centering/whitening constants stay global — the score ranks λ, it
  is not an unbiased error estimate. The reported score is the fit's own objective
  `(1 − w)·MSE_E + w·MSE_T` on the pooled out-of-fold residuals — the scale
  [`cross_validate`](@ref)'s `pooled_score` and [`select_support`](@ref)'s `score`
  report, so the three are directly comparable.

!!! warning "The two criteria are not on a common scale"
    `:cv` reports the objective itself, while `:gcv` reports `n·RSS/(n − df)²` over
    the assembled rows — and those rows are already row-scaled by `√((1−w)/n_E)`, so
    the `:gcv` number is smaller by roughly the informative row count. Each score
    ranks its own path (both are invariant under a positive uniform factor); a `:gcv`
    score and a `:cv` score must never be compared to each other, or to a `delta`
    calibrated on the other.

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
    _check_metric_provenance(est, :energy, dataset.basis.salc_basis.fingerprint, w)
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
    # The penalty metric, resolved once for the whole path: every solve, every GCV
    # weight, the per-fold solves, and the cold re-solve of the selected point must see
    # the SAME diagonal, or the returned fit and the score attached to it would
    # describe different estimators.
    metric = _metric_vector(est.metric, size(X, 2), "GroupAdaptiveRidge")

    # Warm-started descending path: each IRLS is seeded with the previous (more
    # regularized, already group-sparse) λ's solution.
    # Non-convergent IRLS solves are counted, not warned about one by one: this driver
    # runs one solve per λ and another per (fold, λ), and a warning repeated over that
    # grid buries the one thing the user can act on — how much of the path is affected.
    nonconv_path = Ref(0)
    nonconv_fold = Ref(0)
    nf_total = Ref(0)
    betas = Vector{Vector{Float64}}(undef, nl)
    prev = nothing
    for i = 1:nl
        b = lams[i] == 0.0 ? (X \ y) :
            _solve_gar(XtX, Xty, lams[i], est.column_groups, est.group_weights,
                       est.group_sizes, metric, est.epsilon, est.max_iter, est.tol;
                       beta0 = prev, nonconvergent = nonconv_path)
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
    # The informative row count and intercept charge follow `_gcv_neff` /
    # `effective_dof`: at `w == 1` the energy rows enter with weight exactly zero,
    # so GCV must not count them, and `j0` is estimated from rows it does not
    # count. [Backported from SLCE.jl 54457ca, review M4.]
    neff = n - (w == 1.0 ? length(dataset.y_E) : 0)
    icpt = w == 1.0 ? 0.0 : 1.0
    if criterion === :gcv
        wv = Vector{Float64}(undef, length(Xty))
        normsq = Vector{Float64}(undef, G)
        for i = 1:nl
            if lams[i] == 0.0
                score[i], edof[i] = _gcv_score(X, y, betas[i], 0.0, nothing;
                                               n_eff = neff, intercept = icpt)
            else
                _gar_weights!(wv, betas[i], est.column_groups, est.group_weights,
                              est.group_sizes, metric, est.epsilon, normsq)
                score[i], edof[i] = _gcv_score(X, y, betas[i], lams[i], wv;
                                               XtX = XtX, n_eff = neff,
                                               intercept = icpt)
            end
        end
    else
        units = rowgroups === nothing ? collect(1:n) : rowgroups
        nunits = length(unique(units))
        nf = _cv_fold_count(nunits, nfolds, "select_fit", "resampling units",
                            " Use criterion = :gcv or pass more data.")
        folds = _grouped_folds(units, nf, seed)
        nf_total[] = nf * nl
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
                               est.group_weights, est.group_sizes, metric,
                               est.epsilon, est.max_iter, est.tol; beta0 = prevf,
                               nonconvergent = nonconv_fold)
                end
                prevf = bf
                sse[i] += sum(abs2, yho .- Xho * bf)
            end
        end
        # NO further division. `_assemble_problem` already row-scaled the design by
        # `√((1−w)/n_E)` and `√(w/n_T)`, so a fold's squared holdout residual is
        # ALREADY per-row; summed over folds (each row held out exactly once) `sse`
        # is exactly `(1−w)·SSE_E/n_E + w·SSE_T/n_T`, which is the objective
        # `cross_validate.pooled_score` and `select_support.score` report. Dividing
        # by `neff` on top of that scaled the score by a further `1/neff` — measured
        # on 60 energy configurations, `select_fit(:cv).score` came out 61.9× below
        # the `cross_validate` pooled score it claims to match. `_select_pareto` is
        # invariant under a positive uniform factor, so only the REPORTED number was
        # wrong; the selection was not.
        score .= sse
    end

    # One report for the whole grid. A score at a non-convergent point describes a
    # smoother that was never solved (`gcv` / `effective_dof` rebuild the penalty
    # diagonal from the returned coefficients), so the selection itself is suspect
    # there — which is why the count, not the first offender, is what gets reported.
    if nonconv_path[] + nonconv_fold[] > 0
        @warn "select_fit: $(nonconv_path[]) of $nl λ points" *
              (nf_total[] > 0 ?
               " and $(nonconv_fold[]) of $(nf_total[]) fold solves" : "") *
              " ended on max_iter = $(est.max_iter) rather than on the reweighting's " *
              "stopping rule. Their coefficients are not fixed points of the " *
              "reweighting, so the score attached to them — and therefore the " *
              "selected λ, if it is one of them — describes a smoother that was " *
              "never solved. Raise `max_iter`, loosen `tol`, or raise `epsilon`."
    end

    sel = _select_pareto(score, cost, Float64(delta))
    est_sel = GroupAdaptiveRidge(lams[sel], est.column_groups, est.group_weights,
                                 est.epsilon, est.max_iter, est.tol, est.metric,
                                 est.metric_provenance)
    fsel = fit(SCEFit, dataset, est_sel; torque_weight = w)
    # Re-derive the selected row from the cold re-solve, so `fit` / `threshold` /
    # `n_alive[selected]` / `cost[selected]` are mutually consistent (warm and cold
    # agree only within the IRLS tol — a knife-edge coefficient could differ). Under
    # :gcv the score/edof are functions of the returned fit, so re-derive them too
    # (`path.score[selected] == gcv(path.fit)` exactly); the :cv score aggregates
    # fold models and does not depend on the full-data solve — leave it.
    n_alive[sel], cost[sel], t_sel = alive_stats!(alive, fsel.jphi)
    if criterion === :gcv
        if lams[sel] == 0.0
            score[sel], edof[sel] = _gcv_score(X, y, fsel.jphi, 0.0, nothing;
                                               n_eff = neff, intercept = icpt)
        else
            wv = Vector{Float64}(undef, length(Xty))
            normsq = Vector{Float64}(undef, G)
            _gar_weights!(wv, fsel.jphi, est.column_groups, est.group_weights,
                          est.group_sizes, metric, est.epsilon, normsq)
            score[sel], edof[sel] = _gcv_score(X, y, fsel.jphi, lams[sel], wv;
                                               XtX = XtX, n_eff = neff,
                                               intercept = icpt)
        end
    end
    # The cold re-derivation above mutated `score[sel]`/`n_alive[sel]`/`cost[sel]`
    # AFTER `_select_pareto` ran on the warm table. Warm and cold agree within the
    # IRLS tol, so the selection is expected to stand — but the docstring's rule
    # ("the cheapest λ within delta of the minimum") must hold on the table the
    # caller sees, so re-check and say so if a knife-edge case ever moves it.
    # [Backported from SLCE.jl 54457ca.]
    sel2 = _select_pareto(score, cost, Float64(delta))
    sel2 == sel ||
        @warn "select_fit: the cold re-derivation of the selected row moved the " *
              "Pareto choice — the returned table's rule-based pick is index " *
              "$sel2 (λ = $(lams[sel2])) while the returned fit is index $sel " *
              "(λ = $(lams[sel])). The two solves differ only within the IRLS " *
              "tolerance; tighten `tol` or thin the λ grid near the tie." maxlog = 1
    return SelectionPath(lams, score, criterion, edof, n_alive, cost, Float64(delta),
                         t_sel, sel, fsel)
end

# --- threshold sweep: the (cost, error) front of de-biased refits -----------------

# The auto threshold grid for `select_support`: at most `n` log-rank-spaced points on
# the sorted per-group scaled magnitudes, each threshold the midpoint between
# consecutive *distinct* magnitudes (a midpoint at an exact tie would equal the tied
# value and, under the strict `>` rule, exclude the whole tie — tied groups die
# together instead), plus `0.0` for the full-support anchor. Returned descending
# (sparsest refit first); duplicates and degenerate ranks collapse, so fewer than `n`
# points can come back (always ≥ 1: the anchor).
function _support_thresholds(n::Integer, m_g::Vector{Float64})::Vector{Float64}
    n >= 2 || throw(ArgumentError("npoints must be ≥ 2; got $n"))
    # An empty group vector would reach `log(0)` below and die with a range error
    # naming nothing the caller holds. [Backported from SLCE.jl 54457ca.]
    isempty(m_g) && throw(ArgumentError(
        "the fit has no coefficient groups to threshold (zero-SALC basis)"))
    ms = sort(m_g; rev = true)
    G = length(ms)
    thr = Float64[]
    for r in unique(round.(Int, exp.(range(log(1), log(G); length = n))))
        if r == G
            push!(thr, 0.0)
        elseif ms[r] > ms[r+1]
            push!(thr, (ms[r] + ms[r+1]) / 2)
        end
    end
    return sort!(unique(thr); rev = true)
end

"""
    SupportPath

Result of [`select_support`](@ref): the descending threshold sweep with, per point,
the alive-group count, predicted Monte-Carlo cost `Σ_{g alive} c_g`, the selection
`score` (the fit's own `(1−w)·MSE_E + w·MSE_T` objective on the evaluation dataset),
and the energy / torque RMSEs (`rmse_torque` is `NaN` when the evaluation dataset
carries no torque data); plus the tolerance `delta`, the `selected` index, and the
de-biased refit `fit` at the selected threshold. A Tables.jl source with one row per
threshold (columns `threshold`, `n_alive`, `cost`, `score`, `rmse_energy`,
`rmse_torque`, `selected`).
"""
struct SupportPath
    threshold::Vector{Float64}    # descending (sparsest first)
    n_alive::Vector{Int}
    cost::Vector{Float64}
    score::Vector{Float64}
    rmse_energy::Vector{Float64}
    rmse_torque::Vector{Float64}  # NaN per entry when the evalset has no torque
    delta::Float64
    selected::Int
    fit::SCEFit
end

Tables.istable(::Type{SupportPath}) = true
Tables.columnaccess(::Type{SupportPath}) = true
Tables.columns(p::SupportPath) =
    (; threshold = p.threshold, n_alive = p.n_alive, cost = p.cost, score = p.score,
       rmse_energy = p.rmse_energy, rmse_torque = p.rmse_torque,
       selected = [i == p.selected for i in eachindex(p.threshold)])

function Base.show(io::IO, ::MIME"text/plain", p::SupportPath)
    print(io, "SupportPath (delta = ", p.delta, "; ", length(p.threshold),
          " thresholds):")
    for i in eachindex(p.threshold)
        print(io, "\n  thr = ", round(p.threshold[i]; sigdigits = 4),
              "  score = ", round(p.score[i]; sigdigits = 5),
              "  n_alive = ", p.n_alive[i],
              "  cost = ", round(p.cost[i]; sigdigits = 5))
        i == p.selected && print(io, "   ← selected")
    end
end

"""
    select_support(f::SCEFit; npoints = 25, thresholds = nothing, delta = 0.05,
                   labels = nothing, costs = nothing, evalset = f.dataset,
                   estimator = OLS())
        -> SupportPath

Trace the (predicted Monte-Carlo cost, error) front of **de-biased refits** of `f`
over a sweep of alive thresholds, and select the cheapest point whose score is within
`(1 + delta)` of the front minimum — the same Pareto rule as [`select_fit`](@ref).
This is the second knob of the cost-aware workflow: [`select_fit`](@ref) sweeps the
penalty λ, while this sweeps the support directly. On real data the group-magnitude
spectrum is typically continuous (no clean alive/dead gap), so most of the
cost–error trade lives here.

Per threshold `t`: the alive groups are those with any column satisfying
`|jϕⱼ|·‖X[:, j]‖ > t` on `f`'s assembled design (the [`refit`](@ref) rule, grouped),
the predicted cost is `Σ_{g alive} c_g`, and the point's fit is
`refit(f, estimator; threshold = t)` — note the refit keeps *columns* above `t`, so a
weak column of an alive group may still be dropped (the group's cost is paid either
way). The `score` is the fit's own objective `(1 − w)·MSE_E + w·MSE_T`
(`w = f.torque_weight`) evaluated on `evalset` — pass a held-out `SCEDataset` (built
on the same basis; see dataset slicing) for an honest error axis; the default is the
in-sample training set.

The sweep is either automatic or explicit, and the two are **separate keywords** on
purpose: `npoints` is a point count for the automatic grid (**at most** that many
points, log-rank-spaced on the per-group magnitude spectrum plus the full-support
anchor — duplicate ranks and exact magnitude ties collapse), while `thresholds` is an
explicit vector of absolute thresholds and overrides it. One keyword carrying both
meanings turned `thresholds = 10` — "sweep down to a magnitude of 10" — into a
silent ten-point grid, distinguishable from `thresholds = [10.0]` only by the
literal's type. `labels`/`costs` default to `SCEFitting.salc_groups` /
`SCEFitting.group_costs` of the training basis. [The npoints/thresholds keyword
split is backported from SLCE.jl 2259c54 — the one behavior change of that naming
batch; the Greek keyword renames stay out.]
"""
function select_support(f::SCEFit;
                        npoints::Integer = 25,
                        thresholds::Union{Nothing,AbstractVector{<:Real}} = nothing,
                        delta::Real = 0.05,
                        labels::Union{Nothing,AbstractVector{<:Integer}} = nothing,
                        costs::Union{Nothing,AbstractVector{<:Real}} = nothing,
                        evalset::SCEDataset = f.dataset,
                        estimator::AbstractEstimator = OLS())::SupportPath
    delta >= 0 || throw(ArgumentError("delta must be ≥ 0; got $delta"))
    _reject_precomputed_pilot(estimator)
    isempty(evalset.y_E) && throw(ArgumentError("evalset has no observations"))
    basis = f.dataset.basis
    lab = labels === nothing ? salc_groups(basis) : Vector{Int}(labels)
    G = _validate_labels(lab, length(f.jphi), "select_support")
    cg = costs === nothing ? Float64.(group_costs(basis, lab)) : Float64.(costs)
    length(cg) == G ||
        throw(ArgumentError("costs length $(length(cg)) ≠ number of groups $G"))
    evalset.basis.salc_basis.fingerprint == basis.salc_basis.fingerprint ||
        throw(ArgumentError("evalset was built on a different SCEBasis than the fit"))
    w = f.torque_weight
    if w > 0 && !has_torque(evalset)
        throw(ArgumentError("f is a torque co-fit (torque_weight = $w) but evalset " *
                            "has no torque data"))
    end

    # per-group max scaled magnitude on the assembled training design (refit's rule)
    X, _, _, _, _ = _assemble_problem(f.dataset, w)
    m_g = zeros(Float64, G)
    for j in eachindex(f.jphi)
        m = abs(f.jphi[j]) * norm(view(X, :, j))
        g = lab[j]
        m > m_g[g] && (m_g[g] = m)
    end
    thr = if thresholds === nothing
        _support_thresholds(npoints, m_g)
    else
        isempty(thresholds) && throw(ArgumentError("thresholds must be nonempty"))
        all(t -> isfinite(t) && t >= 0, thresholds) ||
            throw(ArgumentError("thresholds must be finite and ≥ 0"))
        sort!(unique(Float64.(thresholds)); rev = true)
    end

    nt = length(thr)
    n_alive = Vector{Int}(undef, nt)
    cost = Vector{Float64}(undef, nt)
    score = Vector{Float64}(undef, nt)
    rmseE = Vector{Float64}(undef, nt)
    rmseT = fill(NaN, nt)
    fits = Vector{SCEFit}(undef, nt)
    for i = 1:nt
        alive = falses(G)
        for g = 1:G
            m_g[g] > thr[i] && (alive[g] = true)
        end
        n_alive[i] = count(alive)
        cost[i] = sum(cg[g] for g = 1:G if alive[g]; init = 0.0)
        fr = refit(f, estimator; threshold = thr[i])
        fits[i] = fr
        mseE = mean(abs2, evalset.y_E .- (fr.j0 .+ evalset.X_E * fr.jphi))
        rmseE[i] = sqrt(mseE)
        if has_torque(evalset)
            mseT = mean(abs2, evalset.y_T .- evalset.X_T * fr.jphi)
            rmseT[i] = sqrt(mseT)
            score[i] = (1 - w) * mseE + w * mseT
        else
            score[i] = mseE
        end
    end
    sel = _select_pareto(score, cost, Float64(delta))
    return SupportPath(thr, n_alive, cost, score, rmseE, rmseT, Float64(delta), sel,
                       fits[sel])
end

# --- configuration-grouped K-fold cross-validation (generic) -----------------------

"""
    CVResult

Result of [`cross_validate`](@ref): per-fold holdout scores plus the pooled
out-of-fold error. `score` is the fit's own objective `(1 − w)·MSE_E + w·MSE_T`
(`w = torque_weight`) on each fold's held-out configurations; `rmse_energy` and
`rmse_torque` report the two error axes separately whenever the dataset carries the
data, independent of `torque_weight` (`rmse_torque` is `NaN` per entry for an
energy-only dataset). The `pooled_*` fields aggregate the out-of-fold residuals of
all folds — every configuration is held out exactly once, so they are single
whole-dataset numbers, not means of the per-fold columns. A Tables.jl source with
one row per fold (columns `fold`, `n_holdout`, `score`, `rmse_energy`,
`rmse_torque`).
"""
struct CVResult
    nfolds::Int                   # effective fold count (see cross_validate)
    seed::Int
    torque_weight::Float64
    n_holdout::Vector{Int}        # held-out configurations per fold
    score::Vector{Float64}
    rmse_energy::Vector{Float64}
    rmse_torque::Vector{Float64}  # NaN per entry for an energy-only dataset
    pooled_score::Float64
    pooled_rmse_energy::Float64
    pooled_rmse_torque::Float64   # NaN for an energy-only dataset
end

Tables.istable(::Type{CVResult}) = true
Tables.columnaccess(::Type{CVResult}) = true
Tables.columns(r::CVResult) =
    (; fold = collect(eachindex(r.score)), n_holdout = r.n_holdout, score = r.score,
       rmse_energy = r.rmse_energy, rmse_torque = r.rmse_torque)

function Base.show(io::IO, ::MIME"text/plain", r::CVResult)
    print(io, "CVResult (", r.nfolds, " folds, seed = ", r.seed,
          ", torque_weight = ", r.torque_weight, "):")
    for k in eachindex(r.score)
        print(io, "\n  fold ", k, ": n = ", r.n_holdout[k],
              "  score = ", round(r.score[k]; sigdigits = 5),
              "  rmse_E = ", round(r.rmse_energy[k]; sigdigits = 5))
        isnan(r.rmse_torque[k]) ||
            print(io, "  rmse_T = ", round(r.rmse_torque[k]; sigdigits = 5))
    end
    print(io, "\n  pooled: score = ", round(r.pooled_score; sigdigits = 5),
          "  rmse_E = ", round(r.pooled_rmse_energy; sigdigits = 5))
    isnan(r.pooled_rmse_torque) ||
        print(io, "  rmse_T = ", round(r.pooled_rmse_torque; sigdigits = 5))
end

"""
    cross_validate(dataset::SCEDataset, estimator::AbstractEstimator;
                   torque_weight = 0.0, nfolds = 5, seed = 1) -> CVResult

Configuration-grouped `nfolds`-fold cross-validation of
`fit(SCEFit, dataset, estimator; torque_weight)`: each fold's model is fit from
scratch on the training configurations only (energy centering and torque whitening
included — nothing leaks across the split), then scored on the held-out
configurations in prediction space. A configuration is one resampling unit, so its
energy row and its `3·n_atoms` torque rows always land on the same side of the
split. Fold assignment is deterministic in `seed` (seeded-hash round-robin, the
same rule as [`select_fit`](@ref)'s `:cv`); when the dataset has fewer than
`3·nfolds` configurations the fold count is reduced (with a warning) so every fold
keeps ≥ 3 configurations, and fewer than 6 configurations is an error.

Use it to compare estimators, `torque_weight` settings, or a [`refit`](@ref)-style
support on an equal footing: `rmse_energy` and `rmse_torque` report both error axes
regardless of `torque_weight` (a `torque_weight = 0` fit still gets its torque
error measured when the dataset carries torque data). This differs from
`select_fit(...; criterion = :cv)`, which ranks a λ path under **global**
centering/whitening for speed — `cross_validate` is the honest
generalization-error estimate.

A `PrecomputedPilot` (or an `AdaptiveLasso` carrying one) is rejected: its fixed,
full-data coefficient vector does not depend on the training fold, so the holdout
score would leak the held-out data.
"""
function cross_validate(dataset::SCEDataset, estimator::AbstractEstimator;
                        torque_weight::Real = 0.0, nfolds::Integer = 5,
                        seed::Integer = 1)::CVResult
    w = Float64(torque_weight)
    (0.0 <= w <= 1.0) || throw(ArgumentError("torque_weight must be in [0, 1]; got $w"))
    if w > 0 && !has_torque(dataset)
        throw(ArgumentError("torque_weight = $w but the dataset has no torque data"))
    end
    _check_metric_provenance(estimator, :energy,
                             dataset.basis.salc_basis.fingerprint, w)
    nfolds >= 2 || throw(ArgumentError("nfolds must be ≥ 2; got $nfolds"))
    if _carries_precomputed_pilot(estimator)
        throw(ArgumentError("cross_validate does not accept a PrecomputedPilot (or " *
            "an AdaptiveLasso carrying one): its fixed full-data coefficient vector " *
            "does not depend on the training fold, so the holdout score would leak. " *
            "Pass the estimator that produced the pilot instead."))
    end
    nc = length(dataset)
    nf = _cv_fold_count(nc, nfolds, "cross_validate", "configurations")

    folds = _grouped_folds(collect(1:nc), nf, seed)
    hastq = has_torque(dataset)
    n_holdout = Vector{Int}(undef, nf)
    score = Vector{Float64}(undef, nf)
    rmseE = Vector{Float64}(undef, nf)
    rmseT = fill(NaN, nf)
    sseE = 0.0
    sseT = 0.0
    nE = 0
    nT = 0
    for k = 1:nf
        ho = findall(==(k), folds)
        tr = findall(!=(k), folds)
        f = fit(SCEFit, dataset[tr], estimator; torque_weight = w)
        hset = dataset[ho]
        residE = hset.y_E .- (f.j0 .+ hset.X_E * f.jphi)
        mseE = mean(abs2, residE)
        n_holdout[k] = length(ho)
        rmseE[k] = sqrt(mseE)
        sseE += sum(abs2, residE)
        nE += length(residE)
        if hastq
            residT = hset.y_T .- hset.X_T * f.jphi
            mseT = mean(abs2, residT)
            rmseT[k] = sqrt(mseT)
            sseT += sum(abs2, residT)
            nT += length(residT)
            score[k] = (1 - w) * mseE + w * mseT
        else
            score[k] = mseE
        end
    end
    pE = sseE / nE
    if hastq
        pT = sseT / nT
        return CVResult(nf, Int(seed), w, n_holdout, score, rmseE, rmseT,
                        (1 - w) * pE + w * pT, sqrt(pE), sqrt(pT))
    end
    return CVResult(nf, Int(seed), w, n_holdout, score, rmseE, rmseT, pE, sqrt(pE), NaN)
end
