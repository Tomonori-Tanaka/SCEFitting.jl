"""
    Interaction(; nbody, pair_cutoff, lmax, isotropy = false)

SCE interaction specification: maximum body order, the 2-body cutoff (Å),
per-species maximum `l`, and whether to keep only the isotropic (`Lf = 0`) channel.

`pair_cutoff` may be `Inf`, meaning "every resolvable pair" — the whole
Wigner–Seitz cell of the (super)cell — under the default [`MinimumImage`](@ref)
selection (cf. Magesty's `-1` sentinel).
"""
struct Interaction
    nbody::Int
    pair_cutoff::Float64
    lmax::Vector{Int}
    isotropy::Bool
end
function Interaction(; nbody::Integer, pair_cutoff::Real, lmax::AbstractVector{<:Integer},
                     isotropy::Bool = false)
    nbody >= 1 || throw(ArgumentError("nbody must be ≥ 1; got $nbody"))
    (pair_cutoff > 0 && !isnan(pair_cutoff)) ||
        throw(ArgumentError("pair_cutoff must be positive (or Inf for the full " *
                            "Wigner–Seitz cell); got $pair_cutoff"))
    return Interaction(Int(nbody), Float64(pair_cutoff), collect(Int, lmax), isotropy)
end

"""
    SCEBasis(crystal, interaction; backend = NoSymmetry(), tol = 1e-5,
             images = MinimumImage())

Build the SCE basis for `crystal`: analyze symmetry, enumerate cluster orbits, and
construct the symmetry-adapted SALC basis. Pass `backend = SpglibBackend()` (with
`using Spglib`) for real space-group symmetry.

`images` selects how periodic images are admitted (see [`AbstractImageSelection`](@ref)):
[`MinimumImage`](@ref) (the default) keeps only the minimum-image, resolvable pairs of
a plain-PBC supercell; [`AllImages`](@ref) keeps every image within the cutoff (the
generalized-Bloch / spin-spiral seam, finite cutoff only). The `images` value is also
applied to the cluster-edge admissibility, so the neighbor list and the clusters stay
consistent.
"""
struct SCEBasis
    crystal::Crystal
    spacegroup::SpaceGroup
    salcs::SALCBasis
    interaction::Interaction
end

function SCEBasis(crystal::Crystal, interaction::Interaction;
                 backend::AbstractSymmetryBackend = NoSymmetry(), tol::Real = 1e-5,
                 images::AbstractImageSelection = MinimumImage())::SCEBasis
    sg = analyze_symmetry(backend, crystal; tol = tol)
    nl = build_neighbor_list(crystal, interaction.pair_cutoff, images)
    clusters = build_clusters(crystal, nl, sg; nbody = interaction.nbody, selection = images)
    salcs = build_salc_basis(crystal, sg, clusters;
                             lmax_by_species = interaction.lmax, isotropy = interaction.isotropy)
    return SCEBasis(crystal, sg, salcs, interaction)
end

nsalc(b::SCEBasis) = length(b.salcs)

"""
    SCEDataset(basis, configs, energies)
    SCEDataset(basis, configs, energies, torques)

Pair an [`SCEBasis`](@ref) with training data: spin configurations `configs`
(each `3 × n_atoms`, unit columns) and their `energies`. Materializes the energy
design matrix `X_E[config, salc] = evaluate(salc, config)`.

The four-argument form additionally takes per-configuration torques (each
`3 × n_atoms`, the DFT torque `τ_a = e_a × ∂E/∂e_a` on every atom) and builds the
torque design matrix `X_T` for an energy+torque co-fit (see [`fit`](@ref)). Its
rows are flattened config-major, then atom-major, then `xyz`:
`row = 3·n_atoms·(config−1) + 3·(atom−1) + component`. Torque-free datasets leave
`X_T`/`y_T` empty.
"""
struct SCEDataset
    basis::SCEBasis
    configs::Vector{Matrix{Float64}}
    X_E::Matrix{Float64}
    y_E::Vector{Float64}
    X_T::Matrix{Float64}
    y_T::Vector{Float64}
end

"""
    has_torque(dataset) -> Bool

Whether `dataset` carries torque training data (so an energy+torque co-fit is
possible).
"""
has_torque(d::SCEDataset) = !isempty(d.y_T)

function SCEDataset(basis::SCEBasis, configs::AbstractVector, energies::AbstractVector)::SCEDataset
    cfgs = [Matrix{Float64}(c) for c in configs]
    X = _design_energy(basis, cfgs)
    empty_T = Matrix{Float64}(undef, 0, size(X, 2))
    return SCEDataset(basis, cfgs, X, collect(Float64, energies), empty_T, Float64[])
end

function SCEDataset(basis::SCEBasis, configs::AbstractVector, energies::AbstractVector,
                   torques::AbstractVector)::SCEDataset
    cfgs = [Matrix{Float64}(c) for c in configs]
    length(torques) == length(cfgs) ||
        throw(ArgumentError("got $(length(torques)) torque blocks for $(length(cfgs)) configs"))
    X_E = _design_energy(basis, cfgs)
    X_T = _design_torque(basis, cfgs)
    y_T = _flatten_torques(torques, cfgs)
    return SCEDataset(basis, cfgs, X_E, collect(Float64, energies), X_T, y_T)
end

# Energy design matrix: X_E[config, salc] = Φ_salc(config).
function _design_energy(basis::SCEBasis, cfgs::Vector{Matrix{Float64}})::Matrix{Float64}
    salcs = basis.salcs.salcs
    n = length(cfgs)
    m = length(salcs)
    X = Matrix{Float64}(undef, n, m)
    @inbounds for j = 1:m, i = 1:n   # column-major: stride-1 writes down each column
        X[i, j] = evaluate(salcs[j], cfgs[i])
    end
    return X
end

# Torque design matrix: for each SALC column, each config, each atom, the three
# components of τ_a = e_a × ∂Φ/∂e_a stacked config-major / atom-major / xyz.
function _design_torque(basis::SCEBasis, cfgs::Vector{Matrix{Float64}})::Matrix{Float64}
    salcs = basis.salcs.salcs
    n = length(cfgs)
    m = length(salcs)
    nat = isempty(cfgs) ? 0 : size(cfgs[1], 2)
    block = 3 * nat
    X = Matrix{Float64}(undef, n * block, m)
    G = Matrix{Float64}(undef, 3, nat)
    @inbounds for j = 1:m            # column per SALC (column-major writes)
        for ci = 1:n
            c = cfgs[ci]
            fill!(G, 0.0)
            accumulate_grad!(G, salcs[j], c, 1.0)
            row_off = block * (ci - 1)
            for a = 1:nat
                ea = SVector{3,Float64}(c[1, a], c[2, a], c[3, a])
                ga = SVector{3,Float64}(G[1, a], G[2, a], G[3, a])
                t = cross(ea, ga)
                rb = row_off + 3 * (a - 1)
                X[rb + 1, j] = t[1]
                X[rb + 2, j] = t[2]
                X[rb + 3, j] = t[3]
            end
        end
    end
    return X
end

# Flatten per-config torque targets in the same row order as `_design_torque`.
function _flatten_torques(torques::AbstractVector, cfgs::Vector{Matrix{Float64}})::Vector{Float64}
    nat = isempty(cfgs) ? 0 : size(cfgs[1], 2)
    block = 3 * nat
    y = Vector{Float64}(undef, length(cfgs) * block)
    @inbounds for ci in eachindex(torques)
        τ = torques[ci]
        size(τ) == (3, nat) ||
            throw(ArgumentError("torque block $ci must be 3 × $nat (got $(size(τ)))"))
        rb = block * (ci - 1)
        for a = 1:nat, d = 1:3
            y[rb + 3 * (a - 1) + d] = τ[d, a]
        end
    end
    return y
end

"""
    SCEModel

A lightweight predictor: the basis plus the fitted reference energy `j0` and SALC
coefficients `jphi`. `keys` records each coefficient's [`SALCKey`](@ref).

`SCEModel` is a session-scoped object: `jphi[k]` pairs with `basis.salcs.salcs[k]`
positionally (they are built together by `SCEModel(::SCEFit)`), and `keys` is
diagnostic metadata. A future persistence/reload workflow should re-pair `jphi` to
a freshly built basis **by key**, not by position.
"""
struct SCEModel
    basis::SCEBasis
    j0::Float64
    jphi::Vector{Float64}
    keys::Vector{SALCKey}
end

"""
    SCEFit

A fitted SCE model: the dataset, fitted `j0`/`jphi`, the estimator, the
`torque_weight` used, and the (energy) residuals.
"""
struct SCEFit
    dataset::SCEDataset
    j0::Float64
    jphi::Vector{Float64}
    estimator::AbstractEstimator
    residuals::Vector{Float64}
    torque_weight::Float64
end

"""
    fit(SCEFit, dataset, estimator; torque_weight = 0.0) -> SCEFit

Fit the SCE coefficients. The energy design matrix is column-centered, so the
reference energy `j0` is recovered analytically as `mean(y_E − X_E·jϕ)`
(independent of the estimator), and the centered problem is handed to
`solve_coefficients`.

With `torque_weight = w ∈ (0, 1]` and a torque-carrying `dataset` (see
[`SCEDataset`](@ref)), the fit minimizes the per-sample-normalized objective

    L = (1 − w)·MSE_energy + w·MSE_torque

by whitening: the centered energy block is row-scaled by `√((1−w)/n_E)` and the
torque block by `√(w/n_T)`, then stacked. `j0` does not enter the torque block, so
it is still recovered from the energy data alone. `w = 0` (the default) is a
pure-energy fit; requesting `w > 0` without torque data is an error.
"""
function fit(::Type{SCEFit}, dataset::SCEDataset, estimator::AbstractEstimator;
             torque_weight::Real = 0.0)::SCEFit
    X_E = dataset.X_E
    y_E = dataset.y_E
    isempty(y_E) && throw(ArgumentError("dataset has no observations"))
    w = Float64(torque_weight)
    (0.0 <= w <= 1.0) || throw(ArgumentError("torque_weight must be in [0, 1]; got $w"))
    if w > 0 && !has_torque(dataset)
        throw(ArgumentError("torque_weight = $w but the dataset has no torque data; " *
                            "build it with SCEDataset(basis, configs, energies, torques)"))
    end

    xbar = vec(mean(X_E; dims = 1))
    ybar = mean(y_E)
    if w > 0
        n_E = length(y_E)
        n_T = length(dataset.y_T)
        se = sqrt((1 - w) / n_E)
        sm = sqrt(w / n_T)
        X = vcat((X_E .- xbar') .* se, dataset.X_T .* sm)
        y = vcat((y_E .- ybar) .* se, dataset.y_T .* sm)
    else
        X = X_E .- xbar'
        y = y_E .- ybar
    end
    jphi = solve_coefficients(estimator, X, y)
    j0 = ybar - dot(xbar, jphi)
    residuals = y_E .- (j0 .+ X_E * jphi)
    return SCEFit(dataset, j0, jphi, estimator, residuals, w)
end

"""
    SCEModel(f::SCEFit) -> SCEModel

Extract the lightweight predictor from a fit.
"""
SCEModel(f::SCEFit) = SCEModel(f.dataset.basis, f.j0, f.jphi, f.dataset.basis.salcs.keys)

"""
    predict_energy(model, data) -> Float64 or Vector{Float64}

Predict the energy of a spin configuration (`3 × n_atoms` matrix) or a vector of
configurations.
"""
function predict_energy(model::SCEModel, config::AbstractMatrix{<:Real})::Float64
    salcs = model.basis.salcs.salcs
    e = model.j0
    @inbounds for k in eachindex(model.jphi)
        e += model.jphi[k] * evaluate(salcs[k], config)
    end
    return e
end
predict_energy(model::SCEModel, configs::AbstractVector)::Vector{Float64} =
    [predict_energy(model, c) for c in configs]
predict_energy(f::SCEFit, data) = predict_energy(SCEModel(f), data)

"""
    predict_torque(model, config) -> Matrix{Float64}

Predict the per-atom torque `τ_a = e_a × ∂E/∂e_a` of a spin configuration
(`3 × n_atoms`), returned as a `3 × n_atoms` matrix. A vector of configurations
returns a vector of such matrices. The torque is the analytic derivative of the
same energy surface [`predict_energy`](@ref) evaluates, so the two are consistent
by construction (`τ = e × ∇E`).
"""
function predict_torque(model::SCEModel, config::AbstractMatrix{<:Real})::Matrix{Float64}
    salcs = model.basis.salcs.salcs
    nat = size(config, 2)
    G = zeros(Float64, 3, nat)
    @inbounds for k in eachindex(model.jphi)
        accumulate_grad!(G, salcs[k], config, model.jphi[k])
    end
    T = Matrix{Float64}(undef, 3, nat)
    @inbounds for a = 1:nat
        ea = SVector{3,Float64}(config[1, a], config[2, a], config[3, a])
        ga = SVector{3,Float64}(G[1, a], G[2, a], G[3, a])
        t = cross(ea, ga)
        T[1, a] = t[1]
        T[2, a] = t[2]
        T[3, a] = t[3]
    end
    return T
end
predict_torque(model::SCEModel, configs::AbstractVector)::Vector{Matrix{Float64}} =
    [predict_torque(model, c) for c in configs]
predict_torque(f::SCEFit, data) = predict_torque(SCEModel(f), data)

"""
    coef(f) / intercept(f) / nobs(f) -> coefficients / j0 / number of observations
"""
coef(f::SCEFit) = f.jphi
coef(m::SCEModel) = m.jphi
intercept(f::SCEFit) = f.j0
intercept(m::SCEModel) = m.j0
nobs(f::SCEFit) = length(f.dataset.y_E)

"""
    r2_energy(f) / rmse_energy(f) -> in-sample energy R² / RMSE
"""
function r2_energy(f::SCEFit)::Float64
    y = f.dataset.y_E
    ss_res = sum(abs2, f.residuals)
    ss_tot = sum(abs2, y .- mean(y))
    return ss_tot == 0 ? 1.0 : 1 - ss_res / ss_tot
end
rmse_energy(f::SCEFit)::Float64 = sqrt(mean(abs2, f.residuals))

"""
    r2_torque(f) / rmse_torque(f) -> in-sample torque R² / RMSE

Computed over the flattened per-atom torque components. Requires a
torque-carrying dataset. The torque prediction `τ = X_T·jϕ` has no intercept
(`j0` does not enter it), so the R² baseline is the uncentered total sum of
squares `Σ τ²` (the error of the zero predictor), not the mean-centered one.
"""
function r2_torque(f::SCEFit)::Float64
    has_torque(f.dataset) || throw(ArgumentError("dataset has no torque data"))
    y = f.dataset.y_T
    res = y .- f.dataset.X_T * f.jphi
    ss_res = sum(abs2, res)
    ss_tot = sum(abs2, y)            # no-intercept model ⇒ uncentered baseline
    return ss_tot == 0 ? 1.0 : 1 - ss_res / ss_tot
end
function rmse_torque(f::SCEFit)::Float64
    has_torque(f.dataset) || throw(ArgumentError("dataset has no torque data"))
    return sqrt(mean(abs2, f.dataset.y_T .- f.dataset.X_T * f.jphi))
end
