"""
    Interaction(; nbody, pair_cutoff, lmax, isotropy = false)

SCE interaction specification: maximum body order, the 2-body cutoff (Å),
per-species maximum `l`, and whether to keep only the isotropic (`Lf = 0`) channel.
"""
struct Interaction
    nbody::Int
    pair_cutoff::Float64
    lmax::Vector{Int}
    isotropy::Bool
end
Interaction(; nbody::Integer, pair_cutoff::Real, lmax::AbstractVector{<:Integer},
            isotropy::Bool = false) =
    Interaction(Int(nbody), Float64(pair_cutoff), collect(Int, lmax), isotropy)

"""
    SCEBasis(crystal, interaction; backend = NoSymmetry(), tol = 1e-5)

Build the SCE basis for `crystal`: analyze symmetry, enumerate cluster orbits, and
construct the symmetry-adapted SALC basis. Pass `backend = SpglibBackend()` (with
`using Spglib`) for real space-group symmetry.
"""
struct SCEBasis
    crystal::Crystal
    spacegroup::SpaceGroup
    salcs::SALCBasis
    interaction::Interaction
end

function SCEBasis(crystal::Crystal, interaction::Interaction;
                 backend::AbstractSymmetryBackend = NoSymmetry(), tol::Real = 1e-5)::SCEBasis
    sg = analyze_symmetry(backend, crystal; tol = tol)
    nl = build_neighbor_list(crystal, interaction.pair_cutoff)
    clusters = build_clusters(crystal, nl, sg; nbody = interaction.nbody)
    salcs = build_salc_basis(crystal, sg, clusters;
                             lmax_by_species = interaction.lmax, isotropy = interaction.isotropy)
    return SCEBasis(crystal, sg, salcs, interaction)
end

nsalc(b::SCEBasis) = length(b.salcs)

"""
    SCEDataset(basis, configs, energies)

Pair an [`SCEBasis`](@ref) with training data: spin configurations `configs`
(each `3 × n_atoms`, unit columns) and their `energies`. Materializes the energy
design matrix `X_E[config, salc] = evaluate(salc, config)`.
"""
struct SCEDataset
    basis::SCEBasis
    configs::Vector{Matrix{Float64}}
    X_E::Matrix{Float64}
    y_E::Vector{Float64}
end

function SCEDataset(basis::SCEBasis, configs::AbstractVector, energies::AbstractVector)::SCEDataset
    salcs = basis.salcs.salcs
    cfgs = [Matrix{Float64}(c) for c in configs]
    n = length(cfgs)
    m = length(salcs)
    X = Matrix{Float64}(undef, n, m)
    @inbounds for i = 1:n, j = 1:m
        X[i, j] = evaluate(salcs[j], cfgs[i])
    end
    return SCEDataset(basis, cfgs, X, collect(Float64, energies))
end

"""
    SCEModel

A lightweight predictor: the basis plus the fitted reference energy `j0` and SALC
coefficients `jphi` (addressed by `keys`).
"""
struct SCEModel
    basis::SCEBasis
    j0::Float64
    jphi::Vector{Float64}
    keys::Vector{SALCKey}
end

"""
    SCEFit

A fitted SCE model: the dataset, fitted `j0`/`jphi`, the estimator, and residuals.
"""
struct SCEFit
    dataset::SCEDataset
    j0::Float64
    jphi::Vector{Float64}
    estimator::AbstractEstimator
    residuals::Vector{Float64}
end

"""
    fit(SCEFit, dataset, estimator) -> SCEFit

Fit the SCE coefficients. The energy design matrix is column-centered (so `j0` is
recovered analytically as the mean residual, independent of the estimator), then
the centered problem is handed to `solve_coefficients`.
"""
function fit(::Type{SCEFit}, dataset::SCEDataset, estimator::AbstractEstimator)::SCEFit
    X = dataset.X_E
    y = dataset.y_E
    isempty(y) && throw(ArgumentError("dataset has no observations"))
    xbar = vec(mean(X; dims = 1))
    ybar = mean(y)
    Xc = X .- xbar'
    yc = y .- ybar
    jphi = solve_coefficients(estimator, Xc, yc)
    j0 = ybar - dot(xbar, jphi)
    residuals = y .- (j0 .+ X * jphi)
    return SCEFit(dataset, j0, jphi, estimator, residuals)
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
