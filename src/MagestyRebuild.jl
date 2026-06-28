"""
    MagestyRebuild

Clean, extensible, Julia-native rebuild of `Magesty.jl`: fit spin-cluster-expansion
(SCE) models to noncollinear DFT data. The numerical core is reimplemented from
scratch; `Magesty.jl` serves only as a pinned numerical oracle in `test/oracle/`.

Work in progress (v0 vertical slice). See `SPEC.md` for the realized architecture.
"""
module MagestyRebuild

using LinearAlgebra: norm, det, I, eigen, Symmetric, Diagonal, dot, cross, mul!
using StaticArrays
using Statistics: mean
using Random: AbstractRNG, default_rng
import TOML
import Tables

# --- geometry ---
include("geometry/lattice.jl")
include("geometry/crystal.jl")
include("geometry/neighborlist.jl")

# --- symmetry (backend-pluggable; Spglib lives in an extension) ---
include("symmetry/types.jl")
include("symmetry/backend.jl")

# --- basis: numeric kernels (self-contained submodules) ---
include("basis/Harmonics.jl")
include("basis/AngularMomentum.jl")
include("basis/coupledbasis.jl")

# --- clusters: enumeration + symmetry orbit reduction ---
include("clusters/enumerate.jl")
include("clusters/orbits.jl")

# --- SALC basis: symmetry-adapted, time-reversal-even invariants ---
include("basis/salc.jl")
include("basis/salcbasis.jl")

# --- mean-field spin-configuration sampling (docs/specs/mfa-sampling.md) ---
# P0: the single-site engine (potential, vMF / Metropolis draws, sphere quadrature).
# P2/P3: the ExchangeModel carrier + coupled self-consistency (before MFASampler, which
# references the ExchangeModel type). P1/P2/P3: the `MFASampler` + the `sample` verb.
include("sampling/site_engine.jl")
include("sampling/exchange.jl")
include("sampling/mfa_sampler.jl")

# --- fitting + high-level SCE API ---
include("fitting/estimators.jl")
include("sce/model.jl")

# --- tabular results (Tables.jl source) ---
include("sce/coeftable.jl")

# --- Sunny export: conversion math in core, Sunny.System assembly in the extension ---
include("sce/sunny.jl")

# --- mean-field sampler: extract an ExchangeModel from a fitted SCE (depends on the Sunny
# bilinear / single-ion extraction `_sunny_supercell_terms` above).
include("sampling/exchange_from_sce.jl")

# --- persistence (format-agnostic schema, serialized as TOML) + TOML input files ---
include("sce/persist.jl")
include("sce/input.jl")

# --- DFT data sources + VASP I/O ---
include("io/dftsource.jl")
include("io/vasp.jl")

export Lattice, Crystal, num_atoms, cartesian_positions, interplanar_spacing
export NeighborPair, NeighborList, build_neighbor_list
export AbstractImageSelection, MinimumImage, AllImages
export SymOp, SpaceGroup, AbstractSymmetryBackend, NoSymmetry, SpglibBackend,
    analyze_symmetry, n_ops
export ClusterMember, ClusterOrbit, ClusterSet, build_clusters
export SALCKey, SALC, SALCBasis, build_salc_basis, evaluate
export Interaction, SCEBasis, SCEDataset, SCEModel, SCEFit, fit, refit, nsalc, read_input
export predict_energy, predict_torque, has_torque
export AbstractEstimator, OLS, Ridge, ElasticNet, Lasso, AdaptiveLasso, AdaptiveRidge,
    PrecomputedPilot, solve_coefficients
export coef, intercept, nobs, dof, r2_energy, rmse_energy, r2_torque, rmse_torque,
    rss_energy, rss_torque, residuals_energy, residuals_torque
export coeftable, SCECoefficients
export to_sunny
# Mean-field spin-configuration sampling (docs/specs/mfa-sampling.md).
export AbstractSampler, MFASampler, MFASample, ExchangeModel, MultipoleField, sample,
    mfa_temperature_scale, mfa_sublattice_m, thermal_averaged_m, tau_from_magnetization
# DFT data I/O: only the code-agnostic boundary is exported; per-code adapters are
# namespaced submodules (e.g. `MagestyRebuild.VASP.read_poscar`), so adding a code
# touches neither the core nor this export list.
export AbstractDFTSource, AbstractTrainingDatum, SpinDatum, read_configs

end # module MagestyRebuild
