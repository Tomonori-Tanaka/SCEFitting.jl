"""
    SCEFitting

Clean, extensible, Julia-native rebuild of `Magesty.jl`: fit spin-cluster-expansion
(SCE) models to noncollinear DFT data. The numerical core is reimplemented from
scratch; `Magesty.jl` serves only as a pinned numerical oracle in `test/oracle/`.

The v0 feature set (geometry, symmetry, cluster/SALC basis, fitting, prediction,
diagnostics, persistence, Sunny export, introspection) is realized; see `SPEC.md`.
"""
module SCEFitting

using LinearAlgebra: norm, det, I, eigen, Symmetric, Diagonal, dot, cross
using StaticArrays
using Statistics: mean
using Random: AbstractRNG, default_rng
import TOML
import Tables
# Extend the StatsAPI generics rather than shadow them, so `coef` / `fit` / `nobs` /
# `dof` compose with the StatsBase / GLM ecosystem instead of clashing on `using`.
import StatsAPI: coef, fit, nobs, dof

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

# --- fitting + high-level SCE API ---
include("fitting/estimators.jl")
include("sce/model.jl")          # pipeline types + constructors + config validation
include("fitting/design.jl")     # design-matrix assembly (X_E / X_T)
include("fitting/fit.jl")        # fit / refit / predict
include("fitting/diagnostics.jl")  # coef / intercept / residuals / R² / RMSE

# --- tabular results (Tables.jl source) ---
include("sce/coeftable.jl")

# --- external-engine interop: conversion math in core, engine assembly in extensions ---
# Sunny export (Sunny.System assembled in SCEFittingSunnyExt). The general Cartesian
# bilinear / single-ion extraction `_bilinear_terms` lives here too.
include("interop/sunny.jl")

# --- fitted-model introspection: a code-neutral view of the multipole / bilinear terms
# (consumed by downstream packages such as the SCETools.jl samplers); depends on the
# bilinear / single-ion extraction `_bilinear_terms` above.
include("sce/introspect.jl")

# --- I/O: persistence (TOML model schema), TOML input files, and the code-agnostic DFT
# data boundary. Concrete DFT-code adapters (e.g. the VASP reader/writer) live in
# SCETools.jl; the SCE pipeline only ever sees `SpinDatum` / `SCEDataset`.
include("io/persist.jl")
include("io/input.jl")
include("io/dftsource.jl")

# --- Public API (exported) --------------------------------------------------------
# The fitting workflow a user reaches for. Construction internals (cluster / neighbor /
# SALC builders, symmetry analysis) are *public but unexported* — see the block below.

# geometry the user builds
export Lattice, Crystal, n_atoms, cartesian_positions
# how periodic images / symmetry are chosen (passed into `SCEBasis`)
export AbstractImageSelection, MinimumImage, AllImages
export AbstractSymmetryBackend, NoSymmetry, SpglibBackend
# the SCE pipeline
export BasisSpec, SCEBasis, SCEDataset, SCEPredictor, SCEFit, fit, refit, n_salcs, read_setup
export predict_energy, predict_torque, has_torque
# estimators
export AbstractEstimator, OLS, Ridge, ElasticNet, Lasso, AdaptiveLasso, AdaptiveRidge,
    PrecomputedPilot
# fit diagnostics
export coef, intercept, nobs, dof, r2_energy, rmse_energy, r2_torque, rmse_torque,
    rss_energy, rss_torque, residuals_energy, residuals_torque
export coeftable, SCECoefficients
export to_sunny
# Fitted-model introspection: a code-neutral view of the multipole / bilinear terms of a
# fitted SCE, the stable contract downstream packages (e.g. the SCETools.jl mean-field
# samplers) read instead of the SALC-basis internals.
export MultipoleTerm, multipole_terms, bilinear_terms
# DFT data I/O: only the code-agnostic boundary is exported; per-code adapters are
# namespaced submodules (e.g. `SCEFitting.VASP.read_poscar`), so adding a code
# touches neither the core nor this export list.
export AbstractDFTSource, SpinDatum, read_configs

# Deprecated alias: `BasisSpec` was named `Interaction` before v0.2; kept (and exported)
# for one minor version so existing scripts keep working.
const Interaction = BasisSpec
export Interaction

# --- Public, unexported -----------------------------------------------------------
# Reachable as `SCEFitting.<name>` (and documented), but kept out of the flat `using`
# namespace: the `SCEBasis` constructor already drives them for you. Power users and the
# test suite reach them by qualification.
#
#   geometry/neighbors : build_neighbor_list, NeighborPair, NeighborList,
#                        interplanar_spacing
#   symmetry           : analyze_symmetry, n_ops, SymOp, SpaceGroup, AbstractTrainingDatum
#   clusters           : build_clusters, ClusterMember, ClusterOrbit, ClusterSet
#   SALC basis         : build_salc_basis, evaluate_salc, salcs, SALC, SALCKey, SALCBasis
#   estimators         : solve_coefficients

end # module SCEFitting
