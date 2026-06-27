"""
    MagestyRebuild

Clean, extensible, Julia-native rebuild of `Magesty.jl`: fit spin-cluster-expansion
(SCE) models to noncollinear DFT data. The numerical core is reimplemented from
scratch; `Magesty.jl` serves only as a pinned numerical oracle in `test/oracle/`.

Work in progress (v0 vertical slice). See `SPEC.md` for the realized architecture.
"""
module MagestyRebuild

using LinearAlgebra: norm, det, I, eigen, Symmetric, dot
using StaticArrays
using Statistics: mean

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
include("sce/model.jl")

export Lattice, Crystal, num_atoms
export NeighborPair, NeighborList, build_neighbor_list
export SymOp, SpaceGroup, AbstractSymmetryBackend, NoSymmetry, SpglibBackend,
    analyze_symmetry
export ClusterMember, ClusterOrbit, ClusterSet, build_clusters
export SALCKey, SALC, SALCBasis, build_salc_basis
export Interaction, SCEBasis, SCEDataset, SCEModel, SCEFit, fit, predict_energy
export AbstractEstimator, OLS, Ridge, solve_coefficients
export coef, intercept, nobs, r2_energy, rmse_energy

end # module MagestyRebuild
