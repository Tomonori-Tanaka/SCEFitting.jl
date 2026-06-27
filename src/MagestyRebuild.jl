"""
    MagestyRebuild

Clean, extensible, Julia-native rebuild of `Magesty.jl`: fit spin-cluster-expansion
(SCE) models to noncollinear DFT data. The numerical core is reimplemented from
scratch; `Magesty.jl` serves only as a pinned numerical oracle in `test/oracle/`.

Work in progress (v0 vertical slice). See `SPEC.md` for the realized architecture.
"""
module MagestyRebuild

using LinearAlgebra: norm, det, I
using StaticArrays

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

export Lattice, Crystal, num_atoms
export NeighborPair, NeighborList, build_neighbor_list
export SymOp, SpaceGroup, AbstractSymmetryBackend, NoSymmetry, SpglibBackend,
    analyze_symmetry
export ClusterMember, ClusterOrbit, ClusterSet, build_clusters

end # module MagestyRebuild
