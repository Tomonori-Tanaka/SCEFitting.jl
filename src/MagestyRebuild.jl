"""
    MagestyRebuild

Clean, extensible, Julia-native rebuild of `Magesty.jl`: fit spin-cluster-expansion
(SCE) models to noncollinear DFT data. The numerical core is reimplemented from
scratch; `Magesty.jl` serves only as a pinned numerical oracle in `test/oracle/`.

Work in progress (v0 vertical slice). See `SPEC.md` for the realized architecture.
"""
module MagestyRebuild

using LinearAlgebra: norm, det
using StaticArrays

# --- geometry ---
include("geometry/lattice.jl")
include("geometry/crystal.jl")
include("geometry/neighborlist.jl")

# --- basis: numeric kernels (self-contained submodules) ---
include("basis/Harmonics.jl")
include("basis/AngularMomentum.jl")

export Lattice, Crystal, num_atoms
export NeighborPair, NeighborList, build_neighbor_list

end # module MagestyRebuild
