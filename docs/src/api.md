# API reference

```@meta
CurrentModule = MagestyRebuild
```

Every exported type and function, grouped by pipeline stage. The headline workflow is
built from `Crystal` + `Interaction` → [`SCEBasis`](@ref) → [`SCEDataset`](@ref) →
[`fit`](@ref) → [`SCEModel`](@ref).

```@index
```

## Geometry

```@docs
Lattice
Crystal
num_atoms
cartesian_positions
interplanar_spacing
```

## Neighbor list and periodic image selection

```@docs
NeighborPair
NeighborList
build_neighbor_list
AbstractImageSelection
MinimumImage
AllImages
```

## Symmetry

```@docs
SymOp
SpaceGroup
AbstractSymmetryBackend
NoSymmetry
SpglibBackend
analyze_symmetry
n_ops
```

## Cluster orbits

```@docs
ClusterMember
ClusterOrbit
ClusterSet
build_clusters
```

## SALC basis

```@docs
SALCKey
SALC
SALCBasis
build_salc_basis
evaluate
```

## Interaction, basis, dataset, model

```@docs
Interaction
SCEBasis
SCEDataset
SCEModel
SCEFit
nsalc
read_input
```

## Fitting

```@docs
fit
AbstractEstimator
OLS
Ridge
ElasticNet
Lasso
AdaptiveLasso
AdaptiveRidge
PrecomputedPilot
solve_coefficients
```

## Prediction

```@docs
predict_energy
predict_torque
has_torque
```

## Diagnostics

```@docs
coef
intercept
nobs
r2_energy
rmse_energy
r2_torque
rmse_torque
```

## Tabular coefficients

```@docs
coeftable
SCECoefficients
```

## Sunny export

```@docs
to_sunny
```

## DFT data sources

```@docs
AbstractDFTSource
AbstractTrainingDatum
SpinDatum
read_configs
```

### VASP

DFT-code readers are namespaced submodules so that the core and its exports do not grow
as codes are added; the SCE pipeline only ever sees the code-agnostic `SpinDatum` /
`SCEDataset` above. The VASP submodule (`MagestyRebuild.VASP`):

```@docs
MagestyRebuild.VASP.read_poscar
MagestyRebuild.VASP.write_poscar
MagestyRebuild.VASP.Oszicar
```

## Persistence

`save` / `load` are intentionally **not exported** (the names clash with FileIO / JLD2 /
CSV); call them as `MagestyRebuild.save` / `MagestyRebuild.load`.

```@docs
MagestyRebuild.save
MagestyRebuild.load
```
