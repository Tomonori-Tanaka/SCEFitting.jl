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
refit
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
dof
r2_energy
rmse_energy
r2_torque
rmse_torque
rss_energy
rss_torque
residuals_energy
residuals_torque
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

## Mean-field spin-configuration sampling

Generate physically representative finite-temperature spin configurations for SCE
training, from the single-site mean-field (MFA) distribution at a controlled reduced
temperature `τ = T/T_MF`. The `MFASampler` has four constructions of increasing fidelity:
the single global isotropic sampler (`MFASampler(reference)`); the multi-sublattice
isotropic and tensorial samplers from an `ExchangeModel` (Heisenberg / DMI / anisotropic
exchange + single-ion); and the full multipole sampler from a fitted `SCEModel` (all
clusters and `l`, higher-order / many-body). `sample` is the one verb.

```@docs
AbstractSampler
MFASampler
ExchangeModel
MultipoleField
sample
MFASample
mfa_sublattice_m
mfa_temperature_scale
thermal_averaged_m
tau_from_magnetization
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
