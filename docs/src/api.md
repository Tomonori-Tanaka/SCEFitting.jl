# API reference

```@meta
CurrentModule = SCEFitting
```

Every exported type and function, grouped by pipeline stage. The headline workflow is
built from `Crystal` + `Interaction` → [`SCEBasis`](@ref) → [`SCEDataset`](@ref) →
[`fit`](@ref) → [`SCEPredictor`](@ref).

```@index
```

## Geometry

```@docs
Lattice
Crystal
n_atoms
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
evaluate_salc
```

## Interaction, basis, dataset, model

```@docs
Interaction
SCEBasis
SCEDataset
SCEPredictor
SCEFit
n_salcs
salcs
read_setup
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

## Fitted-model introspection

A code-neutral view of a fitted [`SCEPredictor`](@ref)'s multipole terms, the stable public
contract downstream packages (e.g. the mean-field samplers in
[`SCETools.jl`](https://github.com/Tomonori-Tanaka/SCETools.jl)) read instead of the
SALC-basis internals.
`multipole_terms` is the general per-term dump; `bilinear_terms` is the bilinear (`ls=[1,1]`)
and single-ion (`ls=[2]`) extraction as Cartesian `3×3` matrices (reusing the validated
Sunny conversion). The tesseral spherical-harmonic kernel `SCEFitting.Harmonics`
(`Zlm`, `lm_index`) is a stable submodule those consumers may also use.

```@docs
MultipoleTerm
multipole_terms
bilinear_terms
```

## DFT data sources

The **code-agnostic boundary** of the training-data input: the SCE pipeline only ever sees
`SpinDatum` / [`SCEDataset`](@ref). Concrete DFT-code adapters (the VASP POSCAR / OSZICAR
reader and INCAR writer) live in the companion
[SCETools.jl](https://github.com/Tomonori-Tanaka/SCETools.jl) package
(`SCETools.VASP.read_poscar` / `Oszicar` / `write_inputs`), so adding a code touches neither
the core nor its exports.

```@docs
AbstractDFTSource
AbstractTrainingDatum
SpinDatum
read_configs
```

## Persistence

`save` / `load` are intentionally **not exported** (the names clash with FileIO / JLD2 /
CSV); call them as `SCEFitting.save` / `SCEFitting.load`.

```@docs
SCEFitting.save
SCEFitting.load
```
