# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/); this package predates a tagged
release, so everything lives under *Unreleased*.

## [Unreleased]

### Added — torque observable (energy + torque co-fit)

- **Gradient kernel** (`basis/salc.jl`): `accumulate_grad!` sums `jϕ·∂Φ/∂e_a` per
  site (product rule over cluster members), sharing the energy kernel's `μ`-mapping
  and `(4π)^(N/2)` scale.
- **Torque design matrix** `X_T` and per-config torque targets via the four-argument
  `SCEDataset(basis, configs, energies, torques)`; entry `(4π)^(N/2)·(e_a × ∂Φ/∂e_a)`,
  rows flattened config/atom/`xyz`.
- **Energy+torque co-fit**: `fit(SCEFit, dataset, est; torque_weight = w)` minimizes
  `(1−w)·MSE_E + w·MSE_T` by per-block whitening; `j0` stays analytic from the
  energy block.
- **Prediction & metrics**: `predict_torque`, `r2_torque`, `rmse_torque`,
  `has_torque`. `predict_torque = e × ∇(predict_energy)` by construction (validated
  by on-sphere finite differences and the Heisenberg closed form).

### Added — v0 vertical slice (energy fitting end-to-end)

- **Geometry**: `Lattice`, `Crystal`, and a generalized cutoff `NeighborList`
  (`build_neighbor_list`) that replaces a fixed 27-cell image grid; `NeighborPair`
  retains the inter-site lattice translation `R`.
- **Harmonics** (`basis/Harmonics.jl`): real tesseral `Zₗₘ` + tangent-projected
  gradient (Drautz convention).
- **Angular momentum** (`basis/AngularMomentum.jl`): Clebsch–Gordan (Racah),
  `wignerD_real` (least-squares from the package's own `Zₗₘ`), complex→real coupled
  tensors; `CoupledBasis{R}`.
- **Symmetry**: `AbstractSymmetryBackend` with in-tree `NoSymmetry` and a
  `SpglibBackend` whose method lives in `ext/MagestyRebuildSpglibExt`.
- **Clusters**: orbit reduction with `R`-carrying members (`build_clusters`).
- **SALC basis**: orbit–stabilizer projector (isotropic and anisotropic channels),
  deterministic gauge, canonical `SALCKey` column addressing (`build_salc_basis`).
- **Fitting / API**: `Interaction`, `SCEBasis`, `SCEDataset` (energy design
  matrix), `SCEModel`/`SCEFit`, `fit`, `predict_energy`, `OLS`/`Ridge`,
  `coef`/`intercept`/`nobs`/`r2_energy`/`rmse_energy`.
- A runnable Heisenberg-chain example (`examples/heisenberg_chain.jl`) that
  recovers `J`; an oracle test environment (`test/oracle/`) validating against a
  pinned Magesty.jl.

### Notes

- Cluster enumeration is capped at 2-body (the angular-momentum machinery is
  `N`-generic).
- Bit-for-bit agreement with Magesty.jl is explicitly not a goal — see
  `docs/design-notes.md`.
