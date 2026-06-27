# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/); this package predates a tagged
release, so everything lives under *Unreleased*.

## [Unreleased]

### Added — arbitrary body order (N-body clusters)

- **Cluster enumeration** (`clusters/enumerate.jl`): `candidate_clusters` generalized
  to `N`-body pairwise-within-cutoff cliques (`N = 2` stays the directed neighbor pairs).
- **SALC projection** (`basis/salcbasis.jl`): general site stabilizer with induced
  permutations; projection over the combined (ordering × coupling-path × `Mf`) space
  with the action matrix built by contracting against the package's own orthonormal
  coupled tensors (no 6j/9j). Handles, at `N ≥ 3`, coupling-path mixing and — for
  unequal `l` on symmetry-equivalent sites (e.g. `l=(1,1,2)` on a triangle) —
  `l`-ordering mixing. Improper-op parity is automatic (no proper-part special case).
- **Multi-term SALCs** (`basis/salc.jl`): a SALC carries one `SALCTerm` (own `ls` +
  `folded`) per `l`-ordering; `evaluate` and `accumulate_grad!` loop over terms.
- **`SALCKey`** stays injective when a proper-subgroup stabilizer splits a degenerate
  multiset into several ordering orbits (`block` runs across them).
- Validated by ground-truth invariance / time-reversal / linear-independence tests at
  `N = 3` (incl. the multi-term channel) and an energy+torque 3-body recovery;
  cross-checked against Magesty.jl — per-`(body, ls, Lf)` invariant-subspace dimensions
  agree exactly through 3-body, and Magesty's own SALCs independently pass invariance.

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

- Bit-for-bit agreement with Magesty.jl is explicitly not a goal — see
  `docs/design-notes.md`.
