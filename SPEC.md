# MagestyRebuild.jl — specification and architecture

Work-in-progress rebuild of Magesty.jl. This file tracks the **realized**
architecture as milestones land; the full design rationale lives in the
implementation plan.

## Pipeline

```
geometry → symmetry → clusters → basis (SALC) → design matrix → fit → predict
```

## Module layout (`src/`)

Single top-level `module MagestyRebuild`, small files included in dependency
order; the two self-contained numeric kernels (`Harmonics`, `AngularMomentum`)
are nested submodules. Heavy / optional dependencies (Spglib, GLMNet, EzXML) live
in `ext/` package extensions; the core loads without them.

## Extension seams (contracts)

- **DFT sources**: `read_configs(src::AbstractDFTSource) -> Vector{<:AbstractTrainingDatum}`.
- **Estimators**: subtype `AbstractEstimator` + `solve_coefficients(est, X, y) -> jphi`
  (the `(X, y)` is already column-centered and row-scaled — add no intercept, do
  not re-weight). Types live in core; solver methods needing GLMNet live in `ext/`.
- **Symmetry**: `analyze_symmetry(backend::AbstractSymmetryBackend, crystal; tol) -> SpaceGroup`.

## Realized so far

### geometry (M1)
- `Lattice(vectors; pbc)` — columns are lattice vectors; `reciprocal = inv`,
  `interplanar_spacing(lat, i) = 1/‖row_i(reciprocal)‖`.
- `Crystal(lattice, frac_positions, species, species_labels)` — fractional coords
  wrapped to `[0,1)` on periodic axes (inner constructor); `num_atoms`,
  `cartesian_positions`.
- `build_neighbor_list(crystal, cutoff) -> NeighborList` — cutoff-driven image
  range `N_d = ceil(cutoff·‖b_d‖)`; `NeighborPair` retains the integer lattice
  `shift` (R) for later reciprocal-space / spin-spiral rows. Validated against an
  independent over-large-shell brute force (cubic multi-shell + sheared triclinic).

### basis — `Harmonics` submodule (M2)
- `Harmonics.Zlm(l, m, u)` — real tesseral harmonic (Drautz convention, per-site
  `(4π)^(−1/2)`); `Harmonics.grad_Zlm(l, m, u)` — tangent-projected Cartesian
  gradient `∂Z − u(u·∂Z)`. `lm_index(l, m) = l²+l+m+1`, `num_lm(lmax)`. Legendre
  primitive from `LegendrePolynomials.dnPl`.
- Validated by: closed-form standard solid harmonics (`l ≤ 2`), gradient tangency
  + on-sphere central difference, and bit-for-bit agreement with Magesty's
  `TesseralHarmonics` (oracle).

### basis — `AngularMomentum` submodule (M3, M4)
- `clebsch_gordan` (Racah formula); `wignerD_real(l, R)` (real Wigner-D, built
  from the package's own `Zₗₘ` by an exact least-squares fit — handles improper
  rotations); `coupling_paths`, `coeff_tensor_complex`, `c2r_matrix`,
  `complex_to_real_tensor`, `build_real_bases`. `CoupledBasis{R}` / `coupled_bases`.
- Validated: CG vs WignerSymbols + orthonormality; Wigner-D functional identity +
  special cases + Magesty Δl; coupled-tensor realness, Frobenius norm² = 2Lf+1,
  rotational equivariance `f(R·e) = Δ^Lf(R)·f(e)`, and Magesty `build_all_real_bases`.

### symmetry (M5)
- `AbstractSymmetryBackend` + `analyze_symmetry(backend, crystal; tol) -> SpaceGroup`;
  in-tree `NoSymmetry` (P1); `SpglibBackend` type in core, method in
  `ext/MagestyRebuildSpglibExt`. `map_sym` derived in-tree (shared by all backends).

### clusters (M6)
- `ClusterMember` (atoms + per-site lattice `shift` R), `ClusterOrbit`, `ClusterSet`;
  `build_clusters` reduces neighbor-pair / atom candidates to symmetry orbits via a
  canonical key built from the site-image map `(b, τ + W·R)`.

### SALC basis (M7)
- `build_salc_basis` projects each orbit's representative tensor onto the trivial
  irrep of its site stabilizer (per-`Lf` projector `ε_g·wignerD_real(Lf, R_g)`),
  deterministic axis-pivoted gauge, transport to members. `SALCKey` (canonical
  column address) + `SALCBasis` (sorted keys + fingerprint). `evaluate(salc, e)`.
- Validated by the ground-truth tests with non-collinear spins, **all `Lf`**:
  space-group invariance `Φ(g·e)=Φ(e)` and time-reversal evenness. The projector
  represents each op on the `Lf` multiplet by its proper part, so symmetry-
  forbidden odd-`Lf` channels are dropped and allowed ones kept.

### fitting + SCE API (M8, M9)
- `Interaction`, `SCEBasis`, `SCEDataset` (energy design matrix `X_E`),
  `SCEModel`/`SCEFit`, `fit(SCEFit, dataset, estimator)`, `predict_energy`,
  `coef`/`intercept`/`nobs`/`r2_energy`/`rmse_energy`. `AbstractEstimator` with
  `OLS`/`Ridge` (analytic `j0`; centered-`X` `solve_coefficients` contract).
- Validated: OLS recovers an in-span target (R²=1); **a Heisenberg chain fit
  recovers `J = 2√3·jphi` to rtol 1e-8** (the v0 done-line, oracle).

## Not yet implemented (v0 follow-ups)
- Torque design matrix `X_T` (energy-only fitting for now).
- 3-body clusters (the `N`-generic machinery is in place; enumeration capped at 2).
- Extensions for GLMNet estimators / VASP I/O / Sunny export; `Tables.jl` results;
  basis persistence.

## Oracle environment (`test/oracle/`)

A separate Julia env (`[sources]`-deving both `MagestyRebuild` and a pinned
`Magesty.jl`) cross-checks convention-fixed kernels and gauge-invariant
aggregates / predictions against Magesty. The core suite never depends on
Magesty. Run: `julia --project=test/oracle test/oracle/runtests.jl`.

(Sections for M3–M10 added as they land.)
