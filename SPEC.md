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

## Oracle environment (`test/oracle/`)

A separate Julia env (`[sources]`-deving both `MagestyRebuild` and a pinned
`Magesty.jl`) cross-checks convention-fixed kernels and gauge-invariant
aggregates / predictions against Magesty. The core suite never depends on
Magesty. Run: `julia --project=test/oracle test/oracle/runtests.jl`.

(Sections for M3–M10 added as they land.)
