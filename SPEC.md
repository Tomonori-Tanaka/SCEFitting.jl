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
are nested submodules. Heavy / optional dependencies (Spglib, GLMNet) live in
`ext/` package extensions; the core loads without them. Persistence and input
files use the stdlib `TOML` (no external dependency).

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

### clusters (M6) — arbitrary body order
- `ClusterMember` (atoms + per-site lattice `shift` R), `ClusterOrbit`, `ClusterSet`;
  `candidate_clusters` enumerates `N`-body **pairwise-within-cutoff cliques** (`N = 2`
  is the directed neighbor pairs), `build_clusters` reduces them to symmetry orbits
  via a canonical key built from the site-image map `(b, τ + W·R)`.

### SALC basis (M7) — arbitrary body order
- `build_salc_basis` projects each orbit's representative onto the trivial irrep of
  its site stabilizer. At `N ≥ 3` a stabilizer op can permute equivalent sites, mixing
  coupling paths and (for unequal `l`) `l`-orderings, so the projection runs over the
  **combined (ordering × path × `Mf`) space**: each op acts by rotating every site axis
  (`wignerD_real(l_i, R)`) and relabeling axes by `invperm(perm)`, with the action
  matrix read off by contraction against the orthonormal coupled tensors. Deterministic
  axis-pivoted gauge; per-ordering fold into a multi-term SALC (one `SALCTerm` per
  ordering). `SALCKey` (canonical, injective column address — `block` runs across split
  ordering orbits) + `SALCBasis` (sorted keys + fingerprint). `evaluate(salc, e)`.
- Validated by the ground-truth tests with non-collinear spins, **all `Lf`, all body
  orders**: space-group invariance `Φ(g·e)=Φ(e)`, time-reversal evenness, linear
  independence; projector eigenvalues exactly 0/1. Improper-op parity is handled
  automatically (site-axis rotation by full `R` + even `Σl`). Cross-validated against
  Magesty: per-`(body, ls, Lf)` invariant-subspace dimensions agree through 3-body.

### fitting + SCE API (M8, M9)
- `Interaction`, `SCEBasis`, `SCEDataset` (energy design matrix `X_E`, and the
  torque design matrix `X_T` via the four-argument form), `SCEModel`/`SCEFit`,
  `fit(SCEFit, dataset, estimator; torque_weight)`, `predict_energy`/`predict_torque`,
  `coef`/`intercept`/`nobs`/`r2_energy`/`rmse_energy`/`r2_torque`/`rmse_torque`/
  `has_torque`. `AbstractEstimator` with `OLS`/`Ridge` (analytic `j0`; centered-`X`
  `solve_coefficients` contract).
- Validated: OLS recovers an in-span target (R²=1); **a Heisenberg chain fit
  recovers `J = 2√3·jphi` to rtol 1e-8** (the v0 done-line, oracle); torque is the
  exact derivative of the energy surface (on-sphere finite differences, equivariance,
  Heisenberg closed form), energy+torque co-fit recovers an in-span model.

### persistence + TOML input (M10)
- **Persistence** (`sce/persist.jl`): `MagestyRebuild.save(path, x)` and
  `MagestyRebuild.load(SCEBasis | SCEModel, path)` serialize a self-contained,
  human-readable **TOML** document — the crystal, the space-group ops, the interaction,
  and the *full* SALC basis (every member / term / folded tensor); a model adds `j0`
  and per-`SALCKey` coefficients. Reload reconstructs the basis verbatim (no
  re-projection) and re-pairs coefficients to the basis **by key**, not by position.
  The `struct ⇄ Dict` schema (`_to_doc` / `*_from_doc`) is format-agnostic and tested
  without any serializer; TOML is stdlib (no dep) and round-trips `Float64` exactly.
  `save` / `load` are **unexported** (call qualified) to avoid clashing with
  `FileIO`/`JLD2`.
- **TOML input** (`sce/input.jl`): `read_input(path) -> (; crystal, interaction, backend, tol)`
  and `SCEBasis(path::AbstractString; backend, tol)` build a basis from a human-authored
  `input.toml` (`[structure]` inline crystal, `[interaction]`, optional `[symmetry]`);
  keyword arguments override the file's backend/tol. Training data and the estimator
  stay in Julia (mirrors the basis/data separation).
- **Tabular results** (`sce/coeftable.jl`): `coeftable(fit | model) -> SCECoefficients`
  is a **Tables.jl** source — one row per SALC (`body`, `orbit_id`, `ls` as a comma
  string, `Lf`, `block`, `J`) — so it drops into `DataFrame` / `CSV.write` /
  `Arrow.write`. The library owns the internal-storage → labeled-row mapping; the caller
  brings the table/IO package. `j0` is the intercept (`intercept(c)`), not a row.
  Tables.jl is a lightweight core dep.
- Validated: basis / model / fit round-trips (predictions bit-identical, coefficients
  re-paired by key under scrambled order, multi-op space-group ops, empty basis), input
  parsing + defaults + keyword overrides + error paths.

## Not yet implemented (v0 follow-ups)
- Extensions for GLMNet estimators / VASP I/O / Sunny export.

## Oracle environment (`test/oracle/`)

A separate Julia env (`[sources]`-deving both `MagestyRebuild` and a pinned
`Magesty.jl`) cross-checks convention-fixed kernels and gauge-invariant
aggregates / predictions against Magesty. The core suite never depends on
Magesty. Run: `julia --project=test/oracle test/oracle/runtests.jl`.

(Sections for M3–M10 added as they land.)
