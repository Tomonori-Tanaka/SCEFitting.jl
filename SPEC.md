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
- **Estimators**: subtype `AbstractEstimator` + `solve_coefficients(est, X, y; groups)
  -> jphi` (the `(X, y)` is already column-centered and row-scaled — add no intercept,
  do not re-weight; `groups` labels rows from the same sample for grouped resampling).
  Types live in core; solver methods needing GLMNet live in `ext/`.
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
- `build_neighbor_list(crystal, cutoff, selection)` — periodic-image selection
  (`AbstractImageSelection`): `MinimumImage()` (the SCE-fitting default) keeps only the
  minimum-image, plain-PBC-resolvable pairs of the Wigner–Seitz cell — with boundary
  ties (`L/2` faces / edges / `(L/2,L/2,L/2)` corners) and **no `i==j` self-pairs**
  (same spin ⇒ not an independent pair) — over an adaptive, skew-safe image box;
  `cutoff` may be `Inf` (whole WS cell). `AllImages()` is the every-image enumeration
  above (finite cutoff only), the generalized-Bloch / spin-spiral seam.

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
  `candidate_clusters` enumerates `N`-body **edge-admissible cliques** (`N = 2` is the
  directed neighbor pairs), `build_clusters` reduces them to symmetry orbits via a
  canonical key built from the site-image map `(b, τ + W·R)`. An edge is admissible
  under `MinimumImage` iff it sits at its atom-pair minimum-image distance (and the
  clique's atoms are distinct, so a cluster never reuses an atom's image — that would
  alias a lower-body term); under `AllImages` iff within the radial cutoff. The
  `selection` is threaded from `SCEBasis` so the neighbor list and clusters agree.

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
- **TOML input** (`sce/input.jl`): `read_input(path) -> (; crystal, interaction, backend, tol, images)`
  and `SCEBasis(path::AbstractString; backend, tol, images)` build a basis from a human-authored
  `input.toml` (`[structure]` inline crystal, `[interaction]` with optional `images`
  (`"minimum_image"` default / `"all_images"`) and `pair_cutoff = inf` for the full WS
  cell, optional `[symmetry]`);
  keyword arguments override the file's backend/tol. Training data and the estimator
  stay in Julia (mirrors the basis/data separation).
- **Tabular results** (`sce/coeftable.jl`): `coeftable(fit | model) -> SCECoefficients`
  is a **Tables.jl** source — one row per SALC (`body`, `orbit_id`, `ls` as a comma
  string, `Lf`, `block`, `J`) — so it drops into `DataFrame` / `CSV.write` /
  `Arrow.write`. The library owns the internal-storage → labeled-row mapping; the caller
  brings the table/IO package. `j0` is the intercept (`intercept(c)`), not a row.
  Tables.jl is a lightweight core dep.
- **DFT data sources** (`io/dftsource.jl`, `io/vasp.jl`): a code-agnostic boundary —
  `SpinDatum` (energy + spin directions + magmoms + constraining field + the derived
  torque target `τ_a = m_a × B_a`, the physical / Landau–Lifshitz torque) and
  `read_configs(src::AbstractDFTSource) -> Vector{SpinDatum}`, with `SCEDataset(basis, src)`
  going source → dataset. Per-code adapters are **namespaced submodules** kept out of the
  core: `MagestyRebuild.VASP` reads POSCAR/CONTCAR (`read_poscar`/`write_poscar`) and
  constrained-noncollinear OSZICARs (`Oszicar`, with SAXIS rotation). Adding a DFT code is
  one sibling submodule — neither the core nor its export list changes. Cross-checked
  bit-for-bit against Magesty's parsers (oracle).
- Validated: basis / model / fit round-trips (predictions bit-identical, coefficients
  re-paired by key under scrambled order, multi-op space-group ops, empty basis), input
  parsing + defaults + keyword overrides + error paths.

### Sunny.jl export (M11)
- **Conversion core** (`sce/sunny.jl`, dependency-free): `_l1_pair_matrix` / `_l2_onsite_matrix`
  turn a folded tesseral tensor into a Cartesian exchange / single-ion matrix
  (`eₐ'·M·e_b = Σ folded·Z·Z`); `_classify_salc` keeps only `ls=[1,1]` pairs and `ls=[2]`
  single-ion (`ls=[0…]` → `j0`, the rest skipped + reported). `_sunny_supercell_terms`
  folds the directed members into one matrix per undirected supercell bond, and
  `_sunny_primitive` unfolds them onto the chemical primitive cell recovered from the
  pure translations (one Sunny bond per primitive bond, a `clean` flag for the fallback).
  `_reconstruct_energy ≈ predict_energy − j0` gates the whole chain **without Sunny**.
- **Extension** (`ext/MagestyRebuildSunnyExt`, loaded by `using Sunny`):
  `to_sunny(model; spins, g, mode, placement) -> Sunny.System` builds a real `System`
  (`set_exchange_at!`/`set_onsite_coupling_at!` on the supercell, or
  `set_exchange!`/`Bond` on the primitive cell), rescaling `J = M/(SₐS_b)`; its classical
  energy reproduces `predict_energy − j0`. `placement = :explicit` (exact, folded
  dispersion) / `:primitive` (unfolded) / `:auto`.
- Validated: Sunny-free conversion + primitive-fold tests (`test/unit/test_sunny.jl`,
  main suite) and a separate `test/sunny/` environment (real `Sunny.System` energy vs SCE
  for both routes, the primitive system reshaped back to the supercell, mode/spin handling,
  the skip warning).

### GLMNet estimators (M12)
- **Types in core** (`fitting/estimators.jl`): `ElasticNet(; alpha, lambda, standardize,
  nfolds, select, seed, nlambda)` and `Lasso(; …)` (= `alpha = 1`). Named and dispatched
  on without the dependency; argument validation lives here.
- **Extension** (`ext/MagestyRebuildGLMNetExt`, loaded by `using GLMNet`):
  `solve_coefficients(::ElasticNet, X, y; groups)`. GLMNet minimizes
  `(1/2n)‖y−Xβ‖² + λ[(1−α)/2‖β‖₂² + α‖β‖₁]` on the column-centered `(X, y)` with
  `intercept = false` (so `j0` stays analytic) and `standardize` (penalty acts per-column
  at `λ·std`). `lambda = nothing` ⇒ K-fold CV over GLMNet's λ path, `select` =
  `:lambda_min`/`:lambda_1se`; folds are **grouped by configuration** (via `groups`) and
  seeded deterministically (`hash`-ranked, no RNG dependency). A numeric `lambda` skips CV.
- Validated in a separate `test/glmnet/` environment: tiny-λ ≈ OLS, the analytic-`j0`
  centering invariant, CV support recovery + sparsity, `:lambda_1se` shrinkage,
  reproducibility, ElasticNet ≠ Lasso, and an energy+torque grouped-CV co-fit. Core-only
  construction / validation / deferred-backend error in `test/unit/test_fit.jl`.

## Not yet implemented (v0 follow-ups)
- The v0 slice is feature-complete; no estimator/observable/IO follow-ups outstanding.

## Oracle environment (`test/oracle/`)

A separate Julia env (`[sources]`-deving both `MagestyRebuild` and a pinned
`Magesty.jl`) cross-checks convention-fixed kernels and gauge-invariant
aggregates / predictions against Magesty. The core suite never depends on
Magesty. Run: `julia --project=test/oracle test/oracle/runtests.jl`.

(Sections for M3–M10 added as they land.)
