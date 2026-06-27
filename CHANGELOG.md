# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/); this package predates a tagged
release, so everything lives under *Unreleased*.

## [Unreleased]

### Changed — minimum-image periodic resolvability (Wigner–Seitz cell)

- **Periodic-image selection** (`geometry/neighborlist.jl`): `AbstractImageSelection`
  with `MinimumImage` (new default) and `AllImages`. Previously the neighbor list kept
  every periodic image within a *spherical* cutoff, which over-counts beyond `L/2`: a
  farther image of an atom carries the **same spin** as its minimum image, so the two
  interactions are not independently resolvable from a finite supercell (their design
  columns are collinear). `MinimumImage` keeps only the minimum-image, Wigner–Seitz-cell
  pairs — the physically resolvable set, which reaches the body-diagonal corner
  `(L/2,L/2,L/2)` at `√3·L/2`, **not** a sphere of radius `L/2` — with WS-boundary ties
  (faces 2-fold, edges 4-fold, corners 8-fold) kept as distinct members and `i==j`
  self-pairs dropped (same spin ⇒ a constant or a 1-body alias, never an independent
  pair). The image-box search is adaptive (provably sufficient on skewed / non-reduced
  cells). `AllImages` retains the old every-image behavior as the **generalized-Bloch /
  spin-spiral seam**, where `e^{iq·R}` resolves what one supercell cannot.
- **`pair_cutoff = Inf`** (`Interaction`) now means "every resolvable pair" — the whole
  WS cell — under `MinimumImage` (cf. Magesty's `-1` sentinel); `≤ 0` / `NaN` are
  rejected. `SCEBasis(...; images = MinimumImage())` and the `input.toml`
  `[interaction].images` / `pair_cutoff = inf` keys thread the choice through.
- **N-body min-image consistency** (`clusters/enumerate.jl`): a `MinimumImage` clique
  requires every edge at its atom-pair minimum-image distance and all atoms distinct, so
  a cluster is never built from an aliased bond or a reused atom image. The `AllImages`
  edge cutoff now uses the same relative tolerance as the tie band (was an absolute
  `+1e-9`).
- For a cutoff below half the smallest perpendicular cell width, `MinimumImage` and
  `AllImages` coincide (each in-cutoff image is already the minimum), so the existing
  physics is unchanged: the Heisenberg `J = 2√3·jϕ` recovery, the kagome 3-body co-fit,
  and the oracle all use cutoffs in this regime. Bit-for-bit Magesty agreement is not a
  goal; a single-atom cell's self-pair "bond" (1 orbit under Magesty/`AllImages`) is now
  0 orbits under the default `MinimumImage` — a deliberate refinement.
- Validated (`test/unit/test_imageselection.jl`): cutoff `< L/2` equivalence, the cubic
  8-fold corner and `L/2` 2-fold face ties, alias rejection above `L/2`, `AllImages`
  rejects `Inf`, hexagonal / skewed-cell search-box sufficiency, no-self-pair and
  distinct-atom clusters, a full-WS design matrix with full column rank under symmetry,
  `input.toml` `inf` + `images`, and a `pair_cutoff = Inf` persistence round-trip.

### Added — VASP I/O and a code-agnostic DFT-source seam

- **DFT-source boundary** (`io/dftsource.jl`): `AbstractDFTSource` +
  `read_configs(src) -> Vector{SpinDatum}`, the `SpinDatum` training datum (energy,
  `3×n` unit spin directions, moment magnitudes, constraining field, and the derived
  torque target `τ_a = −m_a × B_a`, eV), and `SCEDataset(basis, src | data; use_torque)`.
  This is the *only* thing the SCE pipeline consumes — the originating DFT code is
  irrelevant once you hold a `SpinDatum`/`SCEDataset`.
- **VASP adapter** (`io/vasp.jl`, `module MagestyRebuild.VASP`): `read_poscar` /
  `write_poscar` (POSCAR/CONTCAR ↔ `Crystal`; scaling incl. negative-volume,
  Direct/Cartesian, Selective dynamics, VASP4/5) and `Oszicar`, an `AbstractDFTSource`
  over constrained-noncollinear OSZICARs (energy `F=`/`E0`, `MW_int`/`M_int` moments,
  `lambda*MW_perp` field, SAXIS `Rz(α)·Ry(β)` rotation). Code-specific I/O is a
  **namespaced submodule** kept out of the core: the core and its export list do not grow
  as DFT codes are added; only the code-agnostic boundary is exported.
- A torque-carrying `SCEDataset` rejects all-zero torque targets (no constraining field
  found) so unconstrained data cannot be silently fit as "torque = 0".
- Validated (`test/unit/test_vaspio.jl`): POSCAR (Direct/Cartesian/VASP4/negative-volume/
  selective-dynamics/round-trip), OSZICAR (energy kinds, `mint`, SAXIS, multi-step,
  missing field, errors), the torque formula, and source → `SCEDataset`. The from-scratch
  POSCAR/OSZICAR readers are cross-checked **bit-for-bit against Magesty.jl**'s parsers in
  the oracle (synthetic files, since real VASP outputs are not vendored).

### Added — tabular coefficient output (Tables.jl)

- **`coeftable(f)` / `coeftable(model)` → `SCECoefficients`** (`sce/coeftable.jl`): a
  Tables.jl source — one row per SALC (`body`, `orbit_id`, `ls` as a comma string,
  `Lf`, `block`, `J`) — so the fitted coefficients drop into `DataFrame` /
  `CSV.write` / `Arrow.write`. The library owns the internal-storage → labeled-row
  mapping (the `J` column pairs with the basis keys positionally); the caller brings
  the table/IO package. `j0` is the intercept (`intercept(c)`), not a row. Tables.jl
  is a lightweight core dependency, the same seam that would later open tabular
  training-data ingestion.

### Added — persistence + TOML input files

- **Persistence** (`sce/persist.jl`): `MagestyRebuild.save(path, x)` and
  `MagestyRebuild.load(SCEBasis | SCEModel, path)` serialize a self-contained,
  human-readable **TOML** document — the crystal, the space-group ops, the
  interaction, and the *full* SALC basis (every member / term / folded tensor); a
  model adds `j0` and per-`SALCKey` coefficients. Reload rebuilds the basis verbatim
  (no re-projection) and re-pairs coefficients to the basis **by key**, not by
  position (a `SCEModel` saved by one build reloads correctly into a structurally
  identical basis). The `struct ⇄ Dict` schema (`_to_doc` / `_from_doc`) is
  format-agnostic and unit-tested without any serializer; TOML is the stdlib (no
  dependency) and round-trips `Float64` exactly.
- **TOML input** (`sce/input.jl`): `read_input(path)` and the new
  `SCEBasis(path::AbstractString; backend, tol)` constructor build a basis from a
  human-authored `input.toml` (`[structure]` inline crystal, `[interaction]`,
  optional `[symmetry]`); keyword arguments override the file's backend/tol.
  Training data and the estimator are kept out of the file (loaded/chosen in Julia),
  mirroring the basis/data separation. `read_input` is exported.
- Considered JSON vs TOML deliberately: stdlib TOML round-trips `Float64` exactly and
  expresses the full nested SALC document, so the persistence artifact and the input
  file share one zero-dependency format. The schema layer stays format-agnostic, so a
  JSON (or other) backend remains a thin future addition.
- Validated (`test/unit/test_persist.jl`, `test/unit/test_input.jl`): basis / model /
  fit round-trips (predictions bit-identical, coefficients re-paired by key under a
  scrambled on-disk order, multi-op space-group ops, empty basis), input parsing +
  defaults + keyword overrides + error paths.

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
