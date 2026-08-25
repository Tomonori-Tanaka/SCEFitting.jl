# SCEFitting.jl — specification and architecture

Work-in-progress rebuild of Magesty.jl. This file tracks the **realized**
architecture as milestones land; the full design rationale lives in the
implementation plan.

## Pipeline

```
geometry → symmetry → clusters → basis (SALC) → design matrix → fit → predict
```

## Module layout (`src/`)

Single top-level `module SCEFitting`, small files included in dependency
order; the two self-contained numeric kernels (`Harmonics`, `AngularMomentum`)
are nested submodules. Heavy / optional dependencies (Spglib, GLMNet, Sunny)
live in `ext/` package extensions; the core loads without them. Persistence and
input files use the stdlib `TOML` (no external dependency).

```
src/
├── SCEFitting.jl        top-level module: includes, exports, `public` tier
├── geometry/            lattice.jl, crystal.jl, neighborlist.jl
├── symmetry/            types.jl, backend.jl (Spglib method in ext/)
├── basis/               Harmonics.jl, AngularMomentum.jl (submodules),
│                        coupledbasis.jl, salc.jl, salcbasis.jl
├── clusters/            enumerate.jl, orbits.jl
├── fitting/             estimators.jl, design.jl, fit.jl, diagnostics.jl,
│                        selection.jl (MC-cost groups, GCV, λ-path + Pareto)
├── sce/                 model.jl (pipeline types), coeftable.jl,
│                        bilinear.jl (tesseral → Cartesian extraction),
│                        introspect.jl (public introspection)
├── interop/             sunny.jl (primitive unfold + `to_sunny` entry point)
└── io/                  persist.jl, input.jl, dftsource.jl
```

The include order in the SCE layer is `sce/coeftable.jl` → `sce/bilinear.jl` →
`sce/introspect.jl` → `interop/sunny.jl`: the bilinear extraction is a core
capability consumed by both the introspection and the Sunny interop.

## Extension seams (contracts)

- **DFT sources**: `read_configs(src::AbstractDFTSource) -> Vector{<:AbstractTrainingDatum}`.
- **Estimators**: subtype `AbstractEstimator` + `solve_coefficients(est, X, y; groups)
  -> jphi` (the `(X, y)` is already column-centered and row-scaled — add no intercept,
  do not re-weight; `groups` labels rows from the same sample for grouped resampling).
  Types live in core; solver methods needing GLMNet live in `ext/`.
- **Symmetry**: `analyze_symmetry(backend::AbstractSymmetryBackend, crystal; tol) -> SpaceGroup`.
  `tol` is the backend's symprec: a **Cartesian distance in Å**, compared against
  positions and translations in Å, never in fractional coordinates.

## Realized so far

### geometry (M1)
- `Lattice(vectors; pbc)` — columns are lattice vectors; `reciprocal = inv`,
  `interplanar_spacing(lat, i) = 1/‖row_i(reciprocal)‖`.
- `Crystal(lattice, frac_positions, species, species_labels)` — fractional coords
  wrapped to `[0,1)` on periodic axes (inner constructor); `n_atoms`,
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
  primitive from `LegendrePolynomials.dnPl`. The tesseral `l ≤ 2` Cartesian-conversion
  constants `N1 = √(3/4π)`, `A2 = √(15/16π)`, `B2 = √(5/16π)` are defined once here
  (used by `sce/bilinear.jl` and downstream consumers).
- Validated by: closed-form standard solid harmonics (`l ≤ 2`), gradient tangency
  + on-sphere central difference, and agreement with Magesty's
  `TesseralHarmonics` to a stated tolerance (oracle; not bit-for-bit — the two
  use different Legendre primitives, so the last ulp differs by construction).

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
  `ext/SCEFittingSpglibExt`. `map_sym` derived in-tree (shared by all backends).

### clusters (M6) — arbitrary body order
- `ClusterMember` (atoms + per-site lattice `shift` R), `ClusterOrbit`, `ClusterSet`;
  `candidate_clusters` enumerates `N`-body **edge-admissible cliques** (`N = 2` is the
  directed neighbor pairs), `build_clusters` reduces them to symmetry orbits via a
  canonical key built from the site-image map `(b, τ + W·R)`. An edge is admissible
  under `MinimumImage` iff it sits at its atom-pair minimum-image distance (and the
  clique's atoms are distinct, so a cluster never reuses an atom's image — that would
  alias a lower-body term); under `AllImages` iff within the radial cutoff. The
  `selection` is threaded from `SCEBasis` so the neighbor list and clusters agree.
  The relative same-distance band for tie/cutoff decisions is user-facing as
  `SCEBasis(...; tie_tol)` (default `1e-8`, hard cap `1e-2`): it rides on
  `NeighborList.tol`, which the cluster-edge admission reads back, so one value
  governs both sides. Widening it is the remedy the closure refusal names for
  relaxed/noisy coordinates whose symmetry residual splits minimum-image ties.

### SALC basis (M7) — arbitrary body order
- `build_salc_basis` projects each orbit's representative onto the trivial irrep of
  its site stabilizer. At `N ≥ 3` a stabilizer op can permute equivalent sites, mixing
  coupling paths and (for unequal `l`) `l`-orderings, so the projection runs over the
  **combined (ordering × path × `Mf`) space**: each op acts by rotating every site axis
  (`wignerD_real(l_i, R)`) and relabeling axes by `invperm(perm)`, with the action
  matrix read off by contraction against the orthonormal coupled tensors. Deterministic
  axis-pivoted gauge; per-ordering fold into a multi-term SALC (one `SALCTerm` per
  ordering). `SALCKey` (canonical, injective column address — `block` runs across split
  ordering orbits) + `SALCBasis` (sorted keys + fingerprint). `evaluate_salc(salc, e)`.
  The build ends with an exact **function-space reduction per orbit**: each SALC's
  aggregated (shift-blind) monomial coefficients — keyed `(atoms, ls, index)`, the
  linear space `evaluate_salc` actually spans under cell-periodic evaluation — are
  checked, and combinations that aggregate to zero or go linearly dependent within
  the orbit (WS-boundary ties / merged near-tie shells folding distinct instances
  onto one atom set) are dropped with a warning naming orbit, channel, and reason
  (surviving keys keep their `block` numbers; gaps are legal). The independence
  guarantee is per orbit, distinct-atom members only: cross-orbit aliasing (e.g. a
  trivial space group splitting tied images into separate orbits) and row-deficient
  training data can still leave the design rank deficient — the `OLS` rank warning is
  the gate there. Repeated-atom (`AllImages` self-image) members, which the reduction
  also cannot certify, are refused outright at the `SCEDataset` door with an
  `UnclassifiableBasis`: on the reference cell both ends carry the same spin, so those
  columns are redundant by construction. Building such a basis stays legal — it is the
  tiling template a downstream consumer expands onto a supercell.
- Validated by the ground-truth tests with non-collinear spins, **all `Lf`, all body
  orders**: space-group invariance `Φ(g·e)=Φ(e)`, time-reversal evenness, linear
  independence; projector eigenvalues exactly 0/1. Improper-op parity is handled
  automatically (site-axis rotation by full `R` + even `Σl`). Cross-validated against
  Magesty: per-`(body, ls, Lf)` invariant-subspace dimensions agree through 3-body.
- A **second projection engine**, `_orbit_salcs_decors`, generalizes the same
  construction to explicit decoration labels (sorted `SiteDecor` multisets). It is
  what the pointed site-moment channel needs — that channel's mark is the
  displacement decor `SiteDecor(disp = (1, 0))`, whose factor `|u|² R₀₀` evaluates
  to 1 on the marked atom and 0 elsewhere, so it is not expressible as an `ls`
  tuple. Coupling runs over the slots spin-first, so the total spin rank `L_S` is a
  good quantum number of each coupling path, and projection is per `(L_S, Lf)`
  block; `isotropy` screens `L_S` here where the pure-spin engine screens `Lf`
  (identical on a pure-spin label, where `L_S ≡ Lf`). It is a required keyword,
  not the positional flag upstream SLCE.jl's engine takes in that slot (`soc`,
  with the opposite polarity), so a verbatim upstream call fails to compile
  instead of silently inverting the screen. A label is a sorted multiset
  and cannot express a per-site rule, so the per-species cap arrives through an
  `admit` predicate applied to an orbit of assignments. Nothing calls the engine
  from a public builder yet; the pure-spin production path is unchanged, and the
  gate is that given the same label and admission rule the two engines agree
  **bitwise** (`test/unit/test_mixedsalc.jl`).
- The joint form `evaluate_salc(salc, e, u)` evaluates spin axes as `Z_{lm}(ê)` and
  displacement axes as `|u|^{2k} R_{lm}(u)` through the 4π-free `SolidHarmonics`
  kernel, scaled by `(4π)^(n_spin/2)`. The spin-only `evaluate_salc(salc, e)`,
  `accumulate_grad!`, and `group_costs` **refuse** a decorated SALC rather than read
  a `DISP` rank as a spin harmonic under the wrong scale. Every `DISP` factor a
  pointed label carries is `(k = 1, l = 0)` (the marked site's angular rank is its
  SPIN part), so the production channel reads the kernel only at `R₀₀ ≡ 1`; the
  `l ≥ 1` path exists for the engine's own gates, which evaluate decorated SALCs to
  pin the rotation convention, scale, and slot order.

### fitting + SCE API (M8, M9)
- `BasisSpec` (validated: `nbody ≥ 1`, symmetric per-body `cutoff` matrices with
  entries `≥ 0` or `Inf`, `lsum ≥ 0` per body order, nonempty `lmax`
  with entries ≥ 0), `SCEBasis` (carries its spec in the `spec` field), `SCEDataset`
  (energy design matrix `X_E`, and the
  torque design matrix `X_T` via the four-argument form; supports `length`,
  configuration slicing `dataset[idx]` — integer/`Bool`/`:` — and `vcat` of
  same-fingerprint parts, all without recomputing design rows; the
  `SpinDatum`/source path rejects a zero moment on a basis-referenced atom),
  `SCEPredictor`/`SCEFit`
  (plus the public constructor `SCEPredictor(basis, j0, jphi)` for synthetic models —
  keys filled in from the basis),
  `fit(SCEFit, dataset, estimator; torque_weight)`, `refit(f, estimator; threshold)`
  (re-solve on the scaled-magnitude support of `f` — the de-biasing step after a sparse
  fit; shares `_assemble_problem` with `fit`), `predict_energy`/`predict_torque`,
  `coef`/`intercept`/`nobs`/`dof`/`r2_energy`/`rmse_energy`/`r2_torque`/`rmse_torque`/
  `rss_energy`/`rss_torque`/`residuals_energy`/`residuals_torque`/`has_torque` (energy and
  torque blocks reported separately; `SCEFit.residuals` stores the energy residual).
  The generics `coef`/`fit`/`nobs`/`dof`/`coeftable`/`islinear`/`predict`/`residuals`/`r2`
  **extend StatsAPI** (imported, not shadowed); `predict`/`residuals`/`r2` are thin
  wrappers defaulting to the energy block (`predict_energy`/`residuals_energy`/`r2_energy`).
  `AbstractEstimator` with `OLS`/`Ridge` (validated: `lambda` finite and ≥ 0)/`AdaptiveRidge`
  (analytic `j0`; centered-`X` `solve_coefficients` contract). `OLS` solves by
  explicit pivoted QR and **warns on a rank-deficient design** (diagonal ratio below
  `1e-10`): the solution is returned unchanged (documented min-norm behavior on
  exact deficiency) but coefficients are flagged non-unique.
- Validated: OLS recovers an in-span target (R²=1); **a Heisenberg chain fit
  recovers `J = 2√3·jphi` to rtol 1e-8** (the v0 done-line, oracle); torque is the
  exact derivative of the energy surface (on-sphere finite differences, equivariance,
  Heisenberg closed form), energy+torque co-fit recovers an in-span model.

### persistence + TOML input (M10)
- **Persistence** (`io/persist.jl`): `SCEFitting.save(path, x)` and
  `SCEFitting.load(SCEBasis | SCEPredictor, path)` serialize a self-contained,
  human-readable **TOML** document — the crystal, the space-group ops, the basis spec
  (document key `"spec"`; schema version 2 renamed it from `"interaction"`, version 3
  stores the resolved truncation — per-body `cutoff` matrices, `lsum`, labels — and
  still reads v2's scalar `pair_cutoff`; version 4 stores SALC members in the
  canonical duplicate-free form — one member per physical cluster instance, up to
  `N!`× smaller — and folds v2/v3 members on load; version 5 stores each SALC key's
  decoration multiset `"decors"` and total spin rank `"L_S"` instead of `"ls"`, and
  maps v2–v4 keys on read through the total, value-preserving relabel
  (per-site `l` → a pure-spin decor, `L_S := Lf`) — so older documents load with
  bit-identical predictions and no migration tool),
  and the *full* SALC basis (every member / term / folded tensor); a model adds `j0`
  and per-`SALCKey` coefficients. Reload reconstructs the basis verbatim (no
  re-projection) and re-pairs coefficients to the basis **by key**, not by position.
  The `struct ⇄ Dict` schema (`_to_doc` / `*_from_doc`) is format-agnostic and tested
  without any serializer; TOML is stdlib (no dep) and round-trips `Float64` exactly.
  `save` / `load` are **unexported** (call qualified) to avoid clashing with
  `FileIO`/`JLD2`.
- **TOML input** (`io/input.jl`): `read_setup(path) -> (; crystal, spec, backend, tol,
  images, tie_tol, moment)` and `SCEBasis(path::AbstractString; backend, tol, images,
  tie_tol)` build a basis from a human-authored `input.toml` (`[structure]` inline
  crystal, `[interaction]` with optional `images` (`"minimum_image"` default /
  `"all_images"`), `tie_tol`, and `cutoff = inf` for the full WS cell, optional
  `[symmetry]`); keyword arguments override the file's backend/tol. The optional
  `[moment]` section is read into a `MomentSpec` (values handed to the keyword
  constructor, which owns all validation; unknown keys and the upstream `soc` spelling
  are refused) and `MomentBasis(path; backend, tol, tie_tol)` builds the pointed basis
  from it, sharing `[interaction].tie_tol`. Training data and the estimator stay in
  Julia (mirrors the basis/data separation).
- **Pointed moment basis** (`basis/momentbasis.jl`): `MomentSpec` / `MomentBasis` /
  `moment_resolvability`, the adiabatic site-moment channel's counterpart of
  `SCEBasis`. The mark is `SiteDecor(disp = (1, 0))` riding the decor engine through
  its `admit` kwarg; mark-aware admission (mark `ê` factor for any species up to
  `lmax_mark`, environment spins only for `sampled` species up to `lmax_env` —
  refused loudly when a spin-bearing species is unsampled); 1/2-body clusters from
  the ordinary enumeration at `cutoff_pair`, stars (`3 ≤ N ≤ nbody`, cap 4) cut on
  their `N−1` mark–environment spokes at `cutoff_star` — stored **per star order**
  (`cutoff_star[N-2]`, read through `_star_cutoff`; a scalar or matrix broadcasts to
  every order, and the one shared star neighbor list is built at the elementwise
  envelope) — with env–env edges free, with
  `candidate_clusters`' all-orderings
  multiplicity convention (`_orbits_from_members` is `build_clusters`' factored
  core); even-Σl (TR) only, capped by `lsum` **per body order** (`lsum[N]`, indexed by
  the order itself and read through `_label_lsum`; a scalar broadcasts, body-keyed
  pairs name orders, unnamed orders stay uncapped — `BasisSpec.lsum`'s spelling and
  resolver, `_resolve_lsum`); `isotropy = true` keeps `L_S = 0` (upstream's
  `soc = false`). `nbody ∈ 1:4` — the enumeration is general in `N`, the cap is where
  the oracles stop; a body order that no label can reach (`Σl = 2⌈(N−1)/2⌉` is its
  floor) is dropped with a warning and `show` reports the realized order.
  `_design_moment`: rows `(config, marked atom)`, marked-column
  substitution (exact — one mark per label), threaded per column, a mark→term index
  value-identical to the full evaluation. `moment_resolvability` (D9′): symbolic
  signature rank in the independent variables `(a, ê_a, e)`, null combinations naming
  dependent columns, mark-class census; refuses repeated-image environment members
  (`UnclassifiableBasis`) rather than overcounting. Validated by independent
  references (`test_momentbasis.jl`): FeGe B20 star closed form = 6.0 × a geometric
  triangle enumeration, shell-sum normalization 2√3, G_i covariance with axes,
  bitwise TR, substitution locality, symbolic ≡ random-design rank. The 1-body
  `[MARK]` columns are the per-orbit intercepts μ₀. The star enumeration keeps a
  star whose environment sites include two minimum images of ONE neighbor (as
  upstream does, for column parity on small cells); such members reduce to
  lower-body functions, and the resolvability gate refuses the basis as
  `UnclassifiableBasis` — a hard door at `MomentDataset`, not an option.
- **Moment dataset / fit / model** (`fitting/momentfit.jl`; included after the io
  layer — its constructor takes `SpinDatum`s): rows `(config, marked atom)`,
  `y = ê·M` under the mode rule (mixed modes allowed; a mode-1 zero-axis marked atom
  is `defined = false`, `y = NaN`, excluded from every fit and recorded). Doors:
  `moments_bare` + `constraint_mode` required per datum; every `SpinDatum` sits at
  the reference geometry by construction. Decomposability gate
  `g = ‖M⊥‖²/|M| = |M| sin²θ ≤ gate_eps` (required keyword; cancellation-free
  `M⊥ = M − y ê` form; `|M| = 0` passes; no `m_min`); per marked-atom orbit
  (`map_sym` min-atom) survival + rms `‖M⊥‖` + antiparallel census `n_anti` reported
  and the `coverage_floor` refusal BEFORE the design build; `moment_resolvability`
  at the same door (unclassifiable → refuse, vanishing columns recorded, dependent
  combinations warned + stored). `X` built for all rows with `defined`/`keep` masks;
  `fit(MomentFit, ds, est = OLS())` solves gated and (for disclosure) ungated rows
  with `groups = row_config`, freezes vanishing columns to exact zero, and does no
  centering — the `l = 0` μ₀ `[MARK]` columns are the per-orbit intercepts, and
  regularizers penalize them like any column (v1 OLS-first). `MomentModel` carries
  basis + gated coefficients (no mode flag); `predict_moment(model, e; axes = e)`
  defaults to the mode-4 identity and is a validating door (`e` unit everywhere,
  `axes` unit + component-bounded on marked columns only). The dataset constructor
  is the training-side door: every datum's `directions` pass `_validate_config`
  (a field-built `SpinDatum` carries no direction check), and a referenced atom
  (marked or environment) with `‖MW‖ ≤ zero_moment_atol` is refused — its
  direction is the ẑ placeholder. **Not persisted** (design record §4.2): `save`
  refuses the moment types by name. Diagnostics on the same rows: `order` (per
  config `|⟨e⟩|`) + `moment_band_profile` (equal-count bins, bin-free line, `r`);
  `moment_local_field` (`‖h₁‖`, `ê·ĥ` over the `cutoff_pair` MinimumImage
  neighbors at the basis's `tie_tol`, axis by the mode rule through the single
  `_moment_axis_matrix`) + `moment_coverage` (upper-tail `h1` quantile,
  `frac_anti`); `moment_simple_floor` (per-orbit intercept + Legendre shell sums
  on the kept rows, `inclusion` report, `nested_bound` only for OLS, bitwise
  pairing door). `salc_groups(::MomentBasis)` keys on the mark class (marked
  atoms + marked sites of the canonical member) and `GroupAdaptiveRidge(mb;
  lambda)` carries unit weights; `fit` reduces it to the active columns with the
  vanishing freeze (`_reduce_to_active`).
- **Tabular results** (`sce/coeftable.jl`): `coeftable(fit | model) -> SCECoefficients`
  is a **Tables.jl** source — one row per SALC (`body`, `orbit_id`, `decors` as a
  comma string — `"1,1,2"` on a pure-spin key, exactly the old `ls` column —, `L_S`,
  `Lf`, `block`, `J`) — so it drops into `DataFrame` / `CSV.write` /
  `Arrow.write`. The library owns the internal-storage → labeled-row mapping; the caller
  brings the table/IO package. `j0` is the intercept (`intercept(c)`), not a row.
  Tables.jl is a lightweight core dep.
- **DFT data sources** (`io/dftsource.jl`): the **code-agnostic boundary only** —
  `SpinDatum` (energy + spin directions + magmoms + constraining field + the derived
  torque target `τ_a = m_a × B_a`, the physical / Landau–Lifshitz torque, plus the
  optional adiabatic-moment trio `moments_bare` / `constraint_axes` /
  `constraint_mode ∈ {1, 4}` — mode 1 requires the axes, axes require a mode, the
  bare moment is validated for finiteness only; the E/T paths never read them) and
  `read_configs(src::AbstractDFTSource) -> Vector{SpinDatum}`, with `SCEDataset(basis, src)`
  going source → dataset. Two in-core formats: the **extended-XYZ training container**
  (`io/extxyz.jl`: `write_extxyz` / `read_extxyz` / `ExtxyzFile`; the structure is
  always stored, per-atom columns 1:1 with the datum's channels, shortest-round-trip
  numbers, the same dialect as SLCE.jl — spin-only files interchange, a joint file is
  refused by name, never flattened) and Magesty's legacy EMBSET set (`read_embset`,
  plus `read_embset_pair` for `EMBSET` + `EMBSET_mint` sibling archives: config count /
  block shape / field blocks bitwise, energies deliberately uncompared). The moment
  channel's axis gates (`check_moment_gates`: mode-1 sign consistency on decomposable
  rows + axis-angle p99) run at extxyz generation, extxyz load and the pair reader.
  The **concrete per-code adapters live in the `SCETools.jl`
  package** (`SCETools.VASP`: `read_poscar`/`write_poscar`, `Oszicar` with SAXIS rotation,
  and the INCAR writer), not in the core. Adding a DFT code is one sibling adapter there —
  neither the core nor its export list changes; the VASP parsers are cross-checked bit-for-bit
  against Magesty in SCETools's oracle. The one in-core concrete format is Magesty's
  code-agnostic EMBSET training set (`io/embset.jl`): `read_embset` + the `EmbsetFile`
  source, cross-checked against `Magesty.read_embset` in this package's oracle.
- Validated: basis / model / fit round-trips (predictions bit-identical, coefficients
  re-paired by key under scrambled order, multi-op space-group ops, empty basis), input
  parsing + defaults + keyword overrides + error paths.

### Sunny.jl export (M11)
- **Bilinear extraction core** (`sce/bilinear.jl`, dependency-free): `_l1_pair_matrix` /
  `_l2_onsite_matrix` turn a folded tesseral tensor into a Cartesian exchange / single-ion
  matrix (`eₐ'·M·e_b = Σ folded·Z·Z`, via the `Harmonics.N1`/`A2`/`B2` tesseral constants);
  `_classify_salc` keeps only `ls=[1,1]` pairs and `ls=[2]`
  single-ion (`ls=[0…]` → `j0`, the rest skipped + reported). `_bilinear_terms`
  folds the directed members into one `BilinearTerms` matrix per undirected supercell bond.
  `_reconstruct_energy ≈ predict_energy − j0` gates the whole chain **without Sunny**.
- **Sunny-specific core** (`interop/sunny.jl`, dependency-free): `_sunny_primitive` unfolds
  the supercell matrices onto the chemical primitive cell recovered from the
  pure translations (one Sunny bond per primitive bond, a `clean` flag for the fallback),
  plus the `to_sunny` entry point (a friendly "load Sunny" error without the extension).
- **Extension** (`ext/SCEFittingSunnyExt`, loaded by `using Sunny`):
  `to_sunny(model; spins, g = 2, mode = :auto, scaling = :auto, placement = :auto)
  -> Sunny.System` builds a real `System`
  (`set_exchange_at!`/`set_onsite_coupling_at!` on the supercell, or
  `set_exchange!`/`Bond` on the primitive cell); its classical
  energy reproduces `predict_energy − j0`. `scaling` chooses how the effective spin
  `S_eff` enters: `:moment` puts `S_eff` into Sunny's `Moment` and rescales
  `J = M/(SₐS_b)` (energy and dispersion exact; half-integer `S_eff` only), `:coupling`
  keeps a placeholder `Moment` (`s₀ = 1`) and folds `S_eff` into the couplings (any
  `S_eff`; dispersion-exact, static energy rescaled), `:auto` picks `:moment` for
  half-integer spins, else `:coupling`. `placement = :explicit` (exact, folded
  dispersion) / `:primitive` (unfolded) / `:auto`.
- Validated: Sunny-free conversion + primitive-fold tests (`test/unit/test_sunny.jl`,
  main suite) and a separate `test/sunny/` environment (real `Sunny.System` energy vs SCE
  for both routes, the primitive system reshaped back to the supercell, mode/spin handling,
  the skip warning).

### Fitted-model introspection (M12)
- **Public, stable contract** (`sce/introspect.jl`): `multipole_terms(model) ->
  Vector{MultipoleTerm}` is the code-neutral view downstream packages (the `SCETools.jl`
  mean-field samplers) read **instead of** the SALC-basis internals (`SALCMember` /
  `SALCTerm`). Each `MultipoleTerm` carries the raw fitted `coef = jϕ`, the cluster
  `atoms` / `shifts`, the per-site `ls`, and the `folded` tensor; the per-`N` scale
  `(4π)^(body/2)` is **left to the consumer** (applied in exactly one place, the
  `_energy_from_terms` reconstruction gate). `bilinear_terms(model)` is the thin public
  wrapper of the general Cartesian bilinear / single-ion extraction `_bilinear_terms`
  (`sce/bilinear.jl`).
- Validated: `test/unit/test_introspect.jl` (the terms reconstruct `predict_energy − j0`).

### Regularized estimators (M12)
- **Analytic, in-core** (`fitting/estimators.jl`, no extension): `AdaptiveRidge(; lambda,
  epsilon, max_iter, tol)` — iterative reweighted ridge (Frommlet & Nuel 2016), an L0
  approximation that refits `(X'X + λ·Diagonal(w)) \ X'y` with `wⱼ = 1/(βⱼ² + ε)` until the
  ∞-norm change drops below `tol`; `lambda = 0` ⇒ OLS. `islinear` ⇒ `true` (a linear
  smoother in the converged-weight sense). `PrecomputedPilot(beta)` — adapter returning a
  fixed coefficient vector (length-checked against `size(X, 2)`), for reuse as an
  `AdaptiveLasso` pilot.
- **GLMNet-backed types in core, solve in extension** (`ext/SCEFittingGLMNetExt`,
  loaded by `using GLMNet`): `ElasticNet(; alpha, lambda, standardize, nfolds, select,
  seed, nlambda)`, `Lasso(; …)` (= `alpha = 1`), and `AdaptiveLasso(; pilot, lambda, gamma,
  epsilon, standardize, …)`. Named and dispatched on without the dependency; argument
  validation lives in core. GLMNet minimizes `(1/2n)‖y−Xβ‖² + λ[(1−α)/2‖β‖₂² + α‖β‖₁]` on
  the column-centered `(X, y)` with `intercept = false` (so `j0` stays analytic) and
  `standardize` (penalty acts per-column at `λ·std`). `lambda = nothing` ⇒ K-fold CV over
  GLMNet's λ path, `select` = `:lambda_min`/`:lambda_1se`; folds are **grouped by
  configuration** (via `groups`) and seeded deterministically (`hash`-ranked, no RNG
  dependency). A numeric `lambda` skips CV.
- `AdaptiveLasso` (Zou 2006) runs `pilot` (any estimator; default `OLS`) for `β̂`, then a
  weighted Lasso with per-column penalty factor `wⱼ = 1/max(|β̂ⱼ|, ε)^γ` (held fixed across
  the fixed-λ or CV path). `gamma = 0` reduces exactly to a plain Lasso. The pilot's own
  solve receives the co-fit `groups`.
- Validated in a separate `test/glmnet/` environment: tiny-λ ≈ OLS, the analytic-`j0`
  centering invariant, CV support recovery + sparsity, `:lambda_1se` shrinkage,
  reproducibility, ElasticNet ≠ Lasso, AdaptiveLasso support recovery, `γ = 0` ≡ plain
  Lasso, pilot pluggability (`Ridge`/`PrecomputedPilot`), and energy+torque grouped-CV
  co-fits. Core-only construction / validation / deferred-backend error, plus the full
  `AdaptiveRidge` / `PrecomputedPilot` solves, in `test/unit/test_fit.jl`.

### Cost-weighted group selection (M13)
- **Estimator** (`fitting/estimators.jl`, in-core): `GroupAdaptiveRidge(column_groups,
  group_weights; lambda, epsilon, max_iter, tol)` — the group extension of
  `AdaptiveRidge`, approximating the weighted group-L0 `λ·Σ_g v_g·1{β_g≠0}` by iterating
  `Dⱼ = mⱼ·v_g/(Σ_{k∈g} m_kβ_k² + p_g·ε)` (all columns of a group share one weight; a
  surviving group's converged penalty is exactly `λ·v_g`). Exact degeneration to
  `AdaptiveRidge` for singleton groups with unit weights and a uniform metric;
  `lambda = 0` ⇒ OLS; `islinear` ⇒ `true`.
  `column_groups` labels **columns** (contiguous `1:G`, validated in the inner
  constructor) — unrelated to the per-row `groups` kwarg of `solve_coefficients`. The
  weight map `_gar_weights!` is the single definition of the group form, shared with the
  GCV diagnostics.
- **Penalty metric** (`fitting/selection.jl` / `fitting/momentfit.jl`, exported):
  `penalty_metric(basis; torque_weight, nconfig, seed)` = `(1−w)·Var[Φⱼ] +
  w·E[Σ_a‖(∂Φⱼ/∂e_a)×e_a‖²]/(3·n_atoms)` and `penalty_metric(mb; free_intercepts,
  nconfig, seed)` = `E[Φⱼ²]`, both over uniform-random reference configurations drawn by
  an in-package generator (SplitMix64 + the Archimedes construction; `Random`'s stream
  carries no cross-version guarantee and this quantity enters every penalized
  coefficient). `nconfig = 2048` by default, sized from the measured relative standard
  error (3.3 % median / 5.3 % worst column; `1/√nconfig`) against the measured
  construction cost (0.56 s on bcc Fe 3×3×3, linear in `nconfig` and column count). Carried on `Ridge` /
  `AdaptiveRidge` / `GroupAdaptiveRidge` as `metric` (`nothing` = uniform) with a
  validating `MetricProvenance` (channel, `torque_weight`, `nconfig`, `seed`, basis
  fingerprint) the `fit` / `refit` / `select_fit` / `cross_validate` doors check;
  `AdaptiveLasso` carries one through its pilot. `SCEFitting.with_lambda(est, λ)` moves
  an estimator along a λ path with the metric intact.
  The metric sits in the **denominator** of the adaptive weight maps, which is what makes
  the estimators invariant under a column rescaling and preserves the `λ·v_g` fixed
  point; the IRLS cold starts and the stopping rule (metric coordinates, penalized
  columns only) move with it. `mⱼ = 0` ⇒ column unpenalized, reserved for structural
  exemptions (the moment channel's μ₀ intercepts, identically vanishing columns) and
  refused when it comes out of a sample; the unpenalized block must be well conditioned
  or the solve is refused. Basis-aware constructors attach it by default.
- **Basis helpers** (`fitting/selection.jl`; public, unexported): `salc_groups(basis)`
  — column → group labels by `(body, orbit_id, decors)`, the granularity at which MC
  contraction entries vanish; `group_costs(basis, labels)` — per-group distinct-entry
  union count over canonical members (additive across the `salc_groups` partition);
  `cost_weights(basis; theta)` — `v_g = √p_g·(c_g/c̄)^θ`, `θ ∈ [0, 1]` tilting the
  penalty from cost-blind to cost-proportional. Convenience constructor
  `GroupAdaptiveRidge(basis; lambda, theta, …)` bundles them.
- **GCV / effective dof** (exported): `effective_dof(f)` = `tr(X(X'X+λD)⁻¹X') + 1` with
  the converged penalty diagonal recomputed from the fitted coefficients
  (`_penalty_diagonal`, one method per linear estimator); `gcv(f)` = `n·RSS/(n−df)²` on
  the assembled problem, `Inf` in the near-interpolating regime `df → n`. The dof trace
  is an eigenproblem on the smaller Gram side (`p ≤ n`: weighted `X'X`, reusing the
  path's cached Gram; `n < p`: the `n×n` dual) — never an `n×p` SVD. With unpenalized
  columns (`Dⱼ = 0`) it is `rank(X_F) + Σᵢ sᵢ/(sᵢ+λ)` over the eigenvalues of
  `W^{-1/2}X_P'(I−P_F)X_P W^{-1/2}` (thin QR, never an `n×n` projector), taken before
  the cached-Gram branch and refusing a rank-deficient `X_F`. Linear estimators
  only; on torque co-fits GCV is optimistic (correlated within-configuration rows) —
  grouped CV is the ground truth there (documented, not an error).
  `effective_dof(::MomentFit)` / `gcv(::MomentFit)` are the pointed counterparts, built
  on the design the fit SOLVED (gate-kept rows, frozen columns removed, estimator
  reduced) and with no `+1` (the μ₀ intercepts are columns of that design).
- **λ path + Pareto** (exported): `select_fit(dataset, est; lambdas, torque_weight,
  criterion = :gcv|:cv, delta, costs, threshold, nfolds, seed) -> SelectionPath` —
  descending warm-started path on a once-assembled Gram; per-λ score, effective dof,
  alive groups (the `refit` scaled-magnitude rule, any-column-per-group) and predicted
  MC cost `Σ_{g alive} c_g`; selection = cheapest λ within `(1+δ)` of the minimum score
  (the cost-aware generalization of `:lambda_1se`; `Inf` scores never eligible, cost
  ties → larger λ). The adaptive iteration never yields exact zeros, so the default
  `threshold = nothing` is a per-λ **relative** alive floor (`_ALIVE_RTOL = 1e-6` of
  that λ's largest scaled magnitude; an absolute number reproduces `refit`'s rule,
  `0.0` degenerates to all-alive). `criterion = :cv` is configuration-grouped K-fold in
  core (deterministic seeded folds; per-fold Gram downdate; fold reduction warns). The
  selected fit is re-solved cold — its path row (`n_alive`/`cost`) is re-derived from
  that cold solve and the effective absolute threshold is returned as
  `path.threshold`, so `refit(path.fit; threshold = path.threshold)` realizes exactly
  the reported support; `SelectionPath` is a Tables.jl source.
- **Threshold front** (exported): `select_support(f; npoints = 25, thresholds =
  nothing, delta, labels, costs, evalset = f.dataset, estimator = OLS()) ->
  SupportPath` — the second knob:
  sweeps the alive threshold at a fixed fit (auto grid = log-rank-spaced points on
  the per-group scaled-magnitude spectrum + the full-support anchor, or an explicit
  vector), de-biases with `refit` per point, scores each refit by the fit's own
  `(1−w)·MSE_E + w·MSE_T` objective on `evalset` (pass a held-out slice for an honest
  axis; fingerprint-checked against the training basis), and applies the same
  Pareto rule. Needed because real-data group-magnitude spectra are continuous (no
  alive/dead gap for the λ path to expose). `SupportPath` is a Tables.jl source
  (`threshold`/`n_alive`/`cost`/`score`/`rmse_energy`/`rmse_torque`/`selected`).
- **Generic CV** (exported): `cross_validate(dataset, estimator; torque_weight,
  nfolds = 5, seed = 1) -> CVResult` — configuration-grouped K-fold assessment of any
  `fit` call: each fold refits from scratch (fold-local centering/whitening, no
  leakage) and scores the held-out configurations in prediction space. Reports the
  per-fold and pooled out-of-fold energy **and** torque RMSEs independently of
  `torque_weight`, plus the `(1−w)·MSE_E + w·MSE_T` score. Deterministic seeded
  folds (`_grouped_folds`); fold reduction warns, `< 6` configs errors; a
  `PrecomputedPilot` is rejected (fold-independent coefficients would leak).
  `CVResult` is a Tables.jl source. Unlike `select_fit(criterion = :cv)` (global
  whitening, λ ranking only), this is the honest generalization-error estimate.
- **Moment CV** (exported): `cross_validate(ds::MomentDataset, estimator; nfolds, seed)
  -> MomentCVResult` — the pointed channel's λ criterion. Folds are grouped by
  configuration (its per-marked-atom rows never split), each fold re-solves on its
  training rows with the SAME frozen column set as the full dataset, and a fold whose
  training rows miss a marked orbit is refused by name (that orbit's μ₀ would be
  unidentified). Training and scoring run on the gate-kept rows; `score_defined` reports
  the gate-rejected rows separately as disclosure. One error axis (the `ê·M` RMSE in
  μ_B), hence its own result type; `pooled_score` aggregates out-of-fold residuals.
  Tables.jl source; `PrecomputedPilot` rejected as above.
- Validated in `test/unit/test_selection.jl`: construction/validation, exact
  `AdaptiveRidge` degeneration, group-sparse recovery + weight monotonicity, label/cost
  hand counts + additivity, `θ` endpoints, dense-hat-matrix trace agreement, the
  underdetermined regime + guards, warm/cold path consistency, the Pareto rule, and the
  end-to-end select → refit workflow.

## Not yet implemented (v0 follow-ups)
- The v0 slice is feature-complete; no estimator/observable/IO follow-ups outstanding.

## Oracle environment (`test/oracle/`)

A separate Julia env (`[sources]`-deving both `SCEFitting` and a pinned
`Magesty.jl`) cross-checks convention-fixed kernels and gauge-invariant
aggregates / predictions against Magesty. The core suite never depends on
Magesty. Run: `julia --project=test/oracle test/oracle/runtests.jl`.

(Sections for M3–M10 added as they land.)
