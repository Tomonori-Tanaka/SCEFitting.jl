# CLAUDE.md

> Shared baseline (numerical-correctness priority, JP-conversation / EN-repo
> policy, Conventional Commits, Julia style, shared subagents) is inherited from
> `~/Packages/CLAUDE.md`. Only package-specific rules live here.

## Project goal

Clean, extensible, Julia-native rebuild of `Magesty.jl`: fit a spin-cluster
expansion (SCE) `E = j0 + Σφ Jφ Φφ({e_a})` to noncollinear DFT data. The numerical
core (tesseral harmonics, Clebsch–Gordan coupling, symmetry-adapted basis, design
matrix, regression) is **reimplemented from scratch**; `Magesty.jl` is used only
as a *pinned numerical oracle* in `test/oracle/`. Priority: numerical correctness
and reproducibility over stylistic concerns. See `SPEC.md` for the realized
architecture and `docs/design-notes.md` for the rationale.

Both SCE observables are fitted: the **energy** `E` and the per-atom **torque**
`τ_a = −e_a × ∂E/∂e_a` (the physical / Landau–Lifshitz torque `m_a × B_eff,a`, the
analytic derivative of the same energy surface). An
energy+torque co-fit minimizes `L = (1−w)·MSE_E + w·MSE_T` for a `torque_weight`
`w ∈ [0,1]`. Cluster enumeration is **arbitrary body order** (`nbody = K`):
pairwise-within-cutoff cliques, with the SALC projection generalized to the combined
(ordering × coupling-path × `Mf`) space so that permutation-equivalent sites — which
at `N ≥ 3` mix coupling paths and, for unequal `l`, `l`-orderings — are handled.

Bit-for-bit agreement with Magesty is **not** a goal — refining methods so results
differ slightly is allowed and is the point of the exercise. Validation is
physical/intrinsic (symmetry invariance, finite differences, known-coupling
recovery); the oracle compares convention-fixed kernels and gauge-invariant
aggregates, never raw gauge-dependent SALC coefficients.

## Numerical / physics conventions

Easy to break silently — confirm before touching the algorithm.

- **Spin directions are unit vectors**; `directions[:, a]` (norm 1) is the
  direction, the moment magnitude `magmom[a]` is stored separately.
- **Spin layout `3 × n_atoms`** (rows x, y, z; columns atoms). Transposing breaks
  the pipeline.
- **Real (tesseral) spherical harmonics `Zₗₘ`** (Drautz, PRB 102 024104), per-site
  factor `(4π)^(−1/2)`; the design matrix carries `(4π)^(N/2)` so an N-body term
  cancels to O(1).
- **Time reversal**: even-`Σl_s` channel only (odd-`Σl` `ls`-assignments are not
  enumerated). The SALC projector rotates each **site** axis by `wignerD_real(l_i, R)`
  (full `R`, improper included) and reads the `Lf` action off by contraction against
  the orthonormal coupled tensors; the per-site `(−1)^{l_i}` parity with even `Σl`
  makes the improper handling automatic (no explicit proper-part case), so
  symmetry-forbidden odd-`Lf` channels are dropped and allowed ones kept.
- **Lattice / reciprocal**: `Lattice.vectors` columns are `aᵢ`; `reciprocal =
  inv(vectors)` rows are `bᵢ` with `aᵢ·bⱼ = δᵢⱼ` (no 2π). Interplanar spacing
  `dᵢ = 1/‖row_i(reciprocal)‖`. Fractional coords are wrapped to `[0,1)` on
  periodic axes (neighbor-list precondition).
- **Periodic resolvability (minimum image / Wigner–Seitz)**: a finite supercell under
  plain PBC can only resolve interactions whose displacement lies in the **Wigner–Seitz
  cell** of the (super)lattice — a polyhedron reaching the body-diagonal corner
  `(L/2,L/2,L/2)` at `√3·L/2`, **not** a sphere of radius `L/2`. A farther periodic image
  of an atom carries the *same spin*, so its interaction is collinear with (an alias of)
  the minimum-image one and is not independently fittable. The default `MinimumImage`
  selection enumerates exactly this set (boundary ties kept; `i==j` self-pairs and
  reused-atom clusters dropped); `cutoff = Inf` is the whole WS cell. `AllImages`
  (every image, `R`-distinguished) is **only** for the future generalized-Bloch /
  spin-spiral path where `e^{iq·R}` resolves the images. For `cutoff < min_d dᵢ / 2` the
  two coincide. Do **not** "fix" a `> L/2` cutoff by folding aliases into a shorter shell
  (double-counts) — that regime is simply unresolvable from one supercell.
- **Energy units**: `Jφ` carry the DFT input unit (eV); `j0` is separate.
- **Torque**: `τ_a = −e_a × ∂E/∂e_a` (the physical / Landau–Lifshitz torque `m_a × B_eff,a`,
  matching the *General spin models* paper; the opposite sign of the energy-rotation-gradient
  `+e×∇E`), design-matrix entry `(4π)^(N/2)·(−e_a × ∂Φ/∂e_a)` = `(4π)^(N/2)·(∂Φ/∂e_a × e_a)`
  (same scale and `μ`-mapping as the energy kernel). The torque design matrix `X_T`
  has no `j0` column. Its rows are flattened config-major, then atom-major, then
  `xyz`. The co-fit whitens the (centered) energy block by `√((1−w)/n_E)` and the
  torque block by `√(w/n_T)`; `j0` stays an energy-only quantity.

## Coupled ("linked") code sites — change one, check all

- `basis/Harmonics.jl` (`Zlm`, `grad_Zlm`) ↔ the on-sphere central-difference and
  closed-form agreement tests (`test/unit/test_harmonics.jl`). Normalization / sign
  drift silently biases `X`.
- **Energy kernel `evaluate_salc` ↔ gradient kernel `accumulate_grad!`** (`basis/salc.jl`):
  identical `μ = idx[i] − ls[i] − 1` mapping, `ls`, `folded`, and `(4π)^(N/2)` scale.
  The gate is the finite-difference self-consistency `predict_torque ≈ −e × ∇E_FD`
  (`test/unit/test_torque.jl`, `test_nbody.jl`): the torque must be the exact (negative
  rotation-) derivative of the energy surface. Change one kernel, re-check the other.
  Both spin-only forms **refuse** displacement-decorated SALCs (the joint energy form
  is `evaluate_salc(salc, e, u)`; a joint gradient does not exist in this pure-spin
  package, and the refusal is the guard against silently reading a `DISP` rank as a
  spin harmonic under the wrong `(4π)` scale). `group_costs` refuses for the same
  reason.
- **Pure-spin ↔ decor SALC engines** (`basis/salcbasis.jl`):
  `_project_and_fold`/`_transport_term`/`_enumerate_ls` (production, oracle- and
  pin-checked bitwise) and `_project_and_fold_decors`/`_transport_term_decors`/
  `_orbit_salcs_decors` (mixed channels — the pointed site-moment channel's engine,
  since its mark is a displacement decor) are two implementations of ONE
  construction. The decor engine must reproduce the pure-spin engine's enumeration
  order exactly — orbit representatives discovered in colex (`Iterators.product`)
  order, lex-min representative, assignments `unique(rep[p])` — or `block` indices
  and the canonical gauge silently change and key-addressed coefficient re-pairing
  breaks. A label is a sorted multiset, so a per-site rule (the pure-spin engine's
  per-species `lmax`) has to be handed in through `admit`, whose verdict must be a
  permutation-orbit invariant. `isotropy` screens `Lf` in the pure-spin engine and
  `L_S` in the decor one; those coincide only on pure-spin labels. In the decor
  engine it is a **required keyword**: upstream's engine takes `soc::Bool` in that
  positional slot with the OPPOSITE polarity (`soc = true` keeps every `L_S`), so a
  positional copy of an upstream call must be a `MethodError`, never a silently
  inverted screen. Keep it a keyword when porting. Gate: "engines
  agree on pure spin" in `test/unit/test_mixedsalc.jl`, incl. the Cs-triangle
  (2 ordering orbits) and C3v-triangle (3 assignments) shapes. Change either engine
  → re-run that gate + the oracle suite.
- **`SALCKey`/`SALC` field surface ↔ ALL test environments**: the unit suite is not
  the only consumer — `test/sunny/runtests.jl`, `test/glmnet/runtests.jl`,
  `test/oracle/runtests.jl`, and `examples/*.jl` read key/SALC fields and are NOT
  exercised by the core suite (`TEST_MODE=all`). Rename a field → grep all of
  `test/` and `examples/`. [Backported from SLCE.jl 199d4eb, where the decor-label
  rename missed sunny/glmnet/examples until review.]
- **Image selection ↔ neighbor list ↔ cluster edges** (`geometry/neighborlist.jl`,
  `clusters/enumerate.jl`, `sce/model.jl`): `SCEBasis` threads one `images` value to
  **both** `build_neighbor_list` and `candidate_clusters`/`build_clusters`; they must
  agree. `MinimumImage` keeps minimum-image pairs (no `i==j`) and admits a clique edge
  only at its atom-pair minimum-image distance with all atoms distinct; `AllImages` keeps
  every in-cutoff image and admits edges within the radial cutoff. The tie/cutoff
  tolerance is relative (`_SAME_DIST_RTOL`) on both sides so a degenerate WS-boundary
  shell is never split; it is user-facing as `SCEBasis(...; tie_tol)` (default the
  same constant, hard cap `_TIE_TOL_MAX = 1e-2`), riding on `NeighborList.tol` which
  `candidate_clusters` reads back — one value, both sides. Widening it is the remedy
  for relaxed/noisy coordinates whose symmetry residual splits ties (the orbit
  builder's closure refusal names it); the SALC reduction then handles the merged
  shell's aggregated redundancy exactly (gate: the perturbed-honeycomb closure
  testset in `test/unit/test_clusters.jl` both refuses at the default band and
  builds the ideal orbit structure at `tie = 1e-3`). The minimum-image search box is adaptive — change it and re-check
  the skewed-cell test. `images` is **not** persisted (the full SALC basis is stored and
  reloaded verbatim), so only `read_setup`/`SCEBasis` carry it — and the same rule
  holds for `tie_tol` (`[interaction].tie_tol` in the TOML schema): both change the
  emitted basis, so both must round-trip through the setup file. **At `N ≥ 3` the clique
  check is on the actual chosen images of *every* pair (not just the anchor edges): a
  cluster is admitted only when all `C(N,2)` edges sit at their atom-pair minimum image
  simultaneously** (the compact-cluster criterion). Having each pair individually
  minimum-image-resolvable is *not* sufficient — the images minimizing `i–j` and `i–k`
  may force `j–k` onto a longer image, which must reject the cluster. WS-boundary ties
  then multiply N-body clusters just as they do pairs. The whole count (candidates and
  symmetry orbits) is pinned against an independent brute force in
  `test/unit/test_ws_nbody.jl`; change the edge rule or the search box and re-check it.
  **Per-body × per-species-pair cutoffs** layer on top: the neighbor list is built at
  the element-wise max over body orders (each pair admitted against its own
  species-pair radius, tie band per pair), and `candidate_clusters` re-checks every
  edge against *that body order's* radius — so `neighbors` is a superset that the
  per-order matrices trim. The per-pair admission has its own brute-force gate in
  `test/unit/test_truncation.jl`; change either side and re-check both pins.
- **SALC construction ↔ the ground-truth invariance test** (`test/unit/test_salc.jl`,
  `test/unit/test_nbody.jl`, `test/oracle/runtests.jl`): every SALC must satisfy
  `Φ(g·e) = Φ(e)` (non-collinear spins, all `Lf`, **all body orders**) and
  `Φ(−e) = Φ(e)`. The combined-space projection action (`_project_and_fold`) and the
  member transport (`_transport_term`) must use the *same* rotation direction and
  axis-relabel convention (`invperm(perm)`); the eigenvalue-exactly-0/1 idempotency
  assertion and the invariance test are the gates. **The construction's last step is
  `_canonicalize_members`** (`basis/salc.jl`): the ordered, anchored images (one per
  site ordering — the space the projection needs) are folded to one member per
  physical instance (sites sorted by `(atom, shift)`, `shifts[1] == 0` re-anchored,
  tensors `permutedims`-aligned and summed per site→`l` assignment — an exact
  regrouping). Everything downstream (evaluate/gradient kernels, introspection,
  persistence, analytic test references like the oracle's Heisenberg
  `Φ = 2√3 Σ_undirected`) assumes the canonical form: one member per undirected
  instance, tensors carrying the whole (formerly per-image) weight. Gates:
  `check_canonical_members` / `split_roundtrip_exact` (shared helpers in
  `test/unit/testutils.jl`, applied by `test_salc.jl` and 3-body `test_nbody.jl`)
  and the pre-v4 fold test in `test_persist.jl`. At `N ≥ 3` also confirm SALCs are
  linearly independent (design-matrix rank = #SALC). **After canonicalization each
  orbit is function-space reduced** (`_reduce_orbit_salcs`, `basis/salcbasis.jl`):
  SALCs are expanded into aggregated shift-blind monomial coefficients keyed
  `(atoms, ls, index)` — the space `evaluate_salc` actually spans, since evaluation
  never reads `shifts` — and combinations that aggregate to zero or go dependent
  within the orbit (WS-boundary ties / merged near-tie shells folding instances onto
  one atom set) are dropped with a warning (bulk MnTe + SOC emitted 51 with rank 37
  before this; silent OLS `max|coef| ~1e7`). The independence guarantee is **per
  orbit, distinct-atom members only** — there the vector algebra is *equivalent* to
  function algebra (orthogonal monomials; precondition: cluster sites carry `l ≥ 1`,
  `_enumerate_ls`, so the atom set is recoverable from the key). NOT covered, and
  caught only by the OLS rank warning (`_OLS_RANK_RTOL`, `fitting/estimators.jl`):
  cross-orbit dependence (a trivial space group — `NoSymmetry` — puts tied images in
  separate orbits; genuine cross-orbit aliasing is unresolvable from the supercell,
  not mergeable, so the reduction deliberately does not fold it), repeated-atom
  members (`AllImages` self-pairs; dropping stays sound, independence uncertified),
  and row-deficient training data. Gate:
  the closed-form CsCl tie-shell testset in `test_salc.jl` (8-fold corner tie: 2
  emitted, the aggregate-zero `Lf = 2` dropped loudly, the survivor ≡ `e₁·e₂`).
  Change the evaluation kernel's member/shift semantics and the reduction's
  aggregation key must follow. **The orbit loop in
  `build_salc_basis` is threaded** (`Threads.@threads`, one task per orbit via
  `_orbit_salcs`): orbits are independent and the output is sorted by `SALCKey`, so the
  basis is byte-for-byte thread-count-independent. The precondition is that the Wigner-D
  cache is **precomputed serially** (`_build_wig_cache`, full bounded `(l, g)` grid) and
  **read-only** inside the loop — `_wig` is a pure lookup, never a `get!`. Any new shared
  mutable state in the per-orbit path (or a lazily-filled cache) reintroduces a race;
  keep it task-local. Gate: `test/unit/test_salc.jl` "build is deterministic / thread-safe"
  plus the oracle, run under `julia -t N>1`.
- **Design-matrix columns are identified by `SALCKey`** (`SALCBasis.keys`, sorted),
  not by construction order. The key layout is `(body, orbit_id, decors, L_S, Lf,
  block)`: `decors` is the sorted `SiteDecor` multiset (this package builds pure-spin
  decors only, so `spin_ls(key)` reads back the v4 `ls` label and `L_S == Lf`), and
  the key must stay **injective**: `block` runs across all
  canonical `l`-orderings that share one sorted `ls` label (a proper-subgroup site
  stabilizer splits a degenerate multiset into several ordering orbits — see
  `test/unit/test_nbody.jl`). `SCEPredictor` re-pairs `jphi` to a basis **by key** on any
  reload (positionally paired only within a session); `fingerprint = hash(sorted keys)`
  guards against cross-basis confusion. This reload is realized by persistence
  (`SCEFitting.load(SCEPredictor, …)` rebuilds `jphi` in basis-key order from the
  per-`SALCKey` coefficients).
- **Persistence schema ↔ the serialized structs** (`io/persist.jl`): `_to_doc` /
  `_from_doc` mirror the fields of `Crystal` / `Lattice` / `SpaceGroup` / `BasisSpec`
  / `SALCKey` / `SALCTerm` / `SALCMember` / `SALC` / `SCEBasis` / `SCEPredictor`. Add or
  rename a field on any of these and both halves (and `test/unit/test_persist.jl`'s
  round-trip) must follow; the space group is rebuilt via `_assemble_spacegroup` from
  the stored fractional ops, and the `UInt64` fingerprint is stored as a string and
  **recomputed** on load (never trusted — `hash` is Julia-version dependent). The TOML
  input reader (`io/input.jl`) mirrors only the *setup* structs (crystal + basis spec
  + symmetry), not the SALCs. Schema v3 stores `BasisSpec` **resolved and dense**
  (per-body `cutoff` matrices, `lsum` with `typemax(Int64)` = uncapped, labels) and
  `_spec_from` still expands legacy v2 scalar-`pair_cutoff` docs — keep that branch
  alive as long as v2 model files circulate. Schema v4 stores SALC members in the
  canonical duplicate-free form; `_basis_from_doc` folds pre-v4 members on load
  (`_canonicalize_members`) — keep that branch alive as long as v3 model files
  circulate, and never write a non-canonical basis.
- **BasisSpec sugar resolution ↔ canonical consumers** (`sce/truncation.jl`,
  `sce/model.jl`, `io/input.jl`): the ergonomic forms (label keys, `"*"` wildcards,
  body-keyed tables, unordered `"A-B"` pair keys, specificity resolution) are expanded
  ONCE, in the `BasisSpec` keyword constructor; everything downstream —
  `SCEBasis`'s fan-out (`_superset_cutoff` → `build_neighbor_list`,
  `cutoff` → `candidate_clusters` per-edge admission, `lsum` → `_enumerate_ls`),
  persistence, `show` — reads only the dense fields. Add a sugar form or change the
  specificity rule → update the TOML reader (`_cutoff_from_input` etc.), the BasisSpec
  docstring, and `test/unit/test_truncation.jl` together.
- **The pointed moment basis rides the decor engine, and three conventions keep it
  honest** (`basis/momentbasis.jl` ↔ `basis/salcbasis.jl` `_orbit_salcs_decors`'s
  `admit` kwarg ↔ `clusters/orbits.jl` `_orbits_from_members` ↔
  `clusters/enumerate.jl` `candidate_clusters`): (1) **member multiplicity is the
  engine's all-orderings convention** — `candidate_clusters` lists every physical
  instance once per site ordering (3! for a 3-body), and the SALC value scales with
  that count, so a pointed enumeration emitting fewer orderings silently rescales
  its columns per orbit (upstream measured half the prototype's 6.0 star oracle) —
  `_pointed_star_candidates` therefore expands every translation class to all 3!
  re-anchored orderings, and any new candidate source must do the same. (2) an
  `admit` predicate handed to `_orbit_salcs_decors` is judged on the lex-min
  representative of each permutation orbit, so its verdict MUST be a
  permutation-orbit invariant — anything built from (decor, species,
  edge-lengths-from-site) data is, because stabilizer perms preserve species and
  are isometries; a rule reading raw site indices is not. (3) the marked-column
  substitution in `_design_moment` is exact only because every pointed label
  carries exactly one mark (a member not marked at the row's atom dies on its
  |u|² = 0 factor before reading the substituted column) — a future label with
  two marks breaks the argument, not just the numbers. `moment_resolvability`
  refuses a member with two ENVIRONMENT spin factors on one reference-cell atom
  (two periodic images of one neighbor) as `UnclassifiableBasis` — the monomial
  signature would overcount the rank (harmonic products on one sphere reduce by
  Clebsch–Gordan; upstream measured 108 symbolic vs 98 actual on the FeGe
  primitive cell). The gates in `test_momentbasis.jl` are INDEPENDENT references
  (a geometric triangle enumeration, the 2√3 shell sum, symbolic ≡ random-design
  rank) — keep them that way; a reference through the SALC machinery gates
  nothing. `MomentSpec.isotropy = true` is upstream's `MomentSpec.soc = false`
  (same screen, opposite name — the ledger row).
- **The moment channel's evaluation axis is keyed by `constraint_mode`, never by
  field presence** (`io/dftsource.jl` `SpinDatum` ctor invariants ↔
  `check_moment_gates` ↔ `io/extxyz.jl` reader/writer ↔ `io/embset.jl`
  `read_embset_pair` ↔ `fitting/momentfit.jl` `MomentDataset` / `predict_moment`):
  mode 4 (direction-pinning
  type) reads `ê` from `directions`, mode 1 (transverse-penalty type) from
  `constraint_axes` — an availability-keyed fallback ("use `mconstr` if present, else
  MW") would silently drop a mode-1 datum with missing axes into the broken `ê_MW`
  coordinate (`‖M‖ → 0` folds the MW direction; upstream measured σ 2.1× on FeRh),
  so the ctor REFUSES mode 1 without axes AND axes without a declared mode. The
  gates (`check_moment_gates`) run at every boundary an axis-carrying datum crosses
  — extxyz generation, extxyz load, the EMBSET pair reader — because archived axes
  are re-verified, never believed. Two subtleties the tests pin (`test_extxyz.jl`):
  a whole-axis flip in mode 1 flips `y` with it (the axis sign is a GAUGE — the sign
  gate must NOT fire, and "fixing" it to fire would refuse every legitimately
  re-gauged archive), and the angle gate is a PERCENTILE (p99 < 5°), because
  collapse rows legitimately carry large single-row angles (FeGe τ0.5 max 5.6° at
  p99 0.14°). The same never-trust-the-flag rule shapes the extxyz reader: spin-only
  vs joint is MEASURED from positions (bitwise across frames), `config_type` is only
  a cross-checked claim — and this package, being pure spin, REFUSES a joint file by
  name (displaced frames, `forces` columns, a `joint` claim; SLCE.jl reads them)
  rather than flattening it to its spins. `moments_bare` (bare `M_int`) and the
  smoothed `magmoms·directions` (`MW_int`) are BOTH stored and neither substitutes
  for the other: the constraint acts on MW (τ and the configuration coordinates),
  the projection target reads M_int, and their ratio is configuration-dependent
  (0.691 ± 0.017 on FeGe τ0.5). VASP vocabulary (OSZICAR/INCAR parsing, λ printing,
  SAXIS, sign conventions) stays in SCETools' generator; the `constraint_mode`
  numbers follow `I_CONSTRAINED_M`, but the KEY is the physical class — another
  code's scheme maps onto class 1 or 4 at its adapter. The extxyz dialect (column
  names, info keys, shortest-round-trip printing) is SLCE.jl's — keep it identical,
  it is the interchange format between the packages.
  The dataset layer (`fitting/momentfit.jl`) carries the rule's consequences:
  `predict_moment`'s runtime default `axes = e` IS the mode-4 identity
  substitution (change one and the train/predict coordinates split); a mode-1
  marked atom with an exactly-zero axis has NO defined target — its rows are
  `defined = false` with `y = NaN` (loud on raw use) and a placeholder design
  row, excluded from BOTH the gated and ungated solves, never imputed; and the
  fit does NO centering and adds NO global intercept because the l = 0 1-body
  `[MARK]` columns already are the per-orbit intercepts μ₀ — wiring the moment
  design into any estimator path that centers columns or appends an intercept
  column double-counts μ₀. `predict_moment` is a spin-reading entry point and
  therefore a validating DOOR (the unit-norm rule with the component bound): `e`
  unit everywhere, `axes` unit on the MARKED columns only (unmarked axes columns
  are never read — a whole-matrix door would refuse the legitimate closed-form
  ê = x̂/ŷ/ẑ readout). `MomentDataset` runs `moment_resolvability` at construction
  and `fit` freezes the vanishing columns to EXACT zero — the same frozen-column
  discipline as the energy side, so `coef != 0` reads structure; weakening either
  half silently reintroduces arbitrary min-norm coefficients for columns no cell
  determines. The moment side is NOT persisted (design record §4.2): `save`
  refuses `MomentBasis` / `MomentFit` / `MomentModel` by name; a later schema
  version adds it — do not "just write the coefficients" into the v5 document.
- **`coeftable` columns ↔ `SALCKey` fields** (`sce/coeftable.jl`): each result row is
  read straight off a `SALCKey` (`body` / `orbit_id` / `decors`→comma string / `L_S` /
  `Lf` / `block`) plus `jphi`; the `J` column pairs with `basis.salc_basis.keys` **positionally**
  (same order as the design matrix). Add or rename a `SALCKey` field → update the row
  builder, the `Tables.Schema`, and `test/unit/test_coeftable.jl`.
- **DFT training-torque target ↔ the model torque convention** (`io/dftsource.jl`): the
  training torque carried by a `SpinDatum` is `τ_a = m_a × B_a` (`B` = constraining field),
  which must stay the *same* physical quantity, sign, and `3×n_atoms` config/atom/`xyz` layout
  as the model's `predict_torque = −e_a × ∂E/∂e_a` (the design-matrix convention) — **both**
  are the physical / Landau–Lifshitz torque `m × B_eff`. Flip one side only and the co-fit
  silently biases; flipping **both** (as done when the package moved from the `+e×∇E`
  energy-rotation-gradient to this `−e×∇E` Landau–Lifshitz convention) leaves `J` unchanged.
  **Both sides are gated, independently, and the gates do not cancel.** Model side:
  `test_torque.jl` "predict_torque = −e × ∇E (finite differences, anisotropic)" builds its
  reference from `predict_energy` by central differences — an energy surface carries no
  torque convention, so flipping `predict_torque` (or the design kernel, caught by the
  `torque_weight = 1` recovery test in the same file) turns it red. Training side:
  `test_dftsource.jl` "torque convention: τ = m × B, closed form" pins a hand-written
  literal with the wrong sign named explicitly in the comment. A SIMULTANEOUS flip of both
  sides is a **gauge**, not a bug — `J` and `predict_energy` come out bit-identical — and it
  cannot pass anyway, since it would have to edit both a closed-form literal and an
  energy-derived reference. What is genuinely NOT gated here, and is the one thing to think
  about by hand, is the semantics of the external file: whether VASP's `lambda*MW_perp`
  block is `+B` or `−B`. That single bit is anchored only by the (joint-family) SLCETools
  oracle (parsers vs Magesty), i.e. against a prior implementation rather than against
  physics. [Backported from SLCE.jl 1495e44.]
  **DFT-code I/O is confined to `AbstractDFTSource` adapters in the SCETools.jl package**
  (`SCETools.VASP`), which produce `SpinDatum`s (rotating moments / field from the `SAXIS`
  frame by `Rz(α)·Ry(β)`); the core consumes only `SpinDatum`/`SCEDataset` and stays
  DFT-code-agnostic. The VASP parsers are cross-checked against Magesty in SCETools's oracle.
  The `SpinDatum` torque sign defined here is the convention source the adapters must match.
- **Sunny export conversion ↔ the energy reconstruction** (`sce/bilinear.jl` — matrices /
  extraction / gate — plus `interop/sunny.jl` — primitive unfold — and
  `ext/SCEFittingSunnyExt.jl`): `_l1_pair_matrix` / `_l2_onsite_matrix` must satisfy
  `eₐ'·M·e_b = Σ folded·Z·Z` (the gate is the `Z₁`/`Z₂` contraction test; the tesseral
  constants `N1`/`A2`/`B2` are defined once, in `basis/Harmonics.jl`); the per-bond
  matrix is `jϕ·(4π)^(N/2)·M` and the two directed members `(a,b,R)`/`(b,a,−R)` fold into
  one matrix on the canonical `a≤b` bond (reverse transposed). The whole chain is checked
  **without Sunny** by `_reconstruct_energy ≈ predict_energy − j0`; the extension then
  rescales by the effective spin (`scaling = :moment` sets `J = M/(SₐS_b)`; `:coupling`
  folds `S_eff` into the couplings at a placeholder `Moment`) so the `Sunny.System` energy
  matches. Only `ls=[1,1]`/`ls=[2]` are
  representable — every other channel must be **reported as skipped**, never silently
  dropped. Change a harmonic normalization or the `(4π)^(N/2)` scale → both the matrix
  formulas and the energy gate move together.
- **Fitted-model introspection ↔ the per-term scale convention** (`sce/introspect.jl`,
  `test/unit/test_introspect.jl`): `multipole_terms` is the **public, stable** view downstream
  packages (the `SCETools.jl` mean-field samplers) read instead of `model.basis.salc_basis.salcs` /
  `SALCMember` / `SALCTerm`. It returns the **raw** fitted `jϕ` as `coef` and leaves the per-N
  scale `(4π)^(body/2)` to the consumer — the scale lives in exactly one place (the
  reconstruction gate `_energy_from_terms`), so do **not** also apply it inside
  `multipole_terms`. `bilinear_terms` is a thin public wrapper of the general
  `_bilinear_terms` extraction (in `sce/bilinear.jl`), so its numerics move with the
  Sunny coupled-site above.
  **`MultipoleTerm.ls` keeps its name and its meaning** (the per-site spin ranks, one
  per atom, `length(ls) == body`, `size(folded) == Tuple(2l+1 …)`): SCEMonteCarlo.jl's
  `TiledHamiltonian` ingest enforces exactly that contract with `throw`s, derives its
  `(4π)^(body/2)` scale from `body`, and hashes the VALUES of `t.ls` into its
  `model_fingerprint` (so a rename is invisible to it, a changed term list invalidates
  every JLD2 checkpoint). The decor vocabulary stops at `multipole_terms`' refusal of
  a decorated basis; it does not leak into this struct.
  Add or rename a `MultipoleTerm` field → update the gate and any downstream consumer
  (in the revived spin family that is SCEMonteCarlo.jl's `TiledHamiltonian` ingest; the
  joint family's SLCETools reads the SLCE equivalent through `mfa/bridge.jl` — the old
  `sce_bridge.jl` name is dead, grep the consumer package rather than trusting a
  filename here). [Pointer repaired per SLCE.jl 1495e44.]
- `solve_coefficients(est, X, y; groups)` receives a **column-centered** `X` (⇒ the
  solver adds no intercept; `j0` is recovered analytically in `fit`). Every estimator —
  in-tree or in an extension — must honor this. `groups` (optional) labels rows from the
  same physical sample (in a co-fit, a configuration's energy row and its
  torque-component rows share a label); a resampling estimator (CV-based `ElasticNet` /
  `Lasso` / `AdaptiveLasso` in `ext/SCEFittingGLMNetExt.jl`) must keep same-label rows
  in the same fold so CV does not leak within-configuration structure. The analytic / adapter
  estimators (`OLS` / `Ridge` / `AdaptiveRidge` / `PrecomputedPilot`) ignore it. The GLMNet
  solve uses `intercept = false` + column `standardize` and selects λ by configuration-
  grouped, seeded CV (`:lambda_min`/`:lambda_1se`); change the centering/whitening in
  `fit` and the penalty scale (`λ·std`) moves with it. `AdaptiveLasso` runs its `pilot`
  through `solve_coefficients` (forwarding `groups`), then a weighted-L1 GLMNet solve with
  fixed `penalty_factor`; `AdaptiveRidge` is a pure-core reweighted-ridge loop sharing the
  centered-`X` contract. Validated in the separate `test/glmnet/` env (GLMNet-backed) and
  `test/unit/test_fit.jl` (core `AdaptiveRidge` / `PrecomputedPilot` solves), never mixing
  the two (GLMNet absent in the core suite).
- **GCV ↔ `_assemble_problem` ↔ `islinear` ↔ the GAR weight map** (`fitting/selection.jl`,
  `fitting/estimators.jl`): `gcv`/`effective_dof` reassemble the design through
  `_assemble_problem` (change the centering/whitening and the score moves with `fit`),
  are gated by `islinear`, and recompute the converged penalty diagonal through the
  **same** functions the solvers iterate — `_gar_weights!` (the single definition of
  `wⱼ = v_g/(‖β_g‖² + p_g·ε)`; `_penalty_diagonal` has one method per linear estimator,
  and `AdaptiveRidge`'s `1/(β² + ε)` must stay in sync with its solve loop). Change a
  weight formula in the solver and the `_penalty_diagonal` method, the design-notes §13
  derivation, and the dense-hat-matrix tests in `test/unit/test_selection.jl` move
  together. `select_fit`'s alive-group rule is the `refit` scaled-magnitude support rule
  (`|jϕⱼ|·‖X[:,j]‖ > threshold`) applied per group — change one side and the other (and
  the E2E cost-recomputation test) follows; `select_support` reuses the same rule for its
  per-point alive/cost columns while delegating each point's fit to `refit` itself
  (column-wise), so all three move together. `select_support`'s evalset score and
  `cross_validate`'s holdout score share the prediction-space convention
  (`y_E − (j0 + X_E·jϕ)`, `y_T − X_T·jϕ`, combined as `(1−w)·MSE_E + w·MSE_T`) —
  change `fit`'s objective normalization and both scores must follow.
  `salc_groups`/`group_costs` assume sorted
  `SALCBasis.keys` and canonical (v4) members; the entry key `(atoms, shifts, ls, index)`
  identifies one distinct contraction entry, and the cost prices a SWEEP — one
  site-program slot per member site of each entry (`Σ length(ls)`), not one per
  entry (that is the energy program, walked once per run; SLCE.jl a596ea3) —
  change either representation and re-check the brute-force union test.
- **`fit` ↔ `refit` share `_assemble_problem`** (`fitting/fit.jl`): the `(X, y, xbar, ybar,
  groups)` centering/whitening assembly lives in one helper so the two build identical
  designs — change the centering or whitening there and **both** move together (the oracle
  pins `fit`'s numerics). `refit` re-solves on the scaled-magnitude support
  `|jϕ_j|·‖X[:,j]‖ > threshold` of an existing fit (a column sub-matrix), so it rejects a
  `PrecomputedPilot`-backed estimator (fixed-length pilot vector ≠ support length).
  **`SCEFit.residuals` is the energy-only residual** `y_E − (j0 + X_E·jϕ)` (not Magesty's
  combined whitened residual); the diagnostics report energy and torque blocks separately
  (`residuals_energy` returns the stored vector, `residuals_torque`/`rss_torque` recompute
  `y_T − X_T·jϕ` and validate `has_torque`; `r2_*`/`rmse_*` build on `rss_*`).

## Upstream divergence ledger (SLCE.jl ↔ this package)

This package is the pure-spin carve-out of SLCE.jl and stays a **standalone**
package (no dependency on it). Ported files are kept close to upstream so patches
apply both ways, but the divergences below are deliberate and must survive a
"re-sync". Before porting anything across, read this table; the polarity column is
the one that bites.

| Divergence | SLCE.jl spelling | This package | Polarity / caution |
|---|---|---|---|
| Decor engine screen | `_orbit_salcs_decors(…, labels, soc::Bool, wcache; lmax_by_species, pmax_by_species, admit)` — 7th **positional**; `soc = true` keeps every `L_S` | `_orbit_salcs_decors(…, labels, wcache; isotropy::Bool)` — **required keyword**; `isotropy = true` keeps `L_S = 0` only | **Opposite meaning in the same slot.** A verbatim upstream call must be a `MethodError` here; never make `isotropy` positional or give it a default |
| Moment-basis screen | `MomentSpec(; soc = false)` — `soc = true` keeps every `L_S` | `MomentSpec(; isotropy = true)` — `isotropy = false` keeps every `L_S` | Same polarity trap as the engine row, one level up; the field is named `isotropy` here and forwarded as the engine's keyword |
| Path screen placement | `_decor_coupled_bases(slots)` builds every path; the screen is applied afterwards | `_decor_coupled_bases(slots, isotropy)` hands `AngularMomentum.build_real_bases` a `keep` predicate so a rejected path never builds its tensor | Same SALCs, different call shape; port logic, not signatures |
| Admission | `_admit_assignment(t, species, …)` — a production, species-resolved rule | only the `admit` hook; callers (tests) transcribe the per-species `lmax` | The pointed builder (D4) will need its own mark-aware rule; upstream's is the reference, not a drop-in |
| Function-space reduction | none | `_function_vector` / `_reduce_orbit_salcs` (pure-spin only; refuses decorated SALCs, message = wiring checklist) | Exists only here; upstream ports nothing back |
| `SolidHarmonics` | values + Euclidean gradient API (`solid_harmonics_grad[!]`, `grad_Rlm`) + `solid_harmonic_poly` (the ASR and lattice-side builders) | **values only** (347 → 240 lines); the value recurrence is upstream's line for line | No force rows here; do not re-port the gradient "because upstream has it" — count what the production path actually reads (`R₀₀ ≡ 1`) |
| Test oracles | `CountingOracle` (852 lines), `_ls_block_stats` (C-2 block diagonality), plus the Cartesian projector since `08743d1` | the ~45-line Cartesian projector only; C-2 deferred to wiring | Counts agree; the projector shares no code with the SALC machinery in either package |

Closed on 2026-08-21 — now identical in both packages, no longer divergences:
content-based `hash(::SiteFactor)` / `hash(::SiteDecor)` (SLCE `4122fba`), the
comma-free `u(k:l)` decors token (`e828524`), and the validating persist reader
(`7546616`, plus the key ↔ members check that followed the second review).

## Tests

| Command | Purpose |
|---|---|
| `julia --project -t 4 -e 'using Pkg; Pkg.test()'` | unit + Aqua (default) |
| `TEST_MODE=all julia --project -t 4 -e 'using Pkg; Pkg.test()'` | unit + Aqua + JET |
| `TEST_MODE=jet julia --project -t 4 -e 'using Pkg; Pkg.test()'` | JET type-stability |
| `julia --project=test/oracle test/oracle/runtests.jl` | from-scratch numerics vs pinned Magesty |
| `julia --project=test/sunny test/sunny/runtests.jl` | real `Sunny.System` energy vs SCE (extension) |
| `julia --project=test/glmnet test/glmnet/runtests.jl` | GLMNet Lasso / elastic-net solve (extension) |

The core suite (`runtests.jl`) dispatches on the `TEST_MODE` env var
(`default`/`all`/`unit`/`aqua`/`jet`) and never depends on Magesty. The oracle,
Sunny, and GLMNet suites are separate environments that carry the heavy/optional
dependency (a pinned Magesty.jl / Sunny / GLMNet) the core deliberately omits.

## Git (this rebuild)

Local git (`add` / `commit` / branch) is pre-authorized for this exploratory
rebuild — no per-action confirmation. Only remote operations (`push`,
registration) require confirmation. Commit directly to `main`; the `feat/v0-core`
working branch was retired once v0 landed (it tracked `main` one-to-one, so it
gave no isolation). Cut a short-lived topic branch only for genuinely risky or
experimental work you do not want on `main` yet.

## References

- `STYLE_GUIDE.md` — package-specific style deltas.
- `SPEC.md` — realized architecture, types, public API.
- `docs/design-notes.md` — why the rebuild diverges from Magesty (the refinements).
- `CHANGELOG.md` — what landed in the v0 slice.
- `examples/heisenberg_chain.jl` — runnable end-to-end (recovers `J`);
  `examples/kagome_threebody.jl` — 3-body / multi-term SALCs, energy+torque co-fit.
- `references/` — supporting literature (notes tracked, PDFs local-only).
