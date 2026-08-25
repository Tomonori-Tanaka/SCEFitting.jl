# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/); this package predates a tagged
release, so everything lives under *Unreleased*.

## [Unreleased]

### Fixed — mechanical audit items: interchange round-off, TOML kinds, IRLS reporting (2026-08-25)

- **Every geometry comparison a reader makes now holds to a band, not exactly.**
  `read_extxyz(...; reference)` measured spin-only vs joint by `all(iszero, pos - refc)`,
  where `refc` is a freshly recomputed `vectors * frac`. This is an interchange format
  (the same dialect as SLCE.jl), so the file's writer is routinely a different build and
  the last bits differ; a 1-ulp (~2e-15 Å) mismatch was reported as *"positions differ
  from the reference crystal's — displaced (joint spin–lattice) data are not
  representable in this pure-spin package"*, a claim about the physics that was false.
  Frame-to-frame, lattice, and file-against-reference now all use one absolute band,
  **`_REF_GEOM_ATOL = 1e-6 Å`**, and each refusal reports the deviation it measured.
  The value is set by the format, not by round-off: ASE's extxyz writer defaults to
  `%16.8f`, a text grid whose own rounding reaches 5e-9 Å, so a band at 1e-8 Å would
  leave a factor of two and a six-decimal producer would reproduce the very verdict the
  band removes. 1e-6 Å clears an eight-decimal grid by ~200× and still sits three orders
  below the ≳ 1e-3 Å displacement the distinction is about. **It is therefore also a
  floor**: a structure displaced by less than 1e-6 Å is not representable through a file
  and reads as sitting at the reference.
- **`SpinDatum` screens every field for finiteness**, not only `moments_bare` and
  `constraint_axes`. The file readers screen numbers as they parse, but an adapter that
  builds a `SpinDatum` directly — the production path — reached the constructor with
  whatever the SCF produced, and one diverged frame turned every fitted coefficient into
  `NaN` with nothing naming the configuration. `directions` is deliberately still left
  to the dataset constructor, which checks norm, pole margin and finiteness together.
- **Every numeric and boolean door of the TOML input is kind-checked.** `Bool <: Real`
  and `Bool <: Integer`, so `cutoff = true` was read as 1.0 Å and `lsum = true` /
  `nbody = true` as 1 — a silently different model. `[moment]` already refused this in
  three places; the asymmetry was not a decision. The whole file now shares
  `_is_toml_int` / `_is_toml_number`:
    - the integer doors require a TOML **integer**, so `nbody = 2.0` and `lmax = [2.0]`
      are refused too, as is a quoted `lmax = "2"` — which used to iterate the string
      and load as the codepoint `[50]`;
    - the boolean doors require a TOML **boolean**, from the other side:
      `[interaction].isotropy = 1` and `[structure].pbc = [1, 1, 1]` are refused;
    - `[structure]`'s own geometry (`lattice`, `positions`, `species`) and
      `[interaction].tie_tol` go through the same tests, naming the offending entry;
    - `[symmetry].tol` is spglib's symprec, a **cartesian distance in Å**, and nothing
      downstream bounded it — `tol = true` meant a 1 Å symprec, which merges
      inequivalent sites and reports a larger group, so the basis is silently a
      different one. It is now required to lie in `(0, 0.1]`.
- **A reweighted ridge that exits on `max_iter` says so.** Both IRLS loops discarded
  the last relative change and reported nothing on hitting the iteration cap. That
  matters beyond the coefficients: `islinear` is `true` for `AdaptiveRidge` /
  `GroupAdaptiveRidge`, so `gcv` and `effective_dof` rebuild the penalty diagonal FROM
  the returned coefficients and score the smoother it implies — a smoother that was
  never solved. A direct solve warns on the spot, deliberately without `maxlog` so a
  later non-convergence is never hidden by an earlier one; `select_fit`, which runs one
  solve per λ and another per (fold, λ), instead **counts** them and reports once —
  how much of the path is affected is the part a user can act on.
- The pointed moment basis's self-image comment claimed the `allunique` guard was a
  second lock, "since a minimum-image neighbor list has no self-pairs". That is a
  statement about the *centre's* neighbours, and the `N!` re-anchoring walks the mark
  around the star, so an environment of the original centre becomes the mark and
  another environment can sit on an image of it — the guard is the only lock. The
  spoke test in `admit` is not a second one either: `_dmin2_matrix` leaves its diagonal
  at `Inf`, so that comparison is identically true for exactly this case. Comments only;
  no behaviour change.
- `with_lambda` gained the tests it had none of — a `fieldnames` sweep asserting that
  λ moves and every other field, the penalty metric and its provenance included, is
  carried through. The package documents this as unobservable through a fit (a dropped
  metric is indistinguishable from a deliberate uniform one), so the contract is
  asserted structurally.

### Fixed — a neighbour list keeps its species-pair radii, not their maximum (2026-08-25)

- **`NeighborList.cutoff` is now the symmetric per-species-pair radius matrix the list
  was built with** (a scalar build broadcasts it), replacing the `Float64` maximum.
  `maximum(nl.cutoff)` recovers the old value.
- **`candidate_clusters` gates every `AllImages` clique edge by its own species-pair
  radius.** An `N ≥ 3` clique is grown from one anchor's neighbours and its remaining
  edges are re-checked here; with only the maximum radius on hand, the anchor's edges
  were gated by the species radii and the rest by that maximum. Two rules in one
  predicate, so **which sites a cluster happened to be reachable from decided its
  fate**. Measured on Fe + 2 Te with `d(Fe,Te) = 2.5 Å`, `d(Te,Te) = 4.0 Å` and radii
  `[6 6; 6 3]` under `NoSymmetry`: the `{Fe,Te,Te}` triangle came out with 2 of its 6
  anchor-variants where its own Te–Te radius asks for 0. The predicate is written in
  species and distances, hence group-invariant, so `build_clusters`' closure assertion
  never fired; and because `_canonicalize_members` *sums* the anchor-variants of an
  instance rather than deduplicating them, an orbit short of its variants carries a
  design column scaled against its siblings — the ratio of one `J` to another.
- **Per-order radii may only trim the list.** A `cutoff` entry above the list's own
  radius for that species pair now raises rather than silently under-generating: an
  edge the list never enumerated cannot be recovered here, so the wider radius would
  reach a cluster's re-checked edges and not its anchor's. Same failure, other
  direction — measured as 2 of 6 again, with the list at 2.5001 Å and `cutoff` at
  6.0 Å.
- Production fits were never affected: `SCEBasis` always passes
  `cutoff = spec.cutoff` with the list built at the element-wise maximum over the
  orders, and the pointed moment basis uses the `nbody ≤ 2` `MinimumImage` path. The
  `MinimumImage` path was already anchor-independent — `_dmin2_matrix` re-derives the
  species admission from the list's own contents. Reachable only through the public
  `build_neighbor_list` (matrix form) → `build_clusters` / `candidate_clusters`.

### Fixed — an aperiodic axis now restricts the space group (2026-08-25)

- **`analyze_symmetry` intersects the backend's group with the crystal's declared
  periodicity.** A backend analyses the cell as a fully periodic 3D crystal — Spglib is
  not told about `Lattice(...; pbc)` and has no way to be — while the neighbour list
  refuses to emit an image along an aperiodic axis. `_build_map_sym` then folded mod 1
  on all three axes unconditionally, so an operation that closes only through the
  artificial periodicity was accepted and surfaced downstream as
  `build_clusters`' closure assertion, blamed on "the image selection and the tie
  tolerance" — a cause it does not have. Reproduced on a three-layer slab straddling
  `z = 0` with `pbc = (true, true, false)`; the same slab placed at the centre of the
  cell, and either placement fully periodic, built fine.
- Operations whose rotation mixes a periodic with an aperiodic axis are refused on the
  rotation alone: they send a lattice translation into a direction with none to receive
  it, so the declared translation lattice is not mapped onto itself.
- **Translations are re-seated, not just filtered.** A backend reports `t` modulo a
  lattice translation, which along an aperiodic axis is not an identification one may
  make — exactly one representative is the operation the finite structure has. A slab
  centred at `z = 1/2` has a mirror there, and Spglib is entitled to report it as the
  mirror at `z = 0` with `t_z = 0`; the assembler now searches for the single integer
  shift along the aperiodic axes that makes every atom match exactly, and keeps the
  operation with that representative. So the centred slab keeps its full group, and
  `_site_image` produces a zero cell shift along the aperiodic axes by construction.
- The kept set is a **subgroup** (the aperiodic match is exact, so it composes; the
  identity satisfies it; a bijection's inverse does too) and is re-validated with
  `_validate_ops`. Using a subgroup never over-reduces: the basis is larger than the
  fully periodic one, never short. A fully periodic crystal returns before any of this
  and is bit-for-bit unchanged.
- **An aperiodic axis is never wrapped, so positions may legitimately span cells along
  it** (`Crystal` wraps only the periodic axes). Two same-species atoms an exact cell
  apart there are indistinguishable to the mod-1 comparison the matcher uses to find a
  candidate, which used to make even the IDENTITY fail as "atoms 1 and 2 share the
  image 1" — under `NoSymmetry()` too, so there was no fallback. The whole-operation
  shift is therefore searched over the candidates atom 1 can map to, nearest zero
  first, instead of being read off the first hit; "no consistent shift" is a drop, not
  an error. That seed list is complete (atom 1 must map to something), so nothing
  legitimate is lost, and a drop is the conservative direction anyway.
- **The load-bearing invariant is checked directly.** `_check_zero_aperiodic_shift`
  asserts, for every kept operation and atom, that `round(W·x_a + t − x_b)` vanishes on
  each aperiodic axis — the property `_site_image` needs. `_validate_ops` cannot stand
  in for it: its `_tclose` folds mod 1 on all three axes, so a re-seated operation and
  the representative it replaced are the same object to it. It is now also run on the
  kept set whether or not anything was dropped, since re-seating is a change on its own.
- **`PERSIST_SCHEMA_VERSION` 5 → 6, and a pre-v6 document whose crystal declares an
  aperiodic axis is refused on load.** Its stored operations are the full group a
  backend reported for the fully periodic cell and its SALCs were projected with those,
  while this build re-derives the restricted subgroup — pairing them would give a basis
  whose group is not the one its columns came from, and nothing downstream would notice
  (the `SALCBasis` fingerprint hashes keys, not operations). A document's `pbc` cannot
  be changed from the loader, so the answer is to rebuild the basis. Fully periodic
  documents of every readable version load unchanged.
- The dropped count is warned once per `analyze_symmetry` call (no `maxlog`: a user
  building several slabs needs to be told about each), and `symbol` gains a
  `" (pbc subgroup)"` suffix so a group that is not the reported one cannot be mistaken
  for it. `guide/basis.md`, the `Lattice` docstring and the `analyze_symmetry` docstring
  say all of this, including the alternative: for a vacuum-padded slab meant as a 3D
  crystal, keep `pbc = (true, true, true)` — the vacuum and the cutoff already keep the
  images apart.

### Fixed — `select_fit(:cv)` reported the objective divided by the row count (2026-08-25)

- **`select_fit(...; criterion = :cv)` no longer divides the pooled out-of-fold sum by
  the informative row count.** `_assemble_problem` already row-scales the design by
  `√((1−w)/n_E)` and `√(w/n_T)`, so a fold's squared holdout residual is already
  per-row; summed over folds (every row held out exactly once) the total IS
  `(1−w)·MSE_E + w·MSE_T` — the objective the docstring promises, and the scale
  `cross_validate`'s `pooled_score` and `select_support`'s `score` report. The extra
  `./ neff` put the reported number a further factor of `neff` below all three:
  measured on 60 energy configurations, `select_fit(:cv).score = 1.576e-10` against a
  `cross_validate` pooled score of `9.759e-9`, a ratio of 61.9 ≈ `n_E`.
- **No selection changes.** `_select_pareto` is invariant under a positive uniform
  factor, and the factor is constant along the λ path, so every path this ever selected
  it still selects. Only the number a caller reads (and any `delta` calibrated against
  a hand-computed error) was wrong.
- The docstring now also says plainly that `:cv` and `:gcv` are **not** on a common
  scale: `:gcv` reports `n·RSS/(n−df)²` over the same already-row-scaled rows, so it
  runs roughly `n_eff` below the objective. Each criterion ranks its own path; the two
  numbers must not be compared to each other.
- Gated two ways in `test_selection.jl`: an exact identity that the assembled
  squared-residual sum IS `(1−w)·MSE_E + w·MSE_T` (hand-written from the definition in
  physical units, touching no part of the selection driver), and a scale check against
  `cross_validate`, which recomputes the same pooled quantity by an unrelated route.

### Added — the DFT cross-observable gate (2026-08-25)

- **`examples/cross_observable_bcc_fe.jl`**, run by CI's *runnable examples* job.
  Every existing torque gate closes a loop the package draws itself: the analytic
  Heisenberg torque pins `predict_torque` and the torque design against `−e × ∇E`, and
  the unit suite pins `predict_torque` against finite differences of `predict_energy`.
  None of them can see the one step that leaves the package — VASP writes a constraining
  field `B_a`, `SpinDatum` turns it into `τ_a = m_a × B_a`, and no internal check can
  tell whether that is the same torque as `−e_a × ∂E/∂e_a` of the energies in the same
  file. A sign or a factor there fits happily and is wrong by that factor.
- The example fits the real bcc Fe 4×4×4 fixture (128 atoms, 50 constrained
  configurations) on each observable alone and predicts the other: R² 0.98051 on the
  torques from an energy-only fit, 0.98642 on the energies from a torque-only fit, and
  a least-squares scale `c = 0.87291` between the two coefficient vectors. Rescaling
  the targets by −1, 2 and 0.5 each trips at least two of the three gates (the `c` band
  trips on all three), so the convention is confirmed with orders of headroom.

### Changed (breaking) — `MomentSpec.lsum` is per body order (2026-08-25)

- **`MomentSpec.lsum` is now a `Vector{Int}`, one entry per body order**, indexed by
  the order itself (`lsum[N]`, from `N = 1`) and read through `_label_lsum`. The
  keyword accepts what `BasisSpec`'s `lsum` accepts, through the same resolver
  (`_resolve_lsum`): a scalar caps every order, body-keyed pairs or a `Dict` cap one
  order each with unnamed orders staying uncapped, `nothing` caps nothing. A positional
  vector is refused — with no offset the body-keyed form is unambiguous.
  BREAKING CHANGE: code that reads the field (`spec.lsum == 4`) must read an entry.
  Code that *writes* the keyword is unaffected, and a scalar spelling produces a
  bit-identical basis: the screen is still `Σl` even and `≤ cap`, with the same cap.
- Why per order: a body order starts at `Σl = 2⌈(N−1)/2⌉` — floors 0, 2, 2, 4 — so one
  global cap starves the high orders. At `lsum = 4` the 2- and 3-body sectors get two
  even levels while the 4-body sector gets one, and the number that opens the 4-body
  sector at its floor is the same number that removes the 3-body sector's `Σl = 6`
  content (17 columns on the recorded 54-atom bcc Fe basis). The four-body measurement
  that motivated this could not express "3-body uncapped, 4-body capped" at all and
  needed a separate control run to undo the confound.
- `[moment].lsum` accepts the body-keyed table (`[moment.lsum]` with bare-integer
  keys), which it previously refused by name. Unlike `[moment.cutoff_star]` the table
  is **partial** — an unnamed order is uncapped, whereas a missing radius would leave a
  star order silently at `cutoff_pair`.
- The empty-sector warning now names the cap that is responsible (`lsum[$body]` and its
  value) the way it already named `cutoff_star[$(body-2)]`.

### Added — the pointed moment basis takes four bodies (2026-08-25)

- **`Threads.@threads :greedy`** on the `MomentBasis` orbit loop, the
  `penalty_metric(::SCEBasis)` column loop and `_design_moment`'s column loop. All three
  emit work in ascending body order, so the expensive items are contiguous at the END —
  and the default schedule (and `:dynamic`, which differs only in thread affinity) cuts
  the range into one contiguous chunk per thread, putting every high-body item in the
  last chunk. Bitwise identical at any schedule; the 4-body build drops 5.73 s → 4.36 s
  on the recorded fixture.

- `MomentSpec(; nbody)` accepts 4. The enumeration, the labels and the star
  candidates were rewritten for general `N` first (no numerical change: the pointed
  pin's members, tensor entries and fit are bit-identical across that rewrite), and
  the cap now sits where the test oracles stop rather than where the code did.
- One of the three hard-coded spots was silent: the constructor's star block was a
  single `if spec.nbody >= 3` that built body order 3 and nothing else, so a raised
  cap alone would have returned the 3-body basis without a word.
- A body order starts at `Σl = 2⌈(N−1)/2⌉` — every environment slot needs `l ≥ 1` and
  time reversal keeps only even `Σl` — so the 4-body sector begins at `Σl = 4` and its
  naive lowest member (a rank-0 mark with three `l = 1` environments) is absent: the
  only rotational invariant of three vectors is the triple product, which is
  time-reversal odd. (It is a pseudoscalar, but that is not the reason — a chiral
  space group allows pseudoscalars; the screen is time reversal.)
- Oracles for the new sector: the star candidates match an independent brute force
  written from the definition (four cells including Wigner–Seitz ties, three cutoffs,
  `N = 3` and `4`), with the `N!` ordering multiplicity read off the production set;
  the absolute column normalization is derived rather than captured —
  `4!·(4π)^{3/2}·(1/√5)·(3/4π)·√(15/8π) = 24·3√(3/2)` against the closed-form
  invariant `(eⱼ·e_l)(e_k·e_l) − ⅓(eⱼ·e_k)` on a P1 fixture chosen so the Reynolds
  projector acts on one dimension; covariance under an arbitrary `SO(3)` rotation of
  spins and axes together; bitwise time reversal; and the resolvability gate at `N = 4`
  on a 27-atom supercell against the numerical rank of a random design.
- Review-panel follow-ups (2026-08-25): a requested body order that no label can
  reach is now **loud** — the build warns, naming the sector's `Σl` floor and whether
  the labels or the geometry are at fault, and `show` reports the body order actually
  built rather than the one requested. The `UnclassifiableBasis` refusal names the two
  cheap remedies (narrow `cutoff_star`, step `nbody` back) beside the expensive one,
  and its guard now also catches an environment landing on an image of the mark
  itself. The `MomentBasis` orbit loop is threaded, as `build_salc_basis` already was;
  the resolvability census is one grouping pass instead of orbits × columns × members;
  `_design_moment` schedules dynamically, its per-column costs being orders of
  magnitude apart. `_combinations` / `_permutations` live beside `_ordered_subsets` in
  the cluster layer.
- The star cutoff rule is unchanged at higher `N`: only the `N−1` mark–environment
  spokes are cut, the environment–environment edges stay free. A star has a
  distinguished centre, so each environment site is fixed by its own spoke and two
  orbits cannot carry the same monomial — the energy side needs its compact-cluster
  rule precisely because its clusters have no centre.
- **`cutoff_star` is now per star order** (BREAKING for anything reading the field:
  it is a `Vector{Matrix{Float64}}`, body order `N` at index `N - 2`, read through
  `_star_cutoff`; in TOML, a body-keyed table under `[moment.cutoff_star]` whose keys
  must cover exactly `3:nbody`). A scalar or a matrix still broadcasts to every order,
  so every existing spelling builds the same basis. This is what makes a four-body
  probe reachable at all: at the 3NN star radius the 4-body sector of bcc Fe 3×3×3 is
  2,774,736 members and 115,614 P1 orbits against 5,400 design rows, while the same
  sector cut to the first shell is 72,576 members and 3 orbits. Per-*spoke* radii
  remain impossible — the label is a decor multiset, so permuting the environment
  sites leaves it unchanged and "first spoke short, second long" has no
  symmetry-invariant meaning; a total-spoke-length or diameter cap would, and is not
  implemented. Below body order 3 an explicit `cutoff_star` is refused rather than
  stored and never read.
- `bench/bench_moment.jl` (`make bench-moment`) times the five stages separately —
  member generation, orbit reduction, the whole build, the resolvability gate, and the
  two per-fit sweeps — plus TTFX from a cold child process. The wall turns out to be
  the SALC projection: enumeration is under 2 % of the build at either order, and
  going from three bodies to four costs +2.2 s of projection for 12 extra columns
  (`bench/BENCH_LOG.md`).

### Added — λ selection for the moment channel (2026-08-24)

- `cross_validate(::MomentDataset, estimator; nfolds, seed)` returns a
  `MomentCVResult`: configuration-grouped K-fold cross-validation of the pointed
  channel, the honest criterion for choosing λ now that the moment fit has a
  penalty worth tuning. Folds are assigned by configuration, so the rows of one
  configuration — one per marked atom, sharing its directions — never split across
  the train/holdout boundary, and each fold re-solves with the same frozen column
  set as the full dataset so the folds compare like for like.
- A fold whose training rows miss a marked orbit entirely is refused by name: that
  orbit's μ₀ intercept would be unidentified on the fold rather than merely
  under-determined. A `PrecomputedPilot` is refused for the usual leak reason.
- Training and scoring run on the gate-kept rows — the decomposability gate says
  where the pointed model is defined, not how well it generalizes — with
  `score_defined` reporting the rejected rows separately, as disclosure.
- `effective_dof(::MomentFit)` and `gcv(::MomentFit)` complete the fast pair. Both
  reconstruct the design the fit actually **solved** (gate-kept rows, vanishing
  columns frozen out, estimator reduced to match) rather than the full `ds.X`;
  reconstructing the full design would pass every length check and still charge
  degrees of freedom to columns the solve froze. There is no `+1` intercept term:
  the pointed design carries its μ₀ columns explicitly, and being unpenalized they
  each cost a full degree of freedom through `tr(H)`.

### Changed — **BREAKING**: a basis-intrinsic penalty metric (2026-08-24)

- `Ridge`, `AdaptiveRidge` and `GroupAdaptiveRidge` take a `metric`: a per-column
  penalty scale `m`, with `metric_provenance` recording what it was built from.
  `λ·Σⱼβⱼ²` is not invariant under rescaling a design column, and SALC column norms
  are set by basis conventions — an orbit's member count, and the ordering
  multiplicity `_canonicalize_members` folds into the tensor — rather than by
  physics. A larger column norm means a smaller coefficient at the same physical
  effect, hence *less* shrinkage, so the plain penalty carried an accidental prior
  in favour of large orbits and high body order, on top of the deliberate one
  `cost_weights` states through `theta`.
- `penalty_metric(basis; torque_weight, nconfig, seed)` builds it as the reference
  norm of the column **as the estimator sees it** — the assembled, centered /
  whitened design at that `torque_weight`:
  `mⱼ(w) = (1 − w)·Var[Φⱼ] + w·(1/3n_atoms)·E[Σ_a ‖(∂Φⱼ/∂e_a) × e_a‖²]`, over
  uniform-random spin configurations. The `1/(3·n_atoms)` is the per-row average
  the assembly's `√(w/n_T)` already applies. `penalty_metric(mb)` is the pointed
  channel's `E[Φⱼ²]`.
- The metric sits in the **denominator** of the adaptive weight maps —
  `wⱼ = mⱼ/(mⱼβⱼ² + ε)`, `v_g/(Σ_{k∈g} m_kβ_k² + p_g·ε)`, penalty diagonal
  `Dⱼ = mⱼwⱼ`. Outside it, the adaptive estimators would not be scale invariant
  (their weight map would see `βⱼ²/cⱼ²`), and the group-L0 fixed point would become
  `λ·v_g·⟨m⟩_g`, so `v_g` would stop being the group-L0 weight the `select_fit`
  Pareto front is built on. Invariance of the objective is not enough for a
  non-convex surrogate, so the IRLS cold starts are built from the metric too and
  the stopping rule is measured in metric coordinates `√mⱼ·βⱼ`, restricted to the
  penalized columns.
- An entry of exactly `0` marks a column **unpenalized**, and that is how the
  moment channel's μ₀ intercepts stop being shrunk — identically for all three
  estimators. A group weight could only have done it for the group form, leaving
  `fit(MomentFit, ds, Ridge(λ))` quietly shrinking the reference moment. Zero is
  reserved for a structural exemption (an intercept, an identically vanishing
  column); a numerically-zero estimate is refused, never floored.
- The unpenalized block must be well conditioned or the solve is refused by name:
  `X'X + λD` is positive definite exactly when it is, and a `Symmetric` solve on a
  singular matrix returns garbage rather than throwing.
- The reference ensemble uses a generator specified inside this package
  (SplitMix64 + the Archimedes sphere construction) rather than `Random`: the
  metric enters every penalized coefficient, so a stream that drifted between Julia
  releases would silently move every recorded penalized fit.
- Basis-aware constructors `Ridge(basis; ...)` / `AdaptiveRidge(basis; ...)` /
  `GroupAdaptiveRidge(basis; ...)` (and the `MomentBasis` forms) attach the metric
  **by default**, and the fitting doors refuse a metric whose provenance does not
  match the fit (wrong channel, wrong basis fingerprint, wrong `torque_weight`) —
  a mismatch is invisible to every numerical gate, since scale invariance holds for
  any `m ∝ c²`, right or wrong.
- Review-panel follow-ups (2026-08-24): `refit` carries the same provenance door as
  the other fitting entry points, and `AdaptiveLasso` follows its **pilot** for both
  the metric and its provenance (the pilot's coefficients set the weighted-L1 penalty
  factors, so a pilot under the wrong metric moves the whole solve).
  `SCEFitting.with_lambda(est, λ)` moves an estimator along a λ path with its metric
  intact — a hand rebuild drops it silently, and a dropped metric is indistinguishable
  from a deliberate uniform one. `AdaptiveRidge` gained the validating inner
  constructor it lacked, so a field-typed call can no longer install a negative penalty
  scale. `MetricProvenance` is a validating struct rather than a bare named tuple, and
  `show` renders it. The free-block conditioning guard now cuts where its message says
  it does (`κ(X_free) ≲ 6.7e7`, on the Gram), and the metric machinery moved to
  `src/fitting/metric.jl`.
- The default `nconfig` is 2048, sized from measurement rather than guessed. The
  relative standard error of `mⱼ` over eight independent seeds is 3.3 % median / 5.3 %
  worst column there, 6.0 % / 9.7 % at 512 and 1.6 % / 2.8 % at 8192 (bcc Fe 2×2×2,
  `lmax = 2`; body order barely moves it); the cost on bcc Fe 3×3×3 is 0.56 s there and
  5 s at 8192, linear in both `nconfig` and the column count. What is left is a wobble
  on the prior, an order of magnitude below the systematic factor the metric removes.
- **Breaking**: every penalized fit changes. `OLS` and `lambda = 0` are untouched
  (bitwise), as is any estimator constructed without a metric. λ recorded against
  an earlier penalized fit no longer means the same thing — bcc Fe `l02`…`l044` and
  the `m_*` moment fits, FeRh, and the KLM series need re-selecting, and
  `AdaptiveLasso(pilot = Ridge(...))` moves through its pilot. Same class of change
  as the `torque_weight = 0` assembly fix.

### Fixed — effective dof with unpenalized columns (2026-08-24)

- `_edof` (behind `effective_dof` / `gcv` / the `select_fit` GCV
  path) now handles a penalty diagonal with **exact zeros**, i.e. columns the
  estimator does not penalize. The previous code formed `X·D^{-1/2}` unconditionally
  and, on the cached-Gram path a λ sweep takes, divided the Gram by zero — producing
  a non-finite score rather than an error. The split is now taken before the `XtX`
  keyword is consulted.
- With `X = [X_F X_P]` (unpenalized / penalized) the effective dof is
  `rank(X_F) + Σᵢ sᵢ/(sᵢ + λ)` over the eigenvalues of
  `W^{-1/2}(X_P'(I − P_F)X_P)W^{-1/2}`, so an unpenalized column always costs its
  full degree of freedom and `df → rank(X_F)` as `λ → ∞`. The projector is applied
  through a thin QR of `X_F`, never formed `n × n`, and the eigenproblem still runs
  on the smaller of the primal and dual sides.
- A rank-deficient unpenalized block is **refused by name**: `A = X'X + λD` is
  positive definite exactly when `X_F` has full column rank, so a dependent
  unpenalized column leaves the fit unidentified and any finite dof reported for it
  would be meaningless.
- No behavior change for a fully penalized diagonal (every current in-tree
  estimator): that branch is untouched and byte-identical.

### Fixed — a self-image (`AllImages`) basis can no longer be fitted silently (2026-08-24)

- `SCEDataset` now refuses a basis whose SALCs carry a member using one
  reference-cell atom twice (an `AllImages` self-image pair `(a, 0)-(a, R)`) with
  `UnclassifiableBasis`, naming the offending keys and the way out. Both ends of
  such a pair carry the same spin on the reference cell, so the function collapses
  to a single-site one (a constant for `Lf = 0`); the per-orbit function-space
  reduction cannot see the collapse — its aggregate key treats the two factors as
  living on separate spheres — so the columns survived the build and the design lost
  rank silently. Measured on a one-atom cubic cell (a = 3.0, nbody = 2, lmax = 1,
  cutoff = 3.2): 27 SALCs, 9 columns identically zero and 3 constant, rank 5 of 27
  after centering, **no build warning**, and at solve time only the generic `OLS`
  rank warning (silent under a regularized estimator).
- The refusal sits at the dataset door, not the build: an `AllImages` basis is also
  the **tiling template** a downstream consumer expands onto a supercell, where the
  images become distinct sites and each self-image pair becomes a genuine bond
  (SCEMonteCarlo's cubic-Heisenberg tutorial). Building, `SCEPredictor`,
  introspection and export stay legal; only fitting on the reference cell is
  refused. To fit the same model, build on a supercell with `MinimumImage` — for a
  monatomic cubic cell, 3x3x3 with the same spec gives the same isotropic
  nearest-neighbor channel.
- That second contract is now documented where it is chosen: the `AllImages`
  docstring, the basis guide, `SPEC.md`, and `UnclassifiableBasis` (which now names
  both gates that raise it, and states the boundary of its contract: the **type** and
  the raising doors are public API, the **`reason` text is not** — never classify a
  failure by matching the message).

### Added — `[moment]` section in the TOML setup file (2026-08-24)

- `input.toml` takes an optional `[moment]` section (`nbody`, `lmax_mark`,
  `lmax_env`, `sampled`, `marked`, `cutoff_pair`, `cutoff_star`, `lsum`,
  `isotropy`; label tables and species-pair tables in the `[interaction]`
  spelling) that `read_setup` returns as `moment::Union{Nothing,MomentSpec}`,
  and `MomentBasis(path; backend, tol, tie_tol)` builds the pointed basis from
  the file the way `SCEBasis(path)` builds the energy basis. A fit script no
  longer carries the moment truncation as Julia constants next to a TOML file
  that carries the energy truncation. The reader only converts; every range /
  consistency rule stays in the `MomentSpec` keyword constructor. Unknown keys
  and the upstream spelling `soc` are refused. `[moment].isotropy` defaults to
  `true` (the `MomentSpec` default), unlike `[interaction].isotropy`.
- `include("io/input.jl")` moved after the moment basis (the reader now returns
  a `MomentSpec`); no behavior change.

### Added — `torque_weight_per_site` and a scale note on the co-fit objective (2026-08-23)

- `fit`'s objective `(1 − w)·MSE_energy + w·MSE_torque` measures the energy
  error on the **total** energy of a cell and the torque error per site
  component, so on a large cell a weight that reads torque-dominated is still
  energy-dominated (on a 1296-site cell `w = 0.99` leaves `0.01·MSE_energy`
  above `0.99·MSE_torque`; the torque error only moves near `w ≥ 0.999`).
  Found on Miyazaki's 36×36 Kondo-lattice data, where every `w ≤ 0.99` fit
  showed the same "tangential-field plateau". The `fit` docstring now says so,
  and `torque_weight_per_site(w_site, n_atoms)` maps a weight stated on the
  per-site energy scale to the `torque_weight` that realizes it
  (`w = w_site·n² / (w_site·n² + 1 − w_site)`). No change to `fit` itself —
  `lambda` scales of penalized estimators are untouched.

### Fixed — `moment_resolvability` no longer forms the dense signature block (2026-08-23)

- `_moment_resolvability` used to assemble the pointed signature expansion as a
  dense matrix `S` (rows = distinct signature keys, columns = pointed SALCs) and
  SVD it. On a supercell the row count scales as atoms × neighbours × harmonic
  components — tens of millions of rows on a 36×36 triangular torus — so
  `MomentDataset` (which runs the gate first) was SIGKILLed at 8–16 GB on every
  basis past ~200 columns there, before a single design row existed. The gate
  now groups the signature rows by the marked atom (blocks of different marked
  atoms share no key), assembles each block densely, accumulates the column
  norms, and folds the block into an upper-triangular `R` by a stacked QR
  (`R ← qr([R; B_a]).R`). `R = Qᵀ S` for an orthogonal `Q`, so `svd(R[:, kept])`
  carries the singular values and right singular vectors of `S[:, kept]` — rank,
  null combinations and the wide-block completion are computed exactly as
  before, and the vanishing test still sums duplicate (key, column) entries
  within a block before taking the norm. Peak memory on the 36×36 torus with
  301 pointed columns: 5.1 GB (basis 4.5 GB) instead of a kill; the existing
  gates (symbolic rank ≡ numerical design rank, null combinations annihilate
  the design, census, cache identity) are unchanged and pass. The
  `_MomentRowKey` lost its `mark_atom` field (the block index is the atom).

### Added — theory chapter on the pointed site-moment expansion (2026-08-22)

- `docs/src/theory/moment.md`: what the moment channel computes and why it is
  built this way — the adiabatic map, covariance versus invariance, the
  three-site star worked end to end, the mark as a decoration, Frobenius
  reciprocity (why the energy-side projector is reused unchanged), the site
  sum rule and the `l_mark = 0` "new information", the Landau origin of the
  `(0,1,1)` channel, the signed-projection regression with its axis rule and
  row gate, the induced-moment vector, scope and caveats, normalization
  bookkeeping, and a symbols-to-code table. Wired into the Theory section
  (`theory/index.md`, `make.jl`) and the home-page table.

### Added — Claude Code development procedure, aligned with Magesty.jl (2026-08-22)

The collaboration scaffolding Magesty.jl uses, ported with this package's
names and one deliberate difference (local commits on `main` need no
per-commit confirmation; remote operations still do):

- `Makefile` with the test tiers (`test-unit` / `test-aqua` / `test-jet` /
  `test-all` / `test-oracle` / `test-sunny` / `test-glmnet` / `test-pin` /
  `test-parity` / `test-examples` / `test-downstream` / `docs`), the CI-parity
  gate `test-ci`, the cold-start `ci-local`, and the `bench-*` targets. Every
  test target pins `JULIA_NUM_THREADS=4`.
- `.claude/agents/`: the two review tiers (`code-reviewer`; the four-axis panel
  `numerical-reviewer` / `maintainability-reviewer` / `performance-reviewer` /
  `api-reviewer`), `test-runner`, `profiler`, `spec-reviewer`, `git-helper`,
  `release-helper` — each tuned to this package's layers, hot paths, coupled
  sites, and test map. `.claude/hooks/no-japanese.sh` (PostToolUse; wired by
  the gitignored `.claude/settings.json`), `.claude/mcp-setup.md` + `.mcp.json`.
  `.claude/` is now tracked except the two settings files.
- `CLAUDE.md`: Core / Implementation / Language rules, the Makefile test table,
  Performance guidelines (hot paths + bench bookkeeping), the spec-folder
  workflow (`docs/specs/_template/`, `docs/specs/README.md`), design notes
  (`DESIGN_NOTES.md`, `docs/design-notes/README.md`), and the working
  principles (free / propose / confirm lists, sub-agent usage, review tiers).
- `CONTRIBUTING.md`, `SECURITY.md`, GitHub issue / PR templates,
  `CompatHelper.yml`, `TagBot.yml`.
- `bench/README.md` pointed at `bench/BENCH_LOG.md` (the `.claude/bench_log.md`
  path it named never existed here); one source comment no longer cites
  `CLAUDE.md` (scaffolding must not be referenced from `.jl` files).

### Added — docs: adiabatic site-moment guide; README links to the published site (2026-08-22)

- `docs/src/guide/moment.md`: the pointed-expansion guide (spec and basis, the
  resolvability gate, the `SpinDatum` trio and the readers that fill it, the
  dataset's axis rule and decomposability gate, fitting and prediction, the three
  diagnostics, scope and limits), executed at every docs build on a synthetic
  planted map. Listed on the home page and in `make.jl`.
- `README.md`: documentation and CI badges, the published site URL
  (<https://tomonori-tanaka.github.io/SCEFitting.jl/dev/>) in place of the stale
  "not yet deployed" note, and the moment channel in the status section.

### Added — `test/parity/`: real-data parity tier vs SLCE.jl (2026-08-21)

S6 of the pointed-moment backport. A separate environment (`SCEFitting` and the
sibling `SLCE.jl` as path sources, Manifest gitignored, **no CI job** — it needs
`~/Packages/SLCE.jl` and data outside the repository, resolved from
`SCE_PARITY_FEGE_DIR` / `SCE_PARITY_FEGE_POSCAR` / `SCE_PARITY_FERH_DIR` with
loud skips). Three kinds of statement, labeled in the file: acceptance numbers at
relative 1 % against the upstream protocol (FeGe `lsum2` 39 columns σ/CV/held-out
0.0398/0.0403/0.0308, `full` 181 columns 0.0334/0.0358/0.0295, FeRh gated σ
0.0086 μB), census integers as change detectors (FeGe kept 12757/12800 and
3838/3840; FeRh mode-1 `n_anti` 3833/7744 on the Rh orbit, 0/7744 on Fe, with a
hand recount `ê·e_MW < 0` from the raw data), and absolute column parity (design
columns matched by `SALCKey`, `‖a − b‖ ≤ 1e-10 ‖b‖` per column and elementwise
`1e-12·max|b|`; dataset targets/gates/masks bitwise; held-out predictions at
1e-10; both extxyz readers bitwise on one file). The FeRh constraint axes come
from an independent `M_CONSTR` parser written in the test. First run
(SCEFitting `9b95523` vs SLCE `e5a05c4`): every acceptance number within
±0.15 %, every census integer exact, and the column parity **bitwise**
(worst relative column deviation 0.0 on all three bases).

### Added — moment-channel diagnostics + mark-class shrinkage (2026-08-21)

Step M5 of the pointed site-moment backport: upstream SLCE.jl `bd03517`
(`MomentDataset.order` + `moment_band_profile`), `a6a6386` (mark-class
`salc_groups(::MomentBasis)`, `GroupAdaptiveRidge(mb; lambda)`,
`_reduce_to_active`, the summarized dependency log, the copied dependency
record) and `3d54abb` (`_moment_axis_matrix`, `_pair_neighbors`, `_legendre`,
`moment_local_field`, `moment_coverage`, `moment_simple_floor`), adapted to
`SpinDatum` (no provenance, 2-arg `MomentModel`, `isotropy` for `soc`).

- **`MomentDataset.order`** — per-config marked-sublattice order parameter
  `|⟨e⟩| = ‖Σ_a e_a‖/n_marked`; **`moment_band_profile(model, ds; nbins)`** /
  `(f; nbins)` — per-config mean residual over the KEPT rows in equal-count
  `|⟨e⟩|` bins plus the bin-free least-squares line and Pearson `r` (the L2-2
  basis-insufficiency signature, to report next to any σ).
- **`salc_groups(mb::MomentBasis)`** — group labels keyed by
  `(body, orbit_id, decors, marked atoms, marked sites)` of the canonical
  representative member: the energy-side key folds an Fe-marked and a Ge-marked
  placement of one pair orbit (they differ only in `block`), and the atom set
  alone folds two mark placements on a member carrying two periodic images of
  one atom. **`GroupAdaptiveRidge(mb; lambda, …)`** with unit weights;
  **`fit(MomentFit, …)`** reduces a `GroupAdaptiveRidge` to the active columns
  with the vanishing freeze (`_reduce_to_active`: emptied groups relabeled
  away, weights follow; a label vector of the wrong length is refused). The
  stored `estimator` stays the caller's un-reduced object.
- **`moment_local_field(mb, configs; axes = configs)`** / `(mb, data)` — per row
  `‖h₁‖` and `ê·ĥ` with `h₁ = Σ_j ê_j` over the marked atom's `cutoff_pair`
  MinimumImage neighbors (one term per tied image, same-atom images excluded,
  `lmax_env = 0` species excluded; tie band = `mb.tie_tol`); the `SpinDatum`
  method resolves the axis by the mode rule through `_moment_axis_matrix`, the
  ONE function the dataset constructor reads. **`moment_coverage(train, new;
  q)`** — upper-tail `h1` threshold + `frac_beyond` + the anti-alignment
  fraction `frac_anti` (the measured collapse coordinate on FeGe).
- **`moment_simple_floor(f, data; lmax)`** — the nested simple-feature floor
  `y ≈ μ_g + Σ_l b_{g,l} Σ_j P_l(ê_i·ê_j)` on exactly the fit's kept rows, with
  `sigma_floor` / `sigma_model`, a per-feature `inclusion` (relative projection
  residual onto the kept design's range, rank-cut SVD — never `Xk \ F`),
  `nested_bound = f.estimator isa OLS`, and a two-half pairing door (bitwise
  target replay + one config's design-row replay) that refuses re-paired data.
- The dependency warning summarizes past 8 combinations (wide P1 blocks); the
  dataset stores a COPY of the cached resolvability record.

Gated in `test_momentfit.jl` (204 → 403): hand-oracle `order`, the profile
recomputed from public fields + an independent normal-equation line + a planted
linear-in-order corruption recovered by the slope (+0.5), fewer-configs-than-bins,
the different-basis refusal; `salc_groups` refining the energy key with a
design-side disjoint-row-support oracle and ≥ 2 groups per pair orbit, gauge
blocks folding (same row support), the same-atoms/different-sites split on a P1
image cell; the GAR freeze reduction with `lambda = 0 ≡ OLS` and a non-uniform
weight relabeling; the **real face-(a) vanishing fixture** (3×3×6 cell,
`Δf = (.5, .5, .25)`, ops `[I, m_y]`, `isotropy = false`: 16 of 38 columns vanish
— exact zero on random periodic data, null report complete, end-to-end freeze
with no ctor injection, the GAR path on the real set) and its P1 face-(b)
control (every reported combination annihilates the numerical design; the
summarized log); the local field on a hand-derivable ±x tied cell (`h₁ = 2ê_other`,
`ê·ĥ = cos θ`, mode-1 `x̂` axis → `sin θ`, zero axis → NaN, `lmax_env = 0` species
off), a hand-checked coverage quantile (`9.901`), and the floor recovering planted
`(a₀, b₀)` with `P₂` reported non-representable on an `lmax_env = 1` basis, plus
every pairing door.

### Changed — `MomentDataset` doors closed after the M4 numerical + saboteur reviews (2026-08-21)

- **The `directions` door**: `MomentDataset` now validates every datum's
  `directions` with the family's unit-column rule (`_validate_config`: finite, 1e-6
  norm band, component bound). A `SpinDatum` built field-by-field carries no
  direction check, and nothing upstream of the dataset did either — a non-unit
  column corrupted the gate's `M⊥` (so `g` was not `|M| sin²θ`) and fed the Legendre
  recursion outside `|z| ≤ 1`. The zero-axis exclusion is now explicitly mode-1
  only (a mode-4 axis IS a validated direction; an exactly-zero `directions`
  column is refused, never silently `defined = false`).
- **The zero-moment placeholder door** (`zero_moment_atol = 1e-10`, the moment
  channel's analogue of `SCEDataset`'s): a referenced atom (marked, or an
  environment site of any pointed member) with `‖MW‖ ≤ atol` is refused by name —
  its `directions` column is the ẑ placeholder the moments constructor fabricates,
  which would enter the design as a fake coordinate while the `|M| = 0 → g = 0`
  convention waved the row through. `M_int = 0` on a marked atom is NOT this case
  and still passes. The applicability-limit paragraph now states the remaining
  soft-environment half precisely.
- **Door order**: `moment_resolvability` runs first (basis-only, cached), so an
  unclassifiable basis surfaces its own refusal before any data door can mask it.
- **Antiparallel rows are logged** (`@info`, mode-1 `ê·e_MW < 0`) with the exact
  lossless remedy — re-gauge `ê → sign(ê·e_MW) ê`, `y → −y` at the source — instead
  of a bare count in the report; the one-setup / one-reference obligation the
  datum cannot carry is stated in the docstring.
- `_design_moment`'s shape checks moved out of the threaded loop (a throw there
  surfaced as a `TaskFailedException`), and its precondition comment names the
  actual doors.

New gates (`test_momentfit.jl`, 138 → 204): the `directions` door on marked and
environment columns in both modes (scaled / near-pole / zero / NaN), the
zero-moment door on a marked Fe, a sampled-environment Ge, and an unreferenced Ge
(unsampled basis: passes), two literal hand-worked rows (`ê = ẑ`, `M = (0.1, 0.2,
3.0)` → `y = 3`, `g = 0.05/√9.05`; `ê = x̂`), the two-orbit μ₀ shift (`c0` on Fe
rows, `c1 ≠ c0` on Ge rows, residuals invariant; the two `[MARK]` columns are
exact orbit indicators), the normal-equation oracle and `min|r| ≤ rmse ≤ max|r|`
on the ungated solve, a `groups` spy estimator (receives exactly
`row_config[keep]` then `row_config[defined]`), the Ge-only marked basis with
mixed modes (`row_atom == repeat(5:8, 2)`, the undefined row at `nm + 3`, Fe
`directions` as environment coordinates), the exact survival ratio
`(20nm − 4)/20nm` with a `coverage_floor = 0.96` refusal naming the undefined
count, and `predict_moment` refusing a bad `e` behind good `axes`.

### Added — `MomentDataset` / `fit(MomentFit, …)` / `MomentModel` / `predict_moment` (2026-08-21)

Step M4 of the pointed site-moment backport (upstream SLCE.jl `837135d`, adapted to
the pure-spin `SpinDatum`: no provenance, so no setup-uniformity / reference-identity
doors; no displacements, so no reference-geometry door — every datum sits at the
reference by construction).

- **`MomentDataset`** (`fitting/momentfit.jl`): rows `(config, marked atom)` with
  `y = ê·M` under the mode rule (4 → `directions`, identity substitution; 1 →
  `constraint_axes`; mixed allowed; a mode-1 zero-axis marked atom is
  `defined = false` with `y = NaN`, excluded from every fit). Decomposability gate
  `g = ‖M⊥‖²/|M| = |M| sin²θ ≤ gate_eps` (required keyword, cancellation-free form,
  `|M| = 0` passes) with per-orbit survival, rms `‖M⊥‖`, the antiparallel census
  `n_anti`, and a `coverage_floor` refusal before the design build.
  `moment_resolvability` runs at the door: unclassifiable refuses, vanishing
  columns are recorded (frozen to exact zero by `fit`), dependent combinations are
  warned and disclosed.
- **`fit(MomentFit, ds, estimator = OLS())`** solves the gated rows and, for
  disclosure, the ungated ones (`groups = row_config`, no centering — the `l = 0`
  `[MARK]` columns are the per-orbit intercepts μ₀); `residuals` / `rmse_moment`
  read either solve.
- **`MomentModel` / `predict_moment(model, e; axes = e)`**: a validating door
  (`e` unit everywhere; `axes` unit and component-bounded on the marked columns
  only — unmarked columns are never read) whose default `axes = e` is the mode-4
  identity.
- **Persistence is refused by name** (design record §4.2, decided here): `save`
  on a `MomentBasis` / `MomentFit` / `MomentModel` throws and writes nothing; a
  later schema version adds the moment side. Upstream cannot save it either.

Gated in `test_momentfit.jl` on the FeGe B20 fixture: hand-arithmetic targets and
gates under both modes, mode-1 ≡ mode-4 identity, planted-model recovery and
held-out prediction equivalence, the gate's keep/reject disclosure with a
hand-oracle transverse rms, the `|M| = 0` pass, zero-axis exclusion (and that the
caller's datum is not mutated), the coverage-floor refusal naming the orbit,
bitwise time reversal in both modes, the per-orbit μ₀ absorbing a constant shift,
the vanishing-column freeze mechanism, the prediction door (including the
component bound), multi-orbit bookkeeping, the hard refusal of an unclassifiable
basis at the dataset door, and the persistence refusal.

Real-data acceptance (design record §3.3, the "right after M3" item): on the
2×2×2 FeGe B20 sets SLCE.jl wrote (τ0.1 + τ0.5 training, τ0.3 held-out, gate
2.2e-4 μB), read here through `read_extxyz(...; reference)`, the L1 pair basis
(13 columns) reproduces upstream's held-out σ = 0.0699 μB to +0.00 % and the L2
pair basis (25 columns) 0.0824 μB to +0.03 %, with identical kept counts
(12757/12800). Script and log: the design record's `step2_assets/fege/`.

### Added — `MomentBasis` gates closed after the saboteur and numerical reviews (2026-08-21)

Every fixture had `marked_atoms == 1:8`, so a slot or row index used as an atom
number was invisible, and only two orbits were ever named. New gates: a Ge-only
marked basis (`marked_atoms == 5:8`) with locality in atom numbers; the
configuration-major row contract; the M3-1 species rule on a real basis
(`lmax_env = [2, 0]`: no environment spin ever sits on a Ge atom, Ge is still
marked); `_moment_labels` against a hand enumeration of its four rules, with and
without `lsum`, plus an `lsum`-capped basis; `isotropy = false` on the 3-body
star basis containing the isotropic keys; the census with exact counts (nn Fe–Fe
4, Fe–Ge 8, Ge–Ge 4, and 4 with Ge unmarked); rank–nullity on the reported null
space and `rtol` monotonicity; vanishing columns against the random design; and a
P1 chain whose re-anchored triangle has a non-minimum-image environment edge, so
only the anchor may carry the mark — the `(0,1,1)` column is exactly zero on the
other two rows.

The numerical review then separated oracle from pin in the FeGe gates: the
geometric enumeration is the independent reference and is now asserted
**scale-free** (the ratio of the SALC column to the geometry sum is one constant
across atoms and random configurations — what caught the missing-orderings bug
upstream); the constants `6.0` and `2√3` were read off a run of the upstream
prototype, so they are pinned separately and labelled as change detectors. Two
doors tightened: `_mark_term_index` asserts the one-mark invariant the whole
substitution argument rests on, and `_design_moment` validates the shape of every
axes matrix. The census docstring now says what the code computes (distinct marked
atoms per orbit, not stabilizer-inequivalent placements — the `≥ 2` reading was a
false positive with one species unmarked). The repeated-image star members are
kept, as upstream keeps them (column parity on small cells); the enumeration site
says so, and the dataset door's `UnclassifiableBasis` refusal on such a basis is
now a test.

### Added — `MomentBasis`: the pointed SALC basis for adiabatic site moments (2026-08-21)

Step M3 of the pointed site-moment backport. Ported from upstream SLCE.jl at
`3d54abb`, which is `6960276` (the basis) plus two later commits this file
carries with it: `bd03517` (the mark→term index fast design path, asserted
value-identical to the full evaluation) and `28dfa24` (the `moment_resolvability`
result cache **and a bug fix**: the null-combination report of a WIDE signature
block — more kept columns than signature rows — came back empty upstream because
the economy SVD lists only `min(r, c)` directions; the orthogonal complement is
now read off a QR completion). One spelling change: the screen is
`MomentSpec(; isotropy = true)`, this package's name for upstream's `soc = false`.
The review panel caught the citation (the first draft said `6960276` alone).

The moment channel's counterpart of `SCEBasis` (`src/basis/momentbasis.jl`): per
marked reference-cell atom the design row is the pointed SALC vector, and one
shared coefficient vector serves every symmetry-equivalent site. The MARK is the
displacement decor `SiteDecor(disp = (1, 0))` (the polar `|u|²` factor evaluated
on an indicator field), so the existing decor engine, evaluation kernels, and
canonical-member machinery carry the basis unchanged through the `admit` kwarg
the engine already had.

- **`MomentSpec`**: mark-aware truncation — the marked site's ê factor is
  allowed for ANY species (`lmax_mark`; an E-inactive species' induced moment is
  what the channel predicts), environment spin factors only for species the
  consumer samples (`lmax_env` + the required `sampled` claim, refused loudly on
  mismatch — M3-1 decision A). Time reversal keeps even-Σl labels only (the mark
  rank counts).
- **3-body stars are mark–environment-bond cut** (M2-5): both mark bonds
  minimum-image within `cutoff_star`, the environment–environment edge free (the
  triangle is pinned by the mark bonds). The enumeration expands every star to
  all 3! site orderings — `candidate_clusters`' multiplicity convention, without
  which the closed-star column came out at half the prototype's 6.0 geometric
  oracle upstream. `build_clusters`' orbit-formation core is factored into
  `_orbits_from_members` (behavior identical, asserted) so the star enumeration
  reuses it.
- **`_design_moment`** evaluates rows `(config, marked atom)` with the
  marked-COLUMN substitution: the spin matrix's marked column is replaced by the
  evaluation axis (identity for mode 4), exact because every pointed label
  carries exactly one mark; a mark→term index skips the dead terms and is
  asserted value-identical to the full evaluation.
- **`moment_resolvability`** (D9′): symbolic signature expansion in the
  independent variables `(a, ê_a, e)` → vanishing columns, numerical rank, and
  null combinations that NAME the dependent columns; plus the mark-class census
  per cluster orbit (face-(b) hazard preregistration). Members putting two
  environment spin factors on one reference-cell atom (two periodic images of
  one neighbor) are refused as `UnclassifiableBasis` (defined here; upstream's
  energy-side gate shares the name) — the monomial signature would overcount
  the rank there.

Gated in `test_momentbasis.jl` against the design-record prototype's independent
geometric references: star closed form = 6.0 × the geometry sum, shell-sum
normalization 2√3, G_i covariance including the axes, bitwise time reversal,
substitution locality, signature rank ≡ independent random-design rank with null
combinations annihilating the actual design, and `isotropy = false ⊇ isotropy =
true` on the keys.

### Added — extended-XYZ training container, axis gates, EMBSET pair reader (2026-08-21)

Step M2 of the pointed site-moment backport (upstream SLCE.jl `d2f9d2f`, adapted to
the pure-spin `SpinDatum`).

- **`write_extxyz` / `read_extxyz` / `ExtxyzFile`** (`io/extxyz.jl`): the canonical
  on-disk format for new constrained-noncollinear training sets, in SLCE.jl's
  dialect so spin-only files interchange between the packages unchanged. The
  structure is always stored (self-containment removes the "which POSCAR pairs
  with which EMBSET" bug class); per-atom columns `mw` / `bcon` / `mint` /
  `mconstr` map 1:1 onto the datum's channels; numbers print shortest-round-trip,
  so every stored value survives bit-exactly (directions / magmoms / torques
  re-derive from the written moment vectors, exact to rounding). **This package
  refuses a joint file by name** — displaced frames, `forces` columns, or a
  `config_type=joint` claim all error pointing at SLCE.jl — and never flattens one
  to its spins; spin-only vs joint is measured from positions, the claim only
  cross-checked. Upstream's provenance keys are accepted and ignored.
- **`check_moment_gates`** (public, unexported): the moment channel's
  axis-consistency gates — mode-1 sign consistency (`sign(ê_MW·ê_c) == sign(y)`
  on rows with `|y| > 5e-3 μ_B`) and the axis-angle 99th percentile (`< 5°`) —
  run at extxyz generation, at every extxyz load, and in the pair reader, so
  archived constraint axes are re-verified, never believed. Ported unchanged.
  Pinned subtleties: a whole-axis flip in mode 1 is a gauge and must not fire;
  the angle gate is a percentile because collapse rows carry large single-row
  angles.
- **`read_embset_pair`** (`io/embset.jl`): the legacy `EMBSET` (smoothed `MW`) +
  `EMBSET_mint` (bare `M`) sibling reader with loud pairing checks (config
  count, block shape, field blocks bitwise); energy lines deliberately uncompared
  (upstream measured ΔE = 0.148 eV between the two writers on the FeRh archive).
  `read_embset` now shares the parsing core (`_read_embset_blocks`) and is
  otherwise unchanged. A pair whose moment blocks are all bitwise equal (the same
  file twice, a byte copy) is refused — it is not an MW / M_int pair.
- **Hardened after the saboteur review** (stricter than upstream's reader, which
  the same review notes shares the misreads): a repeated key in the info line,
  a repeated property name, and a second string column are refused instead of
  letting the last value win; the writer quotes free-text values that contain
  whitespace and refuses a double quote or a line break inside one (the lexer
  has no escape), so the file it writes always reloads; the gate docstring
  states the exact tolerance (`⌊n/100⌋` rows) and that the thresholds are
  unbounded knobs. Tests now pin the mode-4 antiparallel axis (180°, no gauge
  there), the `ceil` percentile at 99/100/101 rows, the zero-axis skip at the
  gate, every cross-frame check with a violating frame, a triclinic cell under
  the reference check, the trio-less writer header, keyword forwarding through
  `ExtxyzFile`, and per-config axes in the pair reader.

### Fixed — `constraint_axes` carries the component bound (2026-08-21)

Review of M1: the axis validation had the `1e-6` norm band but not the
`maximum(abs, u) ≤ 1` component bound every other direction door of this
package applies (`_validate_config`, `Harmonics._validate_unit`) — and its test
pinned a column with `u_z = 1.0000005` as legal. An axis is a direction the
moment channel will hand to the harmonic kernels, whose `dnPl` domain is
`|z| ≤ 1`; a near-pole column inside the norm band would throw a bare
`DomainError` from inside an accumulation. The bound is added, the band reuses
the package's `_DIRECTION_ATOL` instead of a second constant, and the test
refuses the near-pole column and accepts an in-band, in-bound one. Same review,
smaller: under mode 1 an all-zero axes matrix is refused (it satisfies "axes
present" in letter only) and the per-atom zero-column convention is written
where the invariant lives (no moment row for that atom, never a fallback to
`directions`); `constraint_mode = true` is refused (a `Bool` is not a class);
the five-argument direct form converts its arguments again, as the default
constructor it replaced did; the exactly-zero column is stated as this
package's convention rather than as a VASP fact.

### Added — `SpinDatum` carries the adiabatic-moment channel trio (2026-08-21)

Step M1 of the pointed site-moment backport (upstream SLCE.jl `4b611f9`, adapted:
this package's datum is the five-field `SpinDatum`, not upstream's
`TrainingDatum`). Three optional fields, `nothing` when absent:

- **`moments_bare`** — the bare `M_int` vectors, the projection target
  `y_a = ê_a · M_a` of the moment channel. Validated for **finiteness only**: the
  signed readout is what keeps the target analytic where `‖M‖ → 0`, so no
  magnitude or sign rule is imposed. Distinct from `magmoms · directions`
  (the smoothed `MW_int` the constraint acts on); both are stored.
- **`constraint_axes`** — unit columns or exactly-zero columns ("no axis");
  anything in between (near-zero noise, an off-unit axis outside the `1e-6`
  band upstream uses for every direction) is refused, never normalized.
- **`constraint_mode`** — `1` (transverse-penalty type) or `4`
  (direction-pinning type). **The evaluation-axis rule is keyed by the mode,
  never by field presence**: mode 1 requires the axes, axes without a mode are
  refused.

The moments constructor passes the trio through as keywords; the five-argument
direct form leaves it absent, so every existing construction reads unchanged.
Gates: both construction paths carry the fields identically, the derived E/T
fields are bit-identical with and without the trio, and `SCEDataset`'s design
matrices and targets are `==` with and without it.

### Fixed — the persist reader validates a term against its slots (2026-08-21)

`_term_from` accepted any slot list alongside any tensor: a term with more
slots than tensor axes, a slot addressing a site the member does not have, or
an axis whose extent is not `2l + 1` was built as is (`SALCTerm` has no inner
constructor) and failed later, inside a kernel, as a `BoundsError` or a
silently truncated contraction. The reader now refuses all three with a named
`ArgumentError`, on both the v5 `slots` and the v2–v4 `ls` spellings, and
`_member_from` hands it the member's site count. `main`'s reader has the same
two gaps for slot/rank and axis extent (its `ls` terms could mis-state the rank
or the extent just as freely); the site range is new to this branch, since an
`ls`-addressed term could not name a site. Upstream SLCE.jl has the identical
reader and receives the identical patch.

The review panel then found the one **silent** member of the family: nothing
checked a term against its **key**. A DISP slot under a pure-spin key passed
every per-term check and then slipped past the spin-only kernels' refusal,
which reads the key's `decors` — the DISP axis was evaluated as a spin harmonic
under the wrong `(4π)` scale and a plausible wrong energy came back.
`_salc_from` now requires `body == length(decors)`, every member to have `body`
atoms, and every term's slots to reconstruct the key's decoration label exactly
(`_term_decors`; two factors of one channel on a site, or a site with none, is
refused too). Load-time only; same patch upstream.

### Added — the mixed-channel (decor) SALC engine (2026-08-21)

A second projection engine, `_orbit_salcs_decors`, alongside the pure-spin
production one. It takes **explicit decoration labels** (sorted `SiteDecor`
multisets) on a cluster orbit instead of enumerating `l`-tuples, which is what
the pointed site-moment channel needs: its mark is a displacement decor, so its
labels are not expressible as an `ls` tuple. Nothing calls it from a public
builder yet — the production pure-spin path is untouched.

- **Construction**: multiset arrangements grouped into site-permutation orbits;
  spin-first slot coupling with the total spin rank `L_S` read off each coupling
  path (`_path_LS` — a good quantum number, since site permutations act within
  channels and commute with the diagonal rotation); projection per `(L_S, Lf)`
  block; transport in canonical slot order; the `Σl_spin`-even time-reversal
  screen. Both channels rotate through the one polar Wigner cache, which is
  exact **because** of that screen (`det(R)^{Σl_spin} ≡ +1`, so the axial spin
  action equals the polar one).
- **`isotropy` is the decor engine's spelling of the pure-spin screen**: there
  it keeps `Lf == 0`, here `L_S == 0`. On a pure-spin label the two are the same
  statement (`L_S ≡ Lf`); on a mixed label neither implies the other, so
  `AngularMomentum.build_real_bases` gained a `keep` path predicate and the
  screen is applied **before** a rejected path builds its tensor. It is a
  **required keyword** of the decor engine (review finding P1): upstream
  SLCE.jl's engine takes `soc::Bool` in that positional slot with the opposite
  polarity, so a verbatim upstream call compiled here and silently inverted the
  screen; now a positional copy is a `MethodError` and an omitted screen an
  `UndefKeywordError`, both asserted.
- **`admit`**: an optional per-orbit-of-assignments predicate. A label is a
  sorted multiset and cannot say which site may carry which rank; the pure-spin
  engine expresses that through its per-species `lmax`, and a caller that needs
  the same rule (or the pointed basis' mark-aware variant) hands it in. It is
  applied at the canonical representative, before any `block` index is consumed.
- **Joint evaluation** `evaluate_salc(salc, e, u)`: spin axes contribute
  `Z_{lm}(ê)`, displacement axes `|u|^{2k} R_{lm}(u)`, and the scale is
  `(4π)^(n_spin/2)` over the spin slots. A pure-spin SALC evaluates `===` to the
  two-argument form; a decorated one is exactly `0` at `u = 0`.
- **The spin-only kernels now refuse a decorated SALC** — `evaluate_salc(salc,
  e)`, `accumulate_grad!`, `group_costs`, and the function-space reduction
  (`_function_vector` / `_reduce_orbit_salcs`, whose aggregate key reads SPIN
  ranks only and would merge or drop decorated terms silently; its message is
  the checklist for wiring the decor engine in) — rather than reading a `DISP`
  rank as a spin harmonic under the wrong `(4π)` scale. Same refusing-beats-
  mis-scaling rule as `multipole_terms`. Allocation-neutral on the
  `bench_salcbasis` gate (73,263,913, unchanged).
- **The anti-drift gate** is the point of the slice: given the same label and
  the same admission rule, the decor engine must reproduce the production
  engine's SALCs **bitwise** — keys, `block` indices, slots and folded tensors.
  It runs over the O_h single site, the D4h bond, and the two shapes an upstream
  review found broken: the Cs isosceles triangle (`ls = [1,1,2]` splits into two
  ordering orbits sharing one sorted label) and the C3v triangle (one orbit of
  three assignments, where the gauge is column-order dependent).
- **Mixed labels are checked under a genuine 3-cycle**: the mixed invariance
  gate had only the single site (identity permutation) and the bond (an
  involution), so a slot relabelling that is wrong under a cyclic site map — the
  case the C3v anti-drift covers for pure spin only — would have passed. The C3v
  triangle now runs the pointed shape (mark on one site, ranks on the other two)
  and an asymmetric three-way label through all six ops with their site maps,
  after asserting the stabilizer really contains a 3-cycle and that the rotation
  without its site map is not a symmetry.
- **The invariant counts are checked against a Cartesian projector** that shares
  no code with the SALC machinery — no Clebsch–Gordan, no Wigner-D, no spherical
  harmonics — built by averaging the group action over the multilinear forms in
  the components of `ê` and `u` and reading the rank.
- **The joint kernel's value is pinned by a closed form** (added after the
  review panel found that none of the gates above fixed one: the twist
  `(ê₁×ê₂)·(u₁×u₂)` is symmetric under `e ↔ u` and the site swap, `u = 0` is
  exact for any homogeneous kernel, and invariance holds at any scale — a kernel
  that dropped `|u|^{2k}` or scaled by the slot count instead of the spin count
  passed everything). Under the trivial group a label's SALCs are an orthonormal
  basis of the full product space, so by the addition theorems
  `Σ_s Φ_s² = (N!)² · Π_spin (2l+1) · Π_disp |u_a|^{2l+4k}` — no folded tensor,
  Clebsch–Gordan coefficient or Wigner matrix enters; the `(N!)²` is the
  documented `_canonicalize_members` convention (the `N!` ordered images of one
  instance are summed). Six labels, including the bare pointed mark
  (`n_spin = 0`, value `|u|²` exactly) and mark-and-rank on different sites; a
  per-SALC contraction written from the public `Zlm`/`Rlm` sits beside it.
  Mutation-checked: each of the three escapes is killed by this testset alone.
  A second panel then closed what this gate still let through: a spectator-atom
  fixture whose bond atoms are `[2, 3]` (a slot's site index used as its atom
  number was invisible on every `atoms == 1:N` fixture), the total spin rank
  `L_S` and the `isotropy = true` subset asserted on every asymmetric label,
  the refusals exercised on the pointed label shape (where one decor *is* pure
  spin, so an `any`-shaped guard would admit it), and a two-DISP-slot label
  under the C3v 3-cycle.
- **Performance**: the basis build is allocation-neutral to the byte; the design
  matrices cost exactly +92 allocations each, fully attributed to the one new
  `SALCScratch` field over the bench's 46 columns. See `bench/BENCH_LOG.md`.

### Added — the `SolidHarmonics` displacement kernel (2026-08-21)

A pure addition: no existing code path calls it yet. `SolidHarmonics` is the
counterpart of `Harmonics` for polar factors — real solid harmonics `Rₗₘ(u)`
evaluated as homogeneous polynomials in the Cartesian components of `u`, so
they are regular and exactly `0` (for `l ≥ 1`) at `u = 0`. Its normalization
is 4π-free (Racah-type):
`Rₗₘ(û) = √(4π/(2l+1)) · Harmonics.Zlm(l, m, û)` on the unit sphere, and the
rank-1 factors are literally `R₁₋₁, R₁₀, R₁₁ = y, z, x`.

Why a pure-spin package carries a displacement kernel — stated precisely,
because the short version overstates it. The **mark** of the pointed site-moment
channel is the displacement decor `SiteDecor(disp = (1, 0))` evaluated on a
synthetic indicator field, so the mark factor is `|u|² R₀₀` — 1 on the marked
atom, 0 on every other. That makes the mark arithmetic rather than procedural
(a member not marked at the row's atom is killed by an exact zero before it can
contribute), which is what keeps the selection independently checkable. But that
factor is the **only** one the pointed enumeration ever builds: the marked site's
angular rank lives in its SPIN part, environment sites are pure spin, so every
`DISP` factor in a pointed label is `(k = 1, l = 0)` and the kernel is read at
`R₀₀ ≡ 1`. The production channel needs `|u|²`, not a harmonic expansion.

What genuinely needs `l ≥ 1` is the **verification layer**: the decor engine's
own gates evaluate displacement-decorated SALCs at `l = 1, 2` — the chirality
twist against its closed form, mixed space-group invariance under the axial
spin / polar displacement action, and the exact zero at `u = 0`. Those are what
pin the engine's rotation convention, `(4π)` scale, and slot order independently
of the pure-spin engine, and they are the reason to carry the module rather than
special-case `l = 0`. Basis *construction* never calls it at all (`l = 0` axes
skip rotation, `l > 0` axes go through the Wigner cache), and no physical
displacement field enters this package's models.

The value kernel is ported from SLCE.jl `0ba2dc5`; its tests are independent
of it (the Racah relation against `Harmonics.Zlm` up to `l = 16`, homogeneity
`R(λu) = λˡ R(u)`, the closed-form rank-1 and `u = 0` values). **The gradient
API upstream ships alongside it (`solid_harmonics_grad`, `solid_harmonics_grad!`,
`grad_Rlm`, and the `∂A/∂z`, `∂A/∂r²` recurrence state inside the evaluator)
was dropped before merge**: a pure-spin package has no force rows, the
production route reads the kernel at `R₀₀ ≡ 1` alone, and the joint gradient
form refuses a decorated SALC, so it had no caller and no prospect of one. The
value recurrence is the upstream one line for line.

### Changed — SALC terms carry slots, not a per-site `ls` (2026-08-21)

**Breaking** for anything that read `SALCTerm.ls` or built a `SALCTerm`
positionally. `MultipoleTerm` is unchanged, so the downstream Monte-Carlo
contract is untouched.

- **`SALCTerm.ls::Vector{Int}` -> `SALCTerm.slots::Vector{Slot}`**: one entry per
  tensor axis, carrying the member-site index it contracts against plus its
  decoration factor. A pure-spin term is the identity slot list
  `spin_slots(ls)`, so every moved path is bit-identical on the bases this
  package builds.
  - `_canonicalize_members` remaps slot sites under the `(atom, shift)` sort and
    brings the axes to the canonical slot order — which reduces exactly to the
    old axes-follow-sites `permutedims` for identity terms.
  - The evaluation and gradient kernels read `l` and the atom through the slots.
  - `group_costs` and the orbit reduction key on `_term_spin_ls(t)`, the term's
    SPIN-axis ranks — the same values the per-site `ls` gave.
- **Persistence**: v5 term docs store `"slots"` (`[site, channel code, k, l]`
  per axis); a v2-v4 term's `"ls"` maps to the identity spin slot list on read.
  The gate-(a) test now downgrades both halves of a document.
- **`multipole_terms` refuses a displacement-decorated basis** rather than
  dropping or mis-scaling: a `MultipoleTerm` has no displacement factor, and
  consumers derive the `(4π)^(body/2)` scale from the term shape.
- **Verified by a cross-commit bit-identity probe**, not only by the suites: a
  model built and serialized by the pre-refactor commit reproduces byte-equal
  keys, coefficients, energies, torques, fingerprints and `multipole_terms`
  both when rebuilt on the new code and when its file is read through the new
  reader (5 crystals, 643 compared lines, both arms identical).
- **Performance**: the design-matrix path is allocation-neutral to the byte; the
  basis build costs +3.4 % allocations and +3.8 % wall time for the wider label.
  Numbers, attribution and one rejected micro-optimisation are in
  `bench/BENCH_LOG.md`.
- `test/unit/test_selection.jl` gains a direct assertion that on a pure-spin
  basis a term's slot count IS the body order — the identity the Monte-Carlo
  slot pricing rests on, and one that byte-equality of the basis cannot catch.

### Fixed — `SALCKey`'s hash is content-based again (2026-08-21)

Adding `SiteDecor` to the key silently made `SALCBasis.fingerprint` a function
of the **build** rather than of the basis: `hash` on a bare struct or enum falls
back to `objectid`, which mixes in the type's identity, and that identity
differs between precompiled images of the same package. Measured: three
environments of this package produced three different hashes for the same
`SiteDecor(; spin = 1)`, stable across processes within each image.

- `hash(::SiteFactor)` and `hash(::SiteDecor)` now hash the content tuples,
  projected to plain `Int`s (`Channel` is an enum, so it hashes by `objectid`
  too and has to go through `Int`).
- Nothing was functionally broken — `_basis_from_doc` already recomputes the
  fingerprint locally and treats the stored one as provenance, and every
  fingerprint comparison is between two in-session objects. What was lost was
  the property the docstring claims (a *structural* fingerprint) and the ability
  to compare a basis across builds, which is how this was found: a cross-commit
  bit-identity probe reported identical keys, coefficients, energies, torques,
  and multipole terms — and a different fingerprint.
- The fingerprint remains Julia-version dependent (`hash` of `Int`s and tuples
  is), exactly as `_basis_from_doc`'s docstring says.

### Changed — `SALCKey` carries decorations, persistence schema v5 (2026-08-21)

**Breaking** for anything that read `key.ls` or constructed a `SALCKey`
positionally; the numerics are untouched (the pin tier is byte-identical across
this commit, which is the evidence that the relabel is value-preserving).

- **`SALCKey` is now `(body, orbit_id, decors, L_S, Lf, block)`.** `decors` is
  the sorted `SiteDecor` multiset — the generalization of the sorted `ls` label —
  and `L_S` the total coupled **spin** rank. The pure-spin construction emits
  `decors = spin_decors(ls)` with `L_S = Lf`, which is the total,
  value-preserving v4 → v5 map. `SALC` mirrors `decors` / `L_S`.
  - `spin_ls(key)` reads the old `ls` back; `is_pure_spin(key)` says whether
    every decor is a bare spin factor. Both are `public`.
- **Persistence schema version 5**: keys store `"decors"` (per-site
  `[spin_l, disp_k, disp_l]` triples, `0` = channel absent) and `"L_S"` instead
  of `"ls"`. **The writer moved too** — `_key_doc` could not stay on `"ls"` — so
  v5 is what `save` now produces. v2–v4 documents back-read through the relabel
  with **no migration tool**: `test/unit/test_persist.jl` fabricates a v4
  document and asserts identical keys, identical fingerprint, `jphi` equality,
  bit-identical `predict_energy` / `predict_torque`, and term-by-term identical
  `multipole_terms` (the downstream Monte-Carlo program-array proxy).
- **`coeftable` columns are now `body`, `orbit_id`, `decors`, `L_S`, `Lf`,
  `block`, `J`.** A pure-spin row still renders as `"1,1,2"`, so the string a
  reader sees is unchanged; displacement factors would render as `u(k:l)`
  (colon inside the token, so the comma-joined column stays splittable).
  README and SPEC list the new columns (they still said `ls`); `CLAUDE.md`
  records that `MultipoleTerm.ls` keeps its name and contract — SCEMonteCarlo's
  ingest enforces it and fingerprints its values.
- `salc_groups` groups on `(body, orbit_id, decors)` — identical groups on a
  pure-spin basis. Bilinear extraction classifies on `decors` and reports any
  displacement-decorated SALC as `:unsupported` rather than dropping it.
- `MultipoleTerm` and `SALCTerm` are **unchanged** — they still carry `ls`, so
  the downstream contract (SCEMonteCarlo's program arrays) is untouched.
- Docs: the `docs/make.jl` HTML `size_threshold` is raised to 512 KiB. `api.md`
  is one page listing the whole public surface and had reached 193 KiB against
  Documenter's 200 KiB default; raised rather than split, matching upstream.

### Added — decoration labels: the shared key vocabulary (2026-08-21)

- **`src/basis/decor.jl`** — the `isbits` value labels that name what decorates a
  cluster site, ported verbatim from the spin–lattice engine so both packages
  speak one vocabulary:
  - `Channel` (`SPIN < DISP < OCC`, a `UInt8` enum whose order is load-bearing —
    it fixes the canonical factor order, the coupling order, and the persisted
    integer codes, so channels may only ever be appended). `OCC` is **reserved
    and unconstructable**.
  - `SiteFactor(channel, k, l)` with per-channel validation (`SPIN` requires
    `k = 0`, `l ≥ 1`; `DISP` requires degree `2k + l ≥ 1`), `SiteDecor` (at most
    one factor per `(site, channel)` slot), the accessors `has_spin` / `has_disp`
    / `spin_rank` / `disp_degree` / `factors` / `is_pure_spin`, the relabel helper
    `spin_decors`, and `rep_scale` — the single declared source of the
    per-channel group action (`SPIN` axial `det(R)^l`, `DISP` polar).
  - `Slot` (in `src/basis/salc.jl`): the axis → `(site, factor)` map that
    generalizes the v4 axis-`i` ↔ site-`i` identity, with the internal
    `spin_slots` producing the pure-spin identity list.
- **Pure addition.** Nothing existing changed: `SALCKey` still carries its v4
  `ls` field, every construction path is untouched, and the whole suite is
  byte-identical. The names are declared `public` (unexported — `Channel` would
  shadow `Base.Channel`) and documented under *Decoration labels* in `api.md`.
- `test/unit/test_decor.jl` covers channel order, both constructors' validation,
  the ordering law (pure-spin decors sort exactly like the sorted `ls` label),
  the slot map, and the `rep_scale` trait including the even-`Σl_spin` identity
  that makes one polar Wigner cache correct for both channels.

### Added — the pin tier: change detectors over the SALC chain (2026-08-21)

- **`test/pin/`** — a byte-level pin over five real crystals (bcc Fe, B2 FeRh,
  hcp Co, wurtzite GaN, rocksalt MnO), with its own environment, its own driver,
  and a CI job that runs it at 4 and at 1 thread. **Change detectors, not
  correctness evidence** — every expected value was produced by this package;
  `test/unit/test_normalization.jl` is the correctness half of the pair.
  - Four layers, each with its own strictness and its own recapture rule:
    **L0** structure (integers and labels — exact on every platform, and the
    precondition for the two layers indexed by position), **L0′** the sign/support
    pattern of `folded`, **L1** every `folded` entry as a raw IEEE-754 bit
    pattern, **L2** `‖X_E‖_F²` / `r2` / `coef` / held-out energies.
  - **The L0′ threshold is measured**: over 12,040 tensor entries, 7,728 are
    exactly zero and the smallest survivor is 0.239, against a prune threshold of
    1e-10 — a **9.38-decade gap**, with `eps = 5e-6` at its geometric mean. That
    gap is what makes L0′ exact across platforms.
  - **L1 is not declared portable**: the chain runs through `eigen` and BLAS, so
    it is exact on the capture platform and `rtol = 1e-12` elsewhere, with the
    platform recorded in each pin's `[meta]`. Byte-stability across *thread
    counts* IS measured — captured at `-t 4`, exact at `-t 1`.
  - **The runner names the layer that moved and the rule that applies.** There is
    no branch protection and no pin-only-commit job, so that printout is the whole
    enforcement mechanism; the rules are in `test/pin/PIN.md`.
  - L2's targets are a seeded random vector, not energies pushed through this same
    basis: measured, that round trip cancels and leaves `coef` and `r2` numerically
    identical under a uniform rescale of the whole basis.

### Fixed — the benchmark suite could not run (2026-08-21)

- **`bench/bench_solver.jl` had been dead since the export/public split**: it
  calls `solve_coefficients`, which is `public` but not exported, so the script
  died at load with `UndefVarError`. Nothing runs `bench/` in CI, so nothing
  caught it. Fixed with an explicit `using <Pkg>: solve_coefficients`.

### Added — benchmark baseline and regression rule (2026-08-21)

- **`bench/BENCH_LOG.md`** — the measured baseline for all six scripts at their
  stress defaults, plus an explicit rule for what counts as a regression: the two
  gate scripts are `bench_salcbasis` and `bench_design_matrix`; **allocation
  count is the primary gate at zero tolerance** (measured bit-identical across
  three consecutive runs, so the quantity is deterministic) and **wall time is
  secondary at 5%** (measured run-to-run spread at most 1.4%, so ~3.5x headroom).
  The other four scripts are context, with the reason each is not a gate.

### Added — absolute-normalization oracles for the SALC chain (2026-08-21)

- **`test/unit/test_normalization.jl`** — four closed-form gates over the stretch
  of the numeric chain that had no oracle: the Reynolds projector,
  `_canonical_basis`, `_canonicalize_members`, the `folded` contraction, and
  `evaluate_salc`. Everything that reached those stages before was gauge- or
  scale-invariant, so a uniform rescale of the whole basis passed the entire
  suite (measured: doubling every `folded` tensor left the orbit-sum invariance
  check AND a member-resolved invariance check green at 1e-14).
  - **O1** two sites, P1, `l₁ = l₂ = 1`: the Heisenberg channel is
    `Φ(e) = 2√3 (e₁·e₂)` — from the real-harmonic normalization
    (`test_harmonics.jl`), `‖C^{Lf}‖_F² = 2Lf+1` (`test_coupledbasis.jl`), the
    `(4π)^{N/2}` scale, and the `N!` ordering fold.
  - **O2a** one site, `l = 2`, P1: `Σ_b Φ_b(e) Φ_b(e′) = 5 P₂(e·e′)` (spherical
    addition theorem). Gauge-free — invariant under any orthogonal remixing of
    the five columns, so a legitimate `_canonical_basis` gauge change cannot
    break it.
  - **O2b** one pair member, P1, `l₁ = l₂ = l`: exactly `(2l+1)²` columns and
    `Σ_all Φ² = 4(2l+1)²` (Clebsch–Gordan completeness), measured for
    `l = 1…4`. This is the only gate in the suite that reaches `l ≥ 3`
    absolute normalization.
  - **O3** the C3v equilateral triangle: `Φ = 2√3 Σ_{i<j}(eᵢ·e_j)`, with the
    reference side built from the three geometric bonds and never from
    `s.members` — so member multiplicity is under test, which O1/O2 cannot see.
  - **Teeth measured**, not asserted: each oracle has a source mutation that
    kills it and leaves the other three green (table in the file header).

### Fixed — anisotropic bases on aliasing supercells; tie tolerance exposed (2026-08-13)

Second MnTe(0001) cross-check follow-up (the first produced the orbit-closure
gate, `7d91231`). Two defects, both confirmed **family-shared** (upstream SLCE.jl
reproduces them bit-for-bit on the same inputs — 51 SALCs, the same 14-dimensional
deficiency across the same six orbits):

- **Function-space reduction of each orbit's SALCs** (`_reduce_orbit_salcs`,
  `basis/salcbasis.jl`). `evaluate_salc` never reads a member's `shifts`, so a
  supercell that folds distinct cluster instances onto one atom set (WS-boundary
  ties kept whole, or a merged near-tie shell) aggregates their tensors; SALCs
  whose aggregate is the zero function — or linearly dependent within the orbit —
  reached the fit as a silently rank-deficient design. Measured on bulk MnTe with
  SOC (3×3×3, 108 atoms, `P6_3/mmc`): 51 SALCs emitted, rank 37, OLS returning
  `max|coef| ~1e7` with normal-looking energies and R², making every
  coefficient-level readout (`coeftable`, `bilinear_terms`, Sunny export)
  meaningless. The builder now expands each SALC into its aggregated
  `(atoms, ls, index)` monomial coefficients — exact linear algebra in the space
  evaluation actually spans — and drops zero/dependent combinations with a warning
  naming orbit, channel, and reason. Post-fix the MnTe basis is 37 SALCs,
  channel-for-channel equal to Magesty's, full rank. Closed-form gate: the CsCl
  8-fold corner-tie fixture in `test/unit/test_salc.jl` (hand invariant theory:
  the aggregate-zero `Lf = 2` drops, the survivor is `∝ e₁·e₂`).
- **`OLS` warns on a rank-deficient (or severely ill-conditioned) design**
  (`fitting/estimators.jl`). The solve is unchanged (explicit pivoted QR — the
  factorization `X \ y` uses — reused for the check), but a diagonal ratio below
  `1e-10` now warns that the coefficients are non-unique or unstable. This is a
  one-sided conditioning gate (it cannot fire while `κ₂ < 1e10`) and the backstop
  for what the per-orbit reduction cannot certify: cross-orbit dependence
  (including tied images split into separate orbits under a trivial space group),
  repeated-atom `AllImages` members, degenerate training data.
- **`SCEBasis(...; tie_tol)`** exposes the relative same-distance band
  (default unchanged at `1e-8`, hard cap `1e-2`), threaded to both the neighbor
  list and the cluster-edge admission, and carried by the TOML setup
  (`[interaction].tie_tol`, `read_setup`, `SCEBasis(path)`) so a build that
  needed a widened band is reproducible from its own file. This is the remedy the orbit-closure
  refusal names: DFT-relaxed coordinates symmetric to less than the tie band
  (the MnTe slab: ~2.4e-6 Å residual, tie splits ~2e-7 relative) split
  symmetry-partner ties differently and lose candidate closure. With
  `tie_tol = 1e-5` the 117-atom slab now builds: 147 SALCs — exactly Magesty's
  count — full rank, and every basis function invariant under all 108 operations
  (max deviation 1.8e-9 against the pre-gate ~4 meV symmetry break).

Bases persisted before this change reload verbatim, including any redundant
columns; rebuild from the setup to get the reduced basis.

### Changed — revived as the spin-only carve-out of SLCE.jl (2026-08-11)

SLCE.jl (this repository's continuation) grew into a joint spin–lattice package;
the spin-only line is carved back out here as its own package. The cut is
`698a841`, the last commit before the joint rewrite (M0). New UUID
(`4a34ec00-…`; SLCE.jl kept the original), same module name and API. The three
`docs/specs/spin-lattice-ce*` design records leave with SLCE.jl.

### Fixed — full re-audit of the remaining upstream range (2026-08-12)

The d4660a7 miss prompted a hunk-level re-audit of **every** post-carve-out
SLCE.jl commit the first audit had not opened (70 of 98; five parallel passes).
Outcome: nine commits carried spin-applicable material, all small; everything
else confirmed joint/strain/naming/feature-only.

- `2259c54` (partial, **BREAKING**) — `select_support`'s point count and explicit
  thresholds are separate keywords (`npoints` / `thresholds`): an Integer
  `thresholds = 10` used to silently mean a ten-point grid, now a TypeError.
  The batch's Greek-keyword renames stay out.
- `011e3c9` (partial) — the `SCEDataset(basis, src::AbstractDFTSource)`
  convenience constructor forwards `zero_moment_atol` instead of silently
  dropping it.
- `a1ac9af` (partial) — the suite refuses a misspelled `TEST_MODE` (it used to
  run zero tests and report success) and refuses one thread
  (`SCEFITTING_ALLOW_SINGLE_THREAD=1` overrides); CLAUDE.md's test commands say
  `-t 4`.
- `575a4e3` (partial) — CI gains the `examples` job (each example self-gates —
  the fence that caught the J/2 regression — and nothing ran them);
  `checkdocs = :public` (the unexported public surface is API too).
- `e8b8ee4` (adapted) — CI gains the `downstream` job: SCEMonteCarlo.jl's suite
  runs against every push here.
- `1495e44` (partial) — CLAUDE.md names the two independent torque-convention
  gates and repairs the dead `sce_bridge.jl` pointer.
- `fe1c9d6` (partial) — `grad_Zlm_unsafe`'s docstring writes the tangent
  projection as `∂Z − û (û·∂Z)`, matching the kernel.
- `b9230c0` (partial) — design-notes §13: what the IRLS fixed-point argument
  does and does not say (the objective is the log-sum, whose per-group factor is
  a plausible mechanism for the l044 nothing-dies observation).
- `de79b92`/`3e68fc1` (adapted) — the resolvability page corrects the
  "independently resolvable tie members" claim and documents what a WS-boundary
  tie costs on a finite cell (identically-zero SALC columns; measured upstream:
  4 of 7 on the bcc corner tie with `l ≤ 2` spin factors), with the
  bond-reversal `(−1)^{L_f}` pattern as a guide, not a rule.
- `bd632f2` (partial) — five spin-era doc corrections: the SALC projector has no
  proper/improper case distinction (the even-`Σl` screen makes `det(R)^{Σl} ≡ +1`);
  no channel is "folded into `j0`" on the Sunny path; the `:coupling` dispersion
  is exact only for uniform `S_eff`; the `clean` check's four conditions; the
  case-1 energy figure is a whole-cell RMSE; the design-notes cost-front
  percentages are labeled as measured under the old entry-count metric.
- Own (not upstream): `test/oracle/README.md` and `docs/design-notes.md` stop
  claiming the oracle kernels match "bit-for-bit" — the suite asserts
  `atol = 1e-13` / `rtol = 1e-12` (the same falsehood was corrected in SPEC.md
  by `b0593ef`; these two sites had survived in both packages).

### Fixed — the orbit builder refuses a non-group-closed candidate list (2026-08-12)

`build_clusters` silently **skipped** a symmetry image missing from the
candidate list, building the orbit short — and the SALCs projected on a short
orbit are not invariant under the space group they claim. Reachable, not
theoretical: on the MnTe(0001) 3×3 slab (relaxed coordinates, symmetric only to
~2.4e-6 Å; spglib at tol 1e-3 still reports P-6m2), the minimum-image distance
ties split at ~2e-7 relative — beyond the 1e-8 neighbor band — so
symmetry-partner pairs kept different tie images, the candidate set lost group
closure, and the basis carried **54 spurious Lf ≠ 0 SALCs** (201 vs Magesty's
147) whose fitted model broke its own space group by ~4 meV. The builder now
errors loudly, naming the missing image ([backported from SLCE.jl `d4660a7` —
missed by the original audit; found via the slab crosscheck 2026-08-12]). Also
from `d4660a7`: `build_neighbor_list(…, AllImages(); search = …)` refuses the
inert keyword instead of ignoring it. Regression: a hand-built exact C3 group
over a 1e-5-perturbed 3×3 honeycomb reproduces the refusal without Spglib.

Remedies for a refused input: symmetrize the coordinates to the reported group
(then all tie images are kept within the band and closure holds), or use
`images = AllImages()` with finite cutoffs.

### Fixed — backports from SLCE.jl (2026-07-25 … 2026-08-11)

Every post-carve-out SLCE.jl fix was audited; the ones whose defect exists in
spin-only code are ported (each code site cites its upstream SHA):

- `b97bb47` — `grad_Zlm`'s tangency is an identity in the input (`û(û·∂Z)`, the
  `/r²` the preamble always specified). Not bit-neutral: torque/gradient numbers
  move at rounding level (≤ ~4e-15 relative).
- `5b60685` + component bound — `Zlm`/`grad_Zlm` validate by the family rule
  (1e-6 band + `max|component| ≤ 1`, refusing near-pole inputs by name instead
  of a bare `DomainError` from inside the Legendre recursion); the same bound
  guards `_validate_config` at the dataset/predict doors (from `8ee6739`).
- `54457ca` (review M2/M3/M4) — the dataset door's `atol` is capped at 1e-2
  (`_check_atol`); `SCEFit` records the `refit` support and `effective_dof`/
  `gcv` refuse a refit result by name (they reconstruct the FULL design);
  `refit` refuses a `GroupAdaptiveRidge` up front; GCV counts only informative
  rows and charges the intercept only while the energy block carries weight
  (`torque_weight == 1` accounting); `select_fit`'s `:cv` score is per
  informative row; the Pareto rule is re-checked after the cold re-derivation;
  `_support_thresholds` refuses an empty group vector; `wignerD_real` asserts
  its LS residual; adaptive-estimator GCV documented as optimistic (frozen-
  weight df is a lower bound).
- `896180e` — `Ridge(lambda = 0)` routes to OLS's QR path exactly (the normal
  equations on a singular Gram returned ‖β‖ ~ 1e16 silently); `_rank_df` uses
  the `min(size)` tolerance `LinearAlgebra.rank` documents.
- `bde9ded` — `SALCBasis` and `SCEBasis` validate in inner constructors
  (key/order/injectivity contract; species table); the fingerprint is derived,
  never supplied.
- `f8a529d` — `_assemble_spacegroup` validates that the operation list is a
  group of the lattice (integrality, |det| = 1, metric orthogonality, identity,
  duplicates, closure, inverses) and `map_sym` columns must be permutations.
  The `test_nbody.jl` C3v fixture moves to a hexagonal cell where its 3-fold
  rotation is exactly integral.
- `72a9cc1` — the same-distance band `tol` travels on the `NeighborList`;
  cluster-edge admissibility reads it instead of a hard-coded constant.
- `8ee6739`/`b0593ef` — fractional wraps go through `_wrap01` (half-open
  `[0, 1)`; `mod(-1e-18, 1.0) === 1.0` silently dropped the cutoff-sphere
  shell), in `Crystal` and the Sunny primitive unfold.
- `92df67c` — a zero-SALC basis is refused at the `SCEDataset` boundary (it
  could only fit the mean energy; reachable from an ordinary spec via the
  same-atom pair minimum image cannot express).
- `23064a4` — `multipole_terms` gains `keep_zero`; the default term list is a
  function of the coefficient values (exact zeros are skipped), which an
  index-addressed consumer must not rely on.
- `563d994`/`103c05d` — the Heisenberg example/tutorial/README synthetic data
  encoded `J/2` after canonical (v4) members list each bond once; the examples
  now self-gate the recovered `J`.
- `572cbe0` (partial) — `_assemble_problem`'s energy-only branch whitens by
  `1/√n_E`, so the objective is the documented MSE at every weight.
  **Breaking**: λ for an energy-only penalized fit is `n_E` times smaller than
  before; the test λ grids moved with it.
- `a596ea3` — `group_costs` prices the sweep program (one site-program slot per
  member site of each distinct entry), not the energy program; at equal entry
  count a 3-body group was priced at 2/3 of its real sweep cost.
- `b0593ef` (partial) — the oracle README no longer claims a rev pin that never
  existed (`runtests.jl` logs the Magesty checkout + rev), the "bit-for-bit"
  header states the real tolerance contract, CI gains the Sunny/GLMNet
  extension jobs and pins `JULIA_NUM_THREADS`, and the Sunny primitive unfold
  is exercised on a two-sublattice, non-centrosymmetric fixture.

Audited and **not** ported (joint-only or not present in this code):
`9f636e7`, `966a447`, `2d256b3`, `4561b1b`, `92ac417`, `2b9d888`, `9145068`,
`cde644e`, `299597e`, `aba4df0`, `a20d73c`, `ac63d52`, the naming batches, and
all ASR / resolvability / force-channel work.

### Added — oracle fit-parity testset vs Magesty

- `test/oracle/`: an end-to-end **fit parity** testset — one shared EMBSET
  through both packages, bases with identical channel content (Magesty
  `body1 lmax = 2` + `body2 lsum = 4` ⟺ SCEFitting `lmax = [3]`,
  `lsum = [1 => 2, 2 => 4]`; SALC counts asserted equal), OLS at the
  single-block endpoints `torque_weight = 0` and `1`. Held-out energy
  predictions and the intercept must agree; held-out torques agree up to the
  **known global sign convention** (SCEFitting reports the physical/LL
  `τ = m × B`, Magesty reports `τ = −m (e × B)`; both flip target and
  predictor together, so the fitted function is identical). Closes the gap
  between the kernel-level oracle checks and "do the two packages return the
  same model".

### Changed — `SALCScratch`: allocation-free harmonic tables in the design hot loop

- The per-term `Z`/`∇Z` site tables of `evaluate_salc`/`accumulate_grad!` —
  previously fresh `Vector`s per (member, term) call — now live in a reusable
  internal `SALCScratch` workspace (dnPl buffer + pooled per-site tables,
  grown on demand). The design-matrix drivers and the predict paths thread
  one scratch per task/call; the `cache::Vector{Float64}` forms remain as a
  compatibility surface (wrapped into a scratch). Same calls in the same
  order ⇒ **bit-identical** design matrices (checked against a serialized
  pre-change reference, `nbody = 3` included). Bench (bcc-Fe 4³, 100
  configs, lmax 2): energy 543→464 ms and 767→230 MiB; torque 1176→906 ms
  and 1877→186 MiB.

### Fixed — review-pass hardening (whole-package review, 2026-07-18)

- `clebsch_gordan` now throws an `ArgumentError` for momenta beyond the `Float64`
  factorial range (`j1 + j2 + J + 1 > 170`) and on any internal overflow, instead
  of silently returning `Inf`/`NaN`. Far outside the validated small-`l` regime;
  in-range values are unchanged.
- `select_fit`: the selected row's `score`/`edof` (under `criterion = :gcv`) are
  now re-derived from the cold re-solve, like `n_alive`/`cost` already were, so
  `path.score[path.selected] == gcv(path.fit)` exactly. The selection decision
  itself still uses the warm path scores and is unchanged.
- `_bilinear_terms`: removed the dead `a > b` reverse-member transpose branch
  (canonical v4 members sort sites by `(atom, shift)`, so it was unreachable);
  a non-canonical member is now an internal error. Docstrings updated to the
  one-canonical-member-per-bond reality. No behavior change for real models.

### Added — cost-weighted group selection (`GroupAdaptiveRidge` + GCV + Pareto λ path)

- **`GroupAdaptiveRidge(column_groups, group_weights; lambda, epsilon, max_iter,
  tol)`** — in-core group extension of `AdaptiveRidge`: iterative reweighted
  ridge with one shared weight `wⱼ = v_g/(‖β_g‖² + p_g·ε)` per column group,
  approximating the weighted group-L0 penalty `λ·Σ_g v_g·1{β_g ≠ 0}` (at
  convergence a surviving group pays exactly `λ·v_g`). Degenerates exactly to
  `AdaptiveRidge` for singleton groups with unit weights; `lambda = 0` ⇒ OLS;
  `islinear` ⇒ `true`. Motivation: Monte-Carlo sweep cost is paid per
  contraction entry, and entries vanish only when a whole `(body, orbit, ls)`
  SALC group is zero — group-level elimination is the only sparsity that
  reduces MC cost.
- **Basis helpers** (public, unexported): `SCEFitting.salc_groups(basis)`
  (column → group labels by `(body, orbit_id, ls)`),
  `SCEFitting.group_costs(basis[, labels])` (per-group distinct-contraction-
  entry union over the canonical members; additive across groups), and
  `SCEFitting.cost_weights(basis; theta)` with `v_g = √p_g·(c_g/c̄)^θ` —
  `θ ∈ [0, 1]` tilts the penalty from cost-blind (`0`) to cost-proportional
  (`1`). Convenience constructor `GroupAdaptiveRidge(basis; lambda, theta, …)`.
- **`gcv(f)` / `effective_dof(f)`** (exported) — closed-form hat-matrix
  diagnostics for the linear estimators (`OLS` / `Ridge` / `AdaptiveRidge` /
  `GroupAdaptiveRidge`; adaptive members in the standard converged-weight
  sense): `effective_dof` = `tr(X(X'X+λD)⁻¹X') + 1`, `gcv` = `n·RSS/(n−df)²`
  on the assembled problem, `Inf` in the near-interpolating regime. The trace
  is an eigenproblem on the smaller Gram side — usable at `n < p`. `dof(f)`
  (the raw parametric count) is unchanged. On torque co-fits GCV is optimistic
  (correlated within-configuration rows); grouped CV is the documented ground
  truth there.
- **`select_fit(dataset, est; lambdas, criterion = :gcv|:cv, delta, …)` /
  `SelectionPath`** (exported) — warm-started descending λ path on a
  once-assembled Gram, scoring each fit by GCV or configuration-grouped K-fold
  CV (in core, deterministic seeded folds, per-fold Gram downdates), tracking
  alive groups (the `refit` scaled-magnitude support rule) and the predicted MC
  cost `Σ_{g alive} c_g`, then selecting the **cheapest λ within `(1 + delta)`
  of the minimum score** — the cost-aware generalization of `:lambda_1se`.
  Because the reweighted ridge crushes dead groups to tiny nonzero values, the
  default alive rule is a per-λ *relative* floor (1e-6 of the largest scaled
  magnitude); the effective absolute threshold at the selected λ is returned as
  `path.threshold`. The selected fit is re-solved cold (reproducible by a plain
  `fit`, its path row re-derived from that solve); follow with
  `refit(path.fit; threshold = path.threshold)` to de-bias exactly the reported
  support. `SelectionPath` is a Tables.jl source (one row per λ). Intended
  workflow: `GroupAdaptiveRidge(basis; theta)` → `select_fit` → `refit`,
  sweeping `theta` to trace the (MC cost, error) Pareto front.
- **`select_support(f; thresholds, delta, evalset, …)` / `SupportPath`**
  (exported) — the second knob: sweep the alive threshold at a fixed fit and
  trace the (cost, error) front of de-biased `refit`s (one OLS per point),
  scored on an evaluation dataset (pass a held-out slice) with the same
  Pareto rule. On real data the group-magnitude spectrum is continuous — no
  alive/dead gap for the λ path to expose — so this is where most of the
  cost–error trade is realized (validated on the production l044 model:
  38 % MC cost at a held-out torque RMSE better than the full model; 3 % at
  +19 %).
- **`cross_validate(dataset, estimator; torque_weight, nfolds, seed)` /
  `CVResult`** (exported) — generic configuration-grouped K-fold CV of any
  `fit` call: each fold refits from scratch (fold-local centering/whitening —
  nothing leaks across the split) and scores the held-out configurations in
  prediction space. Per-fold **and** pooled out-of-fold energy and torque
  RMSEs are reported independently of `torque_weight` (an energy-only fit
  still gets its torque error measured), alongside the fit's own
  `(1−w)·MSE_E + w·MSE_T` score. Deterministic seeded folds; `CVResult` is a
  Tables.jl source. Complements `select_fit(criterion = :cv)`, which whitens
  globally and only ranks a λ path.

### Changed — canonical SALC members (up to `N!`× smaller basis, persist v4)

- The SALC construction now folds its output into a **canonical, duplicate-free
  member form** (`_canonicalize_members`): the projection/transport still runs
  in the ordered-image space (where a stabilizer operation is a plain axis
  permutation — unchanged numerics), but the emitted members are normalized to
  one per physical cluster instance — sites sorted by `(atom, shift)`, shifts
  re-anchored to `shifts[1] = 0`, tensors axis-permuted and summed per site→`l`
  assignment. Previously every instance appeared once per site ordering (`N!`
  copies at `N` distinct sites), a construction-internal redundancy that leaked
  into evaluation. `Φ(e)`, gradients, fits, and SALC keys/fingerprints are
  unchanged up to floating-point regrouping (the merge pre-sums tensors that
  were previously summed after contraction); nothing is approximated.
  Measured on the production Nd₂Fe₁₄B `l044` model (nbody 3): 405,312 → 70,680
  multipole terms (5.73×), 4,392,744 → 744,636 tensor entries (5.90×) — the
  same factors apply to design-matrix/torque evaluation, model files, and
  downstream Monte-Carlo sweep cost.
- **Persist schema v4**: models/bases are saved in the canonical form (~6×
  smaller files for 3-body bases). v2/v3 documents remain readable — members
  are folded on load, so existing model files get the same speedups without a
  refit. Note the regrouping means energies of a reloaded pre-v4 model can
  differ from the previous build at the last-ulp level (bit-identical
  checkpoint resume across this version boundary is not preserved).

### Added — EMBSET reader (legacy Magesty training sets)

- `read_embset(path; n_atoms = nothing, zero_moment_atol = 1e-10)` and the
  `EmbsetFile` `AbstractDFTSource`, so a legacy Magesty training set drops
  straight into the pipeline: `SCEDataset(basis, EmbsetFile("EMBSET"))`. The
  format is code-agnostic (energy + per-atom moment and constraining-field
  vectors — exactly what `SpinDatum` stores), hence in-core rather than an
  SCETools adapter. Stricter than Magesty's reader on malformed input: the
  atom-index column must match the position in its block (Magesty silently
  ignores it), every numeric field must be finite (`NaN`/`Inf` — a failed SCF —
  is an error, not a silent training row), and every failure names the
  config/atom. More lenient on shape: block detection is token-based, so files
  without `#` separators parse (Magesty requires them). Cross-checked against
  `Magesty.read_embset` in the oracle suite.

### Added — dataset slicing / `vcat`, zero-moment guard

- `SCEDataset` now supports `length` (configuration count), configuration
  slicing `dataset[idx]` (integer vector/range, `Bool` mask, or `:`; duplicate
  indices allowed for bootstrap-style resampling), and `vcat` of datasets built
  on the same basis (checked by SALC fingerprint, so parts built from a
  persisted-and-reloaded basis concatenate; torque-bearing and energy-only
  parts do not mix). Design-matrix rows are sliced, never recomputed — the
  cheap path for train/test splits, filtering, and incremental data addition.
- `SCEDataset(basis, data::Vector{SpinDatum})` (and the `AbstractDFTSource`
  path through it) now **errors** when an atom referenced by the SALC basis
  carries a (near-)zero magnetic moment in some configuration
  (`zero_moment_atol = 1e-10` μ_B): such an atom previously entered the design
  matrix through `SpinDatum`'s placeholder `ẑ` direction and silently biased
  the fit. Unreferenced atoms (species removed with `lmax = 0`, sites outside
  every admitted cluster) are exempt. Matches Magesty's guard.

### Changed — BasisSpec truncation: per-body `lsum`, per-body × per-pair `cutoff`

**Breaking**: `BasisSpec`'s `pair_cutoff` keyword (and field) is replaced by
`cutoff`; passing `pair_cutoff` now throws with a migration hint (a scalar
`cutoff` is the exact equivalent). The `[interaction]` TOML key moves the same
way. Persisted documents bump to schema v3; **v2 files still load** (the legacy
scalar is expanded on read).

- `lsum` — a per-body-order budget on `Σl` over a cluster's sites (Magesty's
  `lsum`), as `lsum = [1 => 0, 2 => 4, 3 => 4]`, a scalar, or omitted (no cap).
  Enforced in the SALC `l`-tuple enumeration (`_enumerate_ls`), with the
  per-site ranges tightened to `lsum − (N − 1)` so oversized `lmax` values are
  never enumerated. This makes Magesty's l044/l064/l066-class models exactly
  expressible (`lmax` alone cannot cut `Σl`).
- `cutoff` — per body order **and** species pair: scalar, one pair table for
  all orders, or body-keyed (`[2 => Inf, 3 => ["Fe-*" => 6.0, "*-*" => 8.0]]`).
  Pair keys are unordered and resolve by specificity (concrete > `"A-*"` >
  `"*-*"`; equal-specificity conflicts error). The neighbor list is built at
  the element-wise max over orders and admits each pair against its own
  species-pair radius; `candidate_clusters` then checks **every** edge of an
  `N`-body cluster against that order's own radius (both image selections),
  with the `_SAME_DIST_RTOL` band applied per pair. `Inf` = no cutoff, `0`
  excludes a pair.
- `lmax` accepts label-keyed forms (`["*" => 3, "B" => 0]`) when the labels are
  supplied (`BasisSpec(labels; ...)` / `BasisSpec(crystal; ...)`); specs are
  stored resolved and dense (`lmax::Vector{Int}`, `lsum::Vector{Int}` with
  `LSUM_UNCAPPED`, `cutoff::Vector{Matrix{Float64}}`, `species_labels`).
  Unknown labels, uncovered species/pairs, and body orders outside `nbody` are
  errors (no silently-ignored sections); `display(spec)` prints the resolved
  truncation table.
- Gates: `test/unit/test_truncation.jl` — sugar/specificity resolution, the
  `lsum ≡ lmax` equivalences, per-pair neighbor-list and cluster admission vs
  an independent brute force, v3 + legacy-v2 persistence round-trips, and the
  new TOML forms. Validated on the production l044 model (Nd₂Fe₁₄B, nbody = 3,
  body2/3 `lsum = 4`) against Magesty's own fit.

### Added — allocation-free harmonics evaluation (cache-threaded variants)

- `Harmonics.Zlm_unsafe(l, m, u, cache)` and `Harmonics.grad_Zlm_unsafe(l, m, u,
  cache)`: passing a reusable `Vector{Float64}` workspace (length ≥ `l + 1`,
  contents irrelevant) hands it to LegendrePolynomials' `dnPl`, whose default
  argument otherwise allocates a fresh work vector on **every** call — previously
  the only allocation on these paths, and the dominant per-call cost of hot-loop
  consumers (SCEMonteCarlo's sweep kernels). Returned values are **bit-identical**
  to the cache-less methods, which are unchanged (gate: the NaN-poisoned-cache
  equivalence + zero-allocation testset in `test/unit/test_harmonics.jl`).
- The package's own hot consumers are wired through: `evaluate_salc` /
  `accumulate_grad!` accept an optional trailing `cache` (grown on demand inside
  the site-table builders), the design-matrix loops (`_design_energy` /
  `_design_torque`) hold one task-local workspace per threaded column, and the
  `predict_energy` / `predict_torque` SALC loops reuse one per call. Design-matrix
  cost roughly halves (`.claude/bench_log.md`): bcc Fe stress case `X_E` 1.96 →
  1.10 s and `X_T` 4.60 → 2.35 s; Nd₂Fe₁₄B (nbody 3, m 103) `X_E` 1.31 → 0.76 s
  and `X_T` 3.25 → 1.77 s. Values are bit-identical throughout (same summation
  order), so fitted coefficients do not change.

### Changed — bench suite: stress-scale defaults + Nd₂Fe₁₄B fixture

- The `bench/` scripts (which had broken on the `Interaction` → `BasisSpec` rename)
  are repaired and their defaults promoted from smoke sizes (16-atom bcc Fe, first
  shell) to the recorded seconds-scale stress baselines: 128-atom bcc Fe with
  multi-shell `cutoff = 6.0` Å (`bench_salcbasis` at `lmax = 3`, `bench_clusters` at
  `nbody = 3`, design/end-to-end at `m = 100` configs), with a `cutoff` positional
  argument added throughout. `bench_end_to_end` now also times `SCEDataset` assembly
  and the energy+torque co-fit.
- New realistic fixture: **Nd₂Fe₁₄B** (`bench/assets/nd2fe14b.toml`, `read_setup`
  schema — structure only, 68 atoms, 9 sublattice species, P4₂/mnm) with
  `bench/bench_nd2fe14b.jl` benching the full pipeline — basis build, design
  matrices, and OLS/Ridge energy+torque co-fits sized like the real training set
  (103 configs, 21 012 torque rows). Complements the high-symmetry bcc fixture with
  the few-ops / many-orbit, multi-species, non-magnetic-species regime.
- Baselines recorded in `.claude/bench_log.md` ("Stress baseline — 2026-07-14").

### Added — verification docs page (angular momentum)

- New **Verification** docs section (`docs/src/verification/angular_momentum.md`): a
  human-readable rendition of the `test/unit/test_angmom.jl` checks, recomputed at every
  docs build — a Clebsch–Gordan known-value table against closed forms (Varshalovich),
  orthonormality *and* completeness sweeps to `j ≤ 3`, the real Wigner-D functional
  identity on fresh directions (proper + improper), and the complex→real unitary closed
  form. Every section ends in an assertion, so a numerical regression fails the strict
  docs build; the unit tests remain the primary gate.
- The page closes with a worked example deriving the **tesseral CG coefficients for
  `l₁ = l₂ = 1`** from the complex ones via the package's own
  `coeff_tensor_complex` → `complex_to_real_tensor` pipeline: exact-form matrix
  displays (entries matched to `±p/√q`, unmatched throws), the `δ/√3` (Heisenberg),
  `ε/√2` (cross product, axial) and quadrupole identifications, and a proper-rotation
  equivariance gate against `wignerD_real`.

### Changed (breaking) — pre-registration API polish

- **The deprecated `Interaction` alias is removed** (it was a silent `const`, not a warning
  deprecation, and the package is unreleased — the window to drop it cheaply is now). Use
  `BasisSpec`.
- **`SCEBasis.interaction` → `SCEBasis.spec`** (`::BasisSpec`): the last remnant of the old
  name — reading `basis.interaction` and getting a `BasisSpec` back was the rename's ghost.
  Follows through `read_setup` (its named tuple now carries `spec`; the human-authored TOML
  section keeps its descriptive `[interaction]` name) and the persist schema (the basis
  document key `"interaction"` → `"spec"`, `schema_version` bumped to **2**; v1 files are
  rejected with a clean schema error).
- `BasisSpec` additionally validates `lmax` (nonempty, entries ≥ 0); `Ridge` validates
  `lambda` (finite, ≥ 0) in an inner constructor.

### Added — StatsAPI completion, `public` tiering, synthetic-predictor constructor

- `coeftable` and `islinear` now extend the **StatsAPI** generics instead of shadowing them
  (the same collision `coef`/`fit` already avoided — `using GLM`/`StatsBase` no longer
  clashes), and thin energy-block defaults `predict` (→ `predict_energy`), `residuals`
  (→ `residuals_energy`), and `r2` (→ `r2_energy`) are exported; the torque block keeps its
  explicit `*_torque` accessors.
- The public-but-unexported tier is now declared with the Julia **`public` keyword**
  (machine-checkable via `Base.ispublic` / Aqua) instead of a comment-only promise; `save`,
  `load`, `salcs`, `islinear`, `Harmonics`, and `AngularMomentum` join the declared list.
- **`SCEPredictor(basis, j0, jphi)`**: a public constructor filling `keys` from the basis,
  so synthetic models (hand-set couplings in tests, demos, downstream packages) no longer
  need the 4-argument form with `basis.salc_basis.keys`.
- `LICENSE` (MIT) and a CI workflow (`.github/workflows/CI.yml`: `TEST_MODE=all` on
  Ubuntu/macOS + a strict Documenter build) — the remaining registration blockers.

### Changed — layering and shared constants

- The bilinear / single-ion extraction (`BilinearTerms`, `_bilinear_terms`,
  `_l1_pair_matrix`, `_l2_onsite_matrix`, `_reconstruct_energy`) moves from
  `interop/sunny.jl` to **`sce/bilinear.jl`**: the public `bilinear_terms` introspection no
  longer depends on the interop layer (the layer inversion is gone); Sunny keeps only the
  primitive-cell unfold and `to_sunny`. Same code, same numbers.
- The tesseral constants `N1 = √(3/4π)`, `A2 = √(15/16π)`, `B2 = √(5/16π)` are defined once
  in `Harmonics` and consumed by the extraction (and by SCETools.jl's inverse mapping), so
  the forward/inverse conversions cannot drift. Bit-identical (same expressions moved).

### Fixed — torque-sign docstrings and test hygiene

- Three docstrings still carried the pre-Landau–Lifshitz torque sign
  (`SpinDatum.torques` in `io/dftsource.jl` — the file CLAUDE.md designates as the
  convention source —, `accumulate_grad!`, `SCEDataset`): all now state the actual
  convention, target `τ_a = m_a × B_a`, model `τ_a = −e_a × ∂E/∂e_a`. **Code was always
  correct**; the risk was a future edit "fixing" the code to match the wrong docs.
- New `test/unit/test_dftsource.jl`: the previously untested exported DFT data boundary —
  a closed-form `m·x̂ × B·ŷ = mB·ẑ` sign gate on the `SpinDatum` torque convention, the
  zero-moment placeholder branch, all validation throws, and a mock-source round trip.
- Shared test helpers (`rand_unit` / `rand_rotation` / `randcfg`) hoisted into
  `test/unit/testutils.jl` — the per-file copies overwrote each other in `Main` (warnings
  on every run; an edited copy would silently win). Bodies kept byte-identical, so all
  seeded fixtures draw the same stream.
- The diagnostics tests now check `rss`/`rmse`/`r2` against a from-scratch
  `y − predict_energy` reference instead of re-executing the accessors' own definitions;
  the oracle's Wigner-D comparison resolves the Magesty transpose convention once instead
  of accepting either per sample; dead `_connect` (superseded by `_connect_all`) deleted.

### Changed — source-tree reorganization (no behavior change)

- The 600-line `sce/model.jl` is split by responsibility: `sce/model.jl` keeps the
  pipeline **types** + constructors + config validation; `fitting/design.jl` the
  design-matrix assembly; `fitting/fit.jl` the `fit` / `refit` / `predict` logic; and
  `fitting/diagnostics.jl` the `coef` / `intercept` / residuals / R² / RMSE block.
- I/O is consolidated under `io/`: `persist.jl` and `input.jl` move there (joining
  `dftsource.jl`). The Sunny export moves to `interop/sunny.jl`, making the
  core/extension seam visible in the tree.
- The general Cartesian bilinear / single-ion extraction is renamed off the Sunny brand:
  `SunnyTerms` → `BilinearTerms`, `_sunny_supercell_terms` → `_bilinear_terms` (both
  internal/unexported), so the public `bilinear_terms` introspection no longer reads as
  Sunny-specific. `to_sunny` / `SunnyPrimitive` stay Sunny-named.

### Changed — diagnostics extend StatsAPI

- `coef`, `fit`, `nobs`, and `dof` are now **methods of the StatsAPI generics** rather than
  package-private functions, so they compose with the StatsBase / GLM ecosystem instead of
  clashing under `using` (`SCEFitting.coef === StatsAPI.coef`). Signatures are unchanged, so
  user code is unaffected. `intercept`, `refit`, `residuals_energy` / `residuals_torque`,
  and the `r2_* / rmse_* / rss_*` pairs stay package-specific (the two-observable energy +
  torque split has no single StatsAPI verb). Adds a lightweight `StatsAPI` dependency.

### Changed (breaking) — public API tiering and the `Interaction` rename

- **`Interaction` → `BasisSpec`.** The old name wrongly suggested a fitted coupling term;
  it is a basis/cluster *specification*. (An interim `const Interaction = BasisSpec` alias
  existed only within this unreleased cycle and is removed — see the entry above.)
  `BasisSpec` now validates in an **inner** constructor.
- **Export surface tiered.** The flat ~60-name export is split into the fitting workflow
  (still exported) and the *construction internals*, which are now **public but
  unexported** — reachable as `SCEFitting.build_clusters` / `SCEFitting.build_salc_basis`
  / `SCEFitting.evaluate_salc` / etc., but no longer dumped into `using SCEFitting`. The
  `SCEBasis` constructor already drives them, so typical user code is unaffected; advanced
  callers and tests qualify. Unexported: `build_neighbor_list`, `NeighborPair`,
  `NeighborList`, `interplanar_spacing`, `analyze_symmetry`, `n_ops`, `SymOp`,
  `SpaceGroup`, `build_clusters`, `ClusterMember`, `ClusterOrbit`, `ClusterSet`,
  `build_salc_basis`, `evaluate_salc`, `salcs`, `SALC`, `SALCKey`, `SALCBasis`,
  `solve_coefficients`, `AbstractTrainingDatum`. (Downstream `SCETools.jl` only uses the
  exported user API, so it is unaffected.)

### Performance — design-matrix / prediction hot path (bit-identical)

- `evaluate_salc` and `accumulate_grad!` now tabulate the per-site tesseral harmonics
  `Z_{lᵢμ}` (and `∇Z`) once per term and read them back in the multi-index loop, instead
  of recomputing them — and the expensive `dnPl` Legendre call inside — for every nonzero
  tensor entry. The per-term kernels sit behind a function barrier that specializes on the
  concrete tensor rank (`SALCTerm.folded` is stored as the rank-erased `Array{Float64}`),
  so the loop is type-stable. Same values and the same multiply/accumulate order ⇒ the
  output is **bit-for-bit identical** (verified against the Magesty oracle). Measured on a
  16-atom bcc 3-body, `lmax=2` basis × 200 configs (single thread): energy design matrix
  ~3.4× faster / 3.4× less allocation, torque design matrix ~5.1× faster / 4.8× less
  allocation. The same kernels back `predict_energy` / `predict_torque`.
- `build_salc_basis`: the coupled bases for each ordering are now built once and reused
  across every final `Lf` (the inner projection rebuilt the full chained-CG construction
  per `Lf`), and `_mfslice` returns a view instead of copying the multiplet slice. Both are
  bit-identical; modest build-time allocation reduction.

### Fixed

- `Lattice` now validates in an **inner** constructor, so the derived `reciprocal`
  (`inv(vectors)`, load-bearing for `interplanar_spacing` and the neighbor-list image
  range) can no longer be supplied independently via the auto-generated field constructor,
  and the singular-cell guard rejects **numerically degenerate** cells by a relative-volume
  threshold (`|det| > eps·‖A‖³`) rather than only exactly-zero volume.
- `build_salc_basis`'s docstring had detached and bound to the internal `_orbit_salcs`
  helper (it was inserted between the docstring and the function), leaving the exported
  `build_salc_basis` undocumented and breaking the `checkdocs = :exports` docs build. The
  docstring is back on `build_salc_basis`.

### Added — training-data boundary validation

- Spin configurations are now validated where they enter the pipeline (`SCEDataset`
  constructors, `predict_energy` / `predict_torque`): each must be `3 × n_atoms` with
  **finite, unit-norm columns** — the contract the harmonic kernels assume (they call
  `Zlm_unsafe`, which skips per-call checks). A non-normalized / NaN / wrong-shape config
  previously produced a silently biased design matrix or prediction; it now throws an
  actionable `ArgumentError` / `DimensionMismatch` naming the offending config and column.
  The `SCEDataset` constructors also check `length(configs) == length(energies)`. The unit-
  norm tolerance is the `atol` keyword (default `1e-6`). New `test/unit/test_validation.jl`.

### Added — thread-parallel SALC basis construction

- `build_salc_basis` now builds the cluster orbits in parallel (`Threads.@threads` over a
  flat orbit work list), each orbit producing its SALCs independently into a disjoint slot;
  the results are concatenated and sorted by `SALCKey` as before. The output is **byte-for-byte
  identical at any thread count** (verified: fingerprint *and* full folded-tensor content hash
  match across 1/4/8 threads) — orbits are independent and the final key-sort fixes the order.
- To make this race-free, the Wigner-D cache (`(l, g) → wignerD_real`) is now **precomputed
  serially** over the full, bounded `(l ≤ lmax, g ≤ n_ops)` grid and is **read-only** during
  the threaded loop (the previous lazy `get!` on a shared `Dict` would have raced). Speedup is
  allocation/GC-bound (the build is allocation-heavy): ~2.1× on 4 threads, ~2.4× on 8 for a
  128-atom FeRh 3-body basis (100 SALCs). Serial (1 thread) is unchanged.
- `test/unit/test_salc.jl` gains a determinism/thread-safety regression test (a rebuild must
  reproduce keys and folded tensors exactly), and `test/unit/test_nbody.jl` extends it to a
  **3-body, two-species, multi-term** orbit (`lmax_by_species = [2, 1]`, degenerate-multiset
  split) — the path where the Wigner cache bound and per-orbit `blockcount` matter most. The
  determinism / threading tests now `@warn` when run on a single thread (the threaded path is
  only exercised under `julia -t N>1`).

### Changed

- Internal cleanup following the threading refactor: dropped the now-dead `crystal` / `sg`
  parameters from `_project_and_fold` / `_transport_term` (unreferenced since the Wigner-D
  cache became a read-only `(l, g)` lookup). No behavior change.

### Added — thread-parallel design-matrix assembly and batch prediction

- The two design-matrix builders (`_design_energy`, `_design_torque`) and the vector
  `predict_energy` / `predict_torque` forms now parallelize over independent columns /
  configurations with `Threads.@threads`. Each task owns whole columns (or output slots),
  so writes are disjoint and the result is **identical at any thread count** (the per-atom
  gradient buffer `G` in the torque builder is now task-local). Speedup scales with
  `JULIA_NUM_THREADS` / `julia -t`; serial (1 thread) is unchanged. This is the cost that
  grows with dataset size — large active-learning batches and bigger supercells.
- New `test/unit/test_threading.jl` pins the threaded output to a race-free serial
  reference and to the independent scalar predict path (a data race would fail it under
  `julia -t N>1`).

### Added — `to_sunny` spin scaling routes (`:moment` / `:coupling`)

- `to_sunny` gains a `scaling` keyword so it can export a **non-half-integer** effective spin
  (e.g. the itinerant ``S_{\text{eff}} = m/(g\mu_B) \approx 1.1`` of bcc Fe), which Sunny's
  `Moment` (half-integer only) cannot carry directly:
  - `:moment` — put `S_eff` into the `Moment` and rescale exchange by `1/(SₐS_b)`; static
    energy *and* dispersion exact, but `S_eff` must be a half-integer (the previous, only
    behavior).
  - `:coupling` — keep `Moment` at a placeholder `s₀ = 1` and fold `S_eff` into the couplings
    (`J = M/(s₀√(SᵢSⱼ))`, single-ion `1/(s₀ Sᵢ)`); works for **any** positive `S_eff`. Only
    the magnon *dispersion* is physical (invariant under the overall spin scale `sᵢ→c sᵢ,
    J→J/c`); the static energy is rescaled.
  - `:auto` (default) picks `:moment` for half-integer spins, else `:coupling`. `mode` also
    gains `:auto` (`:dipole` for half-integer, `:dipole_uncorrected` otherwise). A `:coupling`
    placeholder cannot carry the quantum quadrupole, so single-ion + `:dipole` + `:coupling`
    is rejected. Modeled on the Magesty.jl Sunny export.
  - Validated: the `:moment` and `:coupling` dispersions agree to ``<10^{-7}`` for a
    half-integer spin where both apply, the `:coupling` static energy equals
    `(predict_energy − j0)/S`, and a non-half-integer `S_eff` now builds and disperses.

### Added — bcc Fe worked-example tutorial

- New tutorial [`docs/src/tutorials/case1_bcc_fe.md`](docs/src/tutorials/case1_bcc_fe.md):
  a real (non-synthetic) end-to-end fit for body-centered cubic iron — a 128-atom
  ``4\times4\times4`` supercell, isotropic two-body basis (``Im\bar3m``, 13 SALCs), fit from
  noncollinear spin-DFT energies and torques, with in-sample validation (parity plots),
  the isotropic ``J_{ij}`` read back out of [`bilinear_terms`](@ref) against pair distance,
  and a **live** [`to_sunny`](@ref) magnon dispersion (``\Gamma\text{–}H\text{–}N\text{–}\Gamma\text{–}P\text{–}H``
  and a ``\Gamma\text{–}N`` comparison to neutron data) using the new `:coupling` route at
  ``S_{\text{eff}} = 1.1``. The fit uses `torque_weight = 1.0` (a torque fit, constraining the
  energy gradient that sets the dispersion), reproducing the reference Magesty.jl tutorial's
  ``J_{ij}`` and magnon spectrum. The whole pipeline runs **live** in the docs (the ~10× `SCEBasis`
  speedup below makes the 128-atom build a few seconds). Ships the reference data — `POSCAR`,
  a 50-configuration `EMBSET`, and the experimental `febcc_spinwave.csv` — under
  `docs/src/tutorials/case1_inputs/`; the upstream VASP I/O and mean-field sampling that
  produced it live in the companion `SCETools.jl`. Adds `Sunny` to the docs environment.

### Performance — `SCEBasis` build (clusters + SALC)

- The symmetry-orbit and SALC construction were the dominant cost of building a basis on
  a large supercell (minutes for the 128-atom bcc Fe tutorial cell). Four numerics-
  preserving changes cut it to ~1.5 s (≈10×), and removed the GiB-scale transient
  allocations:
  - `build_clusters` reduces orbits by **growing each orbit from one representative**
    (`O(n_orbits · n_ops)`) instead of computing a canonical key per candidate
    (`O(n_candidates · n_ops)`); the canonical-key inner loop is now allocation-free
    (`Val(N)` static site tuples).
  - `build_salc_basis` connects all orbit members to the representative in **one
    `O(n_ops)` sweep** (was an `O(n_ops)` scan per member), and **memoizes the real
    Wigner-D matrices** `D^l(R_g)` by `(l, g)` across the build.
- Output is **byte-identical** (same orbits, same SALC keys/coefficients): verified by the
  full unit suite and the gauge-invariant Magesty oracle. Benchmarks and before/after
  numbers live in [`bench/`](bench/) and `.claude/bench_log.md`.

### Changed — public API naming consistency (BREAKING)

- A naming-and-usability pass renamed several exported symbols for consistency. All are
  **clean renames** (no deprecation aliases); update call sites accordingly:
  - **Count accessors** unified to an `n_*` prefix: `num_atoms` → `n_atoms`,
    `nsalc` → `n_salcs`. (`nobs` / `dof` are StatsAPI names and are unchanged.)
  - **SALC access** de-stuttered: the `SCEBasis` field `.salcs` (a `SALCBasis`) is renamed
    to `.salc_basis`, and a new exported accessor `salcs(basis)::Vector{SALC}` returns the
    basis functions directly — write `salcs(basis)[k]` instead of the old
    `basis.salcs.salcs[k]`.
  - **Predictor type** `SCEModel` → **`SCEPredictor`** (the lightweight, persistable
    predictor): the type, its constructor `SCEPredictor(fit)`, and `load(SCEPredictor, …)`.
  - **Setup reader** `read_input` → `read_setup` (returns the parsed crystal + interaction +
    symmetry setup; pairs with `read_configs` for training data).
  - **SALC kernel** `evaluate` → `evaluate_salc` (a less generic, collision-resistant name
    for the exported invariant evaluator).
- The persisted-model/basis **TOML schema is unaffected** (the on-disk tags `sce-basis` /
  `sce-model` and all doc keys are unchanged), so model and basis files written by an
  earlier build still load.
- Docs additions: the `SCEFit` (heavyweight, data-bearing) vs `SCEPredictor` (lightweight,
  persistable) roles are now contrasted in their docstrings, the README, and Getting
  Started; the `coeftable` columns (`body` / `orbit_id` / `ls` / `Lf` / `block` / `J`) gained
  a legend in the I/O guide.

### Added — lattice figures in the tutorials

- The Heisenberg-chain and kagome three-body tutorials now open with a generated lattice
  figure (CairoMakie) so the system is shown, not just described. Each figure draws the
  **unit (calculation) cell** and the periodic connectivity from the **actual** computed
  geometry — sites from `cartesian_positions`, bonds from `build_neighbor_list`, and the
  highlighted kagome triangle read back from `build_clusters` (its sites and periodic
  `shifts`). The chain figure shows the four sites closing into a ring (the dashed
  ``1\text{–}4`` bond across the cell boundary); the kagome figure shows the cell's three
  sites and a 3-body cluster that borrows two of its sites as periodic images from the
  neighbouring cell. The pictures cannot drift from what the basis is built on. `CairoMakie`
  is a `docs/` dependency only; the package itself gains no plotting dep.

### Changed — package renamed `MagestyRebuild` → `SCEFitting`

- The package, its module, and the repository directory were renamed from
  `MagestyRebuild` (`Magesty_rebuild.jl`) to **`SCEFitting`** (`SCEFitting.jl`), unifying the
  naming with the companion `SCETools.jl` under a shared `SCE*` family. The UUID is unchanged,
  so the package identity is preserved. The package extensions are now
  `SCEFittingGLMNetExt`, `SCEFittingSpglibExt`, and `SCEFittingSunnyExt`. The persisted-model
  TOML schema tags changed from `magesty-rebuild/sce-{basis,model}` to
  `scefitting/sce-{basis,model}` (`schema_version` is still `1`); any model TOML written by an
  earlier build must be re-saved. Downstream code updates `using MagestyRebuild` to
  `using SCEFitting`. The legacy `Magesty.jl` package (the design-reference original) is
  unaffected and keeps its name.

### Changed — VASP I/O moved to `SCETools.jl`

- The concrete VASP adapter (`SCEFitting.VASP`: `read_poscar`, `write_poscar`, `Oszicar`)
  has been **moved to the `SCETools.jl` package** (`SCETools.VASP`), joining the INCAR writer
  so all VASP I/O lives in one place. The core now keeps only the **code-agnostic DFT-data
  seam** — `AbstractDFTSource`, `AbstractTrainingDatum`, `SpinDatum`, `read_configs`, and
  `SCEDataset(basis, src)` — so the SCE pipeline stays code-agnostic. To read VASP training
  data, `using SCETools` and `SCETools.VASP.read_poscar` / `Oszicar` (the `SCEDataset(basis,
  src)` seam is unchanged). The `test/unit/test_vaspio.jl` unit tests, the
  `examples/vasp_dft_source.jl` example, and the VASP-vs-Magesty oracle cross-check moved to
  SCETools with it.

### Changed — sampling extracted into `SCETools.jl`; fitted-model introspection added

- The mean-field spin-configuration **sampler** (the P0–P4 work documented below) has been
  **moved out of this package** into the new auxiliary package `SCETools.jl`, which depends
  on `SCEFitting`. This package is now focused on building and fitting SCE models;
  generating spin configurations (and, later, active learning) lives in `SCETools.jl`. The
  removed exports are `AbstractSampler`, `MFASampler`, `MFASample`, `ExchangeModel`,
  `MultipoleField`, `sample`, `mfa_temperature_scale`, `mfa_sublattice_m`,
  `thermal_averaged_m`, and `tau_from_magnetization` (now exported by `SCETools`).
- **Added a code-neutral fitted-model introspection surface** so a downstream consumer reads
  the fitted Hamiltonian without reaching into the SALC-basis internals:
  - `multipole_terms(model)` returns a flat `Vector{MultipoleTerm}` — one record per cluster
    member / `l`-ordering of every SALC with a nonzero coefficient, carrying the raw `jϕ`
    coefficient (the per-N scale `(4π)^(body/2)` left for the consumer), the `body`, member
    `atoms` / `shifts`, per-site `ls`, and the `folded` tensor;
  - `bilinear_terms(model)` returns `(; pairs, onsites, skipped)` — the bilinear (`ls=[1,1]`)
    and single-ion (`ls=[2]`) channels as Cartesian `3×3` matrices, reusing the validated
    Sunny export conversion;
  - `num_atoms(model::SCEModel)` is a new method of the exported `num_atoms`.
  The gate is energy reconstruction (`test/unit/test_introspect.jl`): summing the per-term
  tesseral contraction reproduces `predict_energy − j0`.
- The tesseral spherical-harmonic submodule `SCEFitting.Harmonics` (`Zlm`, `lm_index`) is
  documented as a stable surface for downstream packages.

### Added — mean-field spin-configuration sampling: P4 (full multipole / many-body)

> The P0–P4 sampler entries below are retained as history; the code now lives in
> `SCETools.jl` (see the *Changed* entry above).

- **`MultipoleField` and `MFASampler(model::SCEModel; reference)`**: the full multipole
  mean-field sampler over **all** SCE clusters and harmonic orders (higher-order /
  many-body). The mean-field decoupling of any cluster term factorizes
  (`⟨∏ Z⟩ → ∏⟨Z⟩`), giving the generalized molecular field
  `h_a^{lm} = Σ_φ jφ·(4π)^(N/2)·folded · ∏_{b≠a} ⟨Z_{l_b}^{m_b}(e_b)⟩` — built by contracting
  each folded coefficient tensor against the *other* sites' multipole averages (the
  `accumulate_grad!` leave-one-out structure with the site-`a` harmonic left symbolic). The
  order parameters generalize from the magnetization to the **full per-atom multipole
  averages `⟨Z_lm⟩_a`** (`l ≤ lmax`), iterated to self-consistency by sphere quadrature; the
  `l=1` (bilinear) Perron eigenvalue sets `T_MF = ρ/3`, and the single-site Bingham /
  higher-multipole distribution is drawn with the Metropolis engine.
- `MFASampler(model::SCEModel; reference)` keeps every channel (bilinear, single-ion, and
  higher-order); `MFASampler(ExchangeModel(model); reference)` remains the bilinear-only
  truncation. Validated by the exact reduction to the single-global Langevin curve for a
  pure-bilinear model, by coupling-scale invariance, and by matching the single-site
  potential to the conditional mean SCE energy `⟨E | e_a⟩` of a biquadratic model (the
  many-body factorization check).

### Added — mean-field spin-configuration sampling: P3 (tensorial + single-ion)

- **Tensorial `ExchangeModel`**: `ExchangeModel` now carries the full bilinear tensor
  `bilinear[a,b] = S_ab` (Heisenberg + DM + anisotropic exchange) and single-ion anisotropy
  `onsite[a] = A_a` (the `ls=[2]` channel), with an `isotropic` flag for the P2 fast path.
  `ExchangeModel(model::SCEModel)` now extracts **all** of these (only the higher-order /
  higher-`l` SALCs are dropped); `ExchangeModel(Jiso; onsite)` and
  `ExchangeModel(bilinear; onsite)` are the raw constructors.
- **`MFASampler(exch; reference)` — tensorial path**: with DMI / anisotropic exchange or
  single-ion anisotropy the single-site potential `V_a(e) = β(e·g_a + e' A_a e)` (molecular
  field `g_a = Σ_b S_ab m_b ê_b`, `β = 3/(ρτ)`) gains an `l=2` Bingham factor that the
  closed-form vMF cannot represent, so the magnetizations are solved as `m_a = ⟨e·ê_a⟩` by
  sphere quadrature and the configurations are drawn with the **Metropolis** engine
  (per-atom chains, proposal scaled to the peak sharpness). The longitudinal molecular-field
  matrix `A[a,b] = −ê_a' S_ab ê_b` (so `T_MF = ρ/3`) generalizes the P2 isotropic case,
  which still takes the closed-form vMF path. **Noncollinear references** are supported
  (rigid cone axis, decision D2; a warning flags a non-stationary reference, e.g. DMI
  canting a collinear state). Because the single-ion enters with the same `β` as the
  exchange, the anisotropy persists above `T_MF` (an anisotropy-weighted paramagnet).
- An easy-axis single-ion sharpens the cone (and keeps order above the exchange `T_MF`); an
  easy-plane single-ion gives a girdle; the sampled `⟨Z_2m⟩` match the quadrature
  self-consistency. The shared sphere quadrature now caps its auto-sized node count (the
  large molecular fields of the `τ → 0` limit would otherwise blow it up).

### Added — mean-field spin-configuration sampling: P2 (multi-sublattice, isotropic)

- **`ExchangeModel`** (`docs/specs/mfa-sampling.md`): the neutral carrier of the isotropic
  bilinear exchange the mean-field sampler needs. `ExchangeModel(model::SCEModel)` extracts
  the symmetric `Jiso[a,b] = Σ_R J_iso(a,b,R)` from a fitted SCE by reusing the Sunny
  bilinear extraction and keeping the Heisenberg part `tr(M)/3` of each bond; the DMI /
  anisotropic, single-ion, and higher-order channels are dropped (P3/P4) and reported via
  `@warn`. `ExchangeModel(Jiso)` takes a raw symmetric matrix (an external-`Jij` reader is
  a P5 target).
- **`MFASampler(exch::ExchangeModel; reference)`**: the multi-sublattice sampler. The
  per-atom magnetizations `m_a(τ)` are solved from the coupled mean-field self-consistency
  `m_a = L(3(Ā m)_a/τ)`, where the molecular-field matrix `A[a,b] = −Jiso[a,b](ê_a·ê_b)`
  folds the reference directions in (so ferro / antiferro / ferri order all become
  ferromagnetic in the magnitude variables) and `Ā = A/ρ` is normalized by the Perron
  eigenvalue `ρ`, with `T_MF = ρ/3`. Distinct sublattices disorder at distinct rates under
  a single `T_MF`; because `Ā` is scale-free, only the coupling *ratios* matter — `m_a(τ)`
  is invariant under an overall coupling rescaling (decision D4). Each spin is then drawn
  from `vMF(ê_a, κ_a)` with the self-consistent per-atom concentration `κ_a = 3(Ā m)_a/τ`.
- **`mfa_sublattice_m(sampler, τ)`** returns the per-atom `m_a(τ)`. The coupled solve uses
  depth-1 Anderson acceleration (plain iteration suffers critical slowing as `τ → 1⁻`), and
  the construction verifies the reference is a stationary, sign-definite ordered state
  (warns on a non-stationary noncollinear or frustrated reference — the rigid-axis MFA is
  exact only for collinear references, decision D2). `MFASample.m` is now a per-atom vector
  per config.

### Added — mean-field spin-configuration sampling: P1 (single global, isotropic)

- **`MFASampler(reference)`** and the **`sample`** verb (`docs/specs/mfa-sampling.md`):
  the single global, isotropic mean-field sampler, building on the P0 single-site engine.
  Each spin is drawn from a von Mises–Fisher distribution `vMF(ê_a, κ)` about its reference
  direction, with one global concentration `κ = 3m/τ` fixed by the classical-Heisenberg
  mean-field self-consistency `m = L(3m/τ)` (`L` = Langevin function) in the reduced
  temperature `τ = T/T_MF`. Numerically equivalent to Magesty's `MfaSampling`, but with an
  explicit seeded `rng::AbstractRNG` and no Roots.jl dependency (a self-written bisection
  solves the monotone self-consistency).
- **`sample(sampler, n; tau | m, …)`** draws `n` configs at one control value;
  **`sample(sampler; tau | m, nsamples, …)`** sweeps a collection. Both return an
  **`MFASample`** (decision D1): `.configs::Vector{Matrix{Float64}}` plus parallel labels
  `.tau` and `.m`, iterable/indexable as its configs. Keywords `fixed` / `uniform` /
  `randomize` carry over from Magesty's `mfa_sweep`.
- **`AbstractSampler`** is the dispatch seam for the later model-backed samplers (P2+);
  **`thermal_averaged_m`** / **`tau_from_magnetization`** expose the self-consistency and
  its inverse; **`mfa_temperature_scale`** returns `T_MF` (decision D4 — `1.0`, reduced
  units, for the coupling-free global sampler).

### Added — `refit` and the regression-diagnostic accessors

- **`refit(f, estimator = OLS(); threshold = 0.0)`**: re-solve on the **support** of an
  existing fit — the de-biasing step that follows a sparse fit. A column survives when its
  scaled-magnitude contribution `|coef(f)[j]|·‖X[:, j]‖` exceeds `threshold` (`0` keeps the
  nonzero support exactly); coefficients off the support are zeroed and `j0` is recovered
  analytically. An empty support returns an all-zero `jϕ` (with a warning) and `j0 =
  mean(y_E)`. A `PrecomputedPilot` (or an `AdaptiveLasso` carrying one) is rejected — its
  fixed coefficient vector has the original column count, not the refit support length.
- **Diagnostic accessors**: `dof` (`length(coef(f)) + 1`), `rss_energy` / `rss_torque`
  (residual sums of squares), and `residuals_energy` / `residuals_torque` (the raw residual
  vectors). The energy and torque blocks are reported separately throughout, matching the
  existing `r2_*` / `rmse_*` split; `r2_energy` / `rmse_energy` / `r2_torque` / `rmse_torque`
  are refactored to build on the new `rss_*` so each metric has a single source of truth.
- Internal: the `(X, y)` centering / whitening / `groups` assembly is factored out of `fit`
  into `_assemble_problem`, shared by `fit` and `refit` so the two build identical designs
  (the oracle confirms `fit` is numerically unchanged).

### Added — adaptive / L0-approximating estimators

- **`AdaptiveRidge`** (in-tree, no extension): iterative reweighted ridge
  (Frommlet & Nuel 2016) that approximates an L0 penalty. It refits the analytic weighted
  ridge `(X'X + λ·Diagonal(w)) \ X'y` with `wⱼ = 1/(βⱼ² + ε)` until the relative ∞-norm
  change drops below `tol` (or `max_iter` steps), so large coefficients keep a light
  penalty and small ones are driven toward zero. `lambda = 0` reduces to `OLS`;
  `islinear ⇒ true` (a linear smoother in the converged-weight sense).
- **`AdaptiveLasso`** (type in core, GLMNet solve in `ext/SCEFittingGLMNetExt`): the
  one-shot Adaptive Lasso (Zou 2006). A `pilot` estimator (default `OLS`, any estimator
  allowed) supplies `β̂`, then a weighted Lasso is solved with per-column penalty factor
  `wⱼ = 1/max(|β̂ⱼ|, ε)^γ`. `gamma = 0` reduces exactly to a plain `Lasso`. It shares
  `ElasticNet`'s `lambda` / `standardize` / grouped-CV behavior (`lambda = nothing` selects
  λ by configuration-grouped CV with the adaptive weights held fixed).
- **`PrecomputedPilot`** (in-tree adapter): returns a fixed coefficient vector from
  `solve_coefficients` (length-checked against `size(X, 2)`), so an `AdaptiveLasso` can
  reuse a prior fit's `coef(f)` as its pilot instead of re-running a pilot regression.
- All three honor the column-centered `(X, y)` / analytic-`j0` contract. The GLMNet
  plumbing is shared between `ElasticNet` and `AdaptiveLasso` (an optional `penalty_factor`).
  Validated in `test/glmnet/` (support recovery, `γ = 0` ≡ plain Lasso, pilot
  pluggability, grouped-CV co-fit) and `test/unit/test_fit.jl` (core `AdaptiveRidge` /
  `PrecomputedPilot` solves, construction / validation, deferred-backend error).

### Changed — torque sign convention → Landau–Lifshitz / physical torque

- The per-atom torque is now `τ_a = −e_a × ∂E/∂e_a` (the physical / Landau–Lifshitz
  torque `m_a × B_eff,a`), matching the published *General spin models* paper
  (Phys. Rev. Research **8**, 023300 (2026)). Previously it was the energy-rotation-gradient
  `+e_a × ∂E/∂e_a` (the methods-paper Eq. 15 convention), the opposite sign.
- The change flips **both** sides together, so the fit is unchanged: `predict_torque` /
  `_design_torque` now compute `∂Φ/∂e × e = −e × ∂Φ/∂e`, and the DFT training target from a
  constrained OSZICAR is `τ_a = m_a × B_a` (was `−m_a × B_a`). Because the torque design
  matrix `X_T` and the target `y_T` both negate, the co-fit objective `‖X_T J − y_T‖²` is
  invariant — **`j0` and every coefficient `J` are identical**; only the *reported* torque
  sign changes. The energy fit, energies, and `predict_energy` are untouched.
- `predict_torque(model, config)` and the torque returned by `read_configs` / `SpinDatum`
  flip sign for downstream consumers. The finite-difference self-consistency gate
  (`test_torque.jl`, `test_nbody.jl`) and the oracle's analytic Heisenberg torque now check
  the `−e × ∇E` convention. Docs, examples, and the design-matrix convention notes updated.

### Added — Documenter.jl documentation site

- A full browser-viewable documentation site under `docs/` (Documenter.jl): home, a
  getting-started page, a four-part guide (building the basis, data and fitting,
  persistence and I/O, Sunny export), two narrated tutorials (Heisenberg chain, kagome
  three-body) with executed `@example` blocks, a three-part theory section (the SCE
  formalism, periodic resolvability, the architecture), and a complete API reference.
  Build with `make -C docs serve` (or `build` / `open`). Local-only for now
  (`remotes = nothing`, no `deploydocs` until a remote exists).
- **`docs` fixed several exported public APIs that silently had no runtime docstring** —
  a comment or a sibling definition sat between the docstring and the documented binding,
  so Julia never attached it. `coeftable`, `load`, `intercept`, `nobs`, `nsalc`,
  `rmse_energy`, `rmse_torque`, and the `AbstractDFTSource` / `AbstractTrainingDatum`
  abstract types now carry their own docstrings; every exported binding is documented
  (the API reference builds with `checkdocs = :exports`). No behavior change.

### Tested — N-body Wigner–Seitz cluster counting (`N ≥ 3`)

- Added `test/unit/test_ws_nbody.jl`, pinning the count of 3- and 4-body clusters on
  the Wigner–Seitz boundary against an **independent brute-force enumeration** of all
  compact clusters. The candidate set (for `N = 2, 3, 4`, including `pair_cutoff = Inf`)
  and the symmetry-orbit partition are checked to match exactly on cells deliberately
  seeded with face / edge / corner ties (cubic face-atoms, fcc, skewed hexagonal), and
  every emitted member is re-verified to have all `C(N,2)` edges at the minimum image.
- Documents and guards the **compact-cluster / third-edge** criterion: a cluster is
  admitted only when *all* its pairwise edges sit at their atom-pair minimum image
  simultaneously — individually minimum-image-resolvable pairs are not enough (the images
  minimizing `i–j` and `i–k` may push `j–k` onto a longer image). The regression includes
  an equal-spaced 1-D ring where every pair is minimum-image yet no compact triangle
  exists, and an over-/under-merge guard on the orbit reduction under a cubic point group.
  No behavior change — the enumeration was already correct; this makes it guaranteed.

### Added — GLMNet estimators (Lasso / elastic-net)

- **Estimator types in core** (`fitting/estimators.jl`): `ElasticNet(; alpha, lambda,
  standardize, nfolds, select, seed, nlambda)` and the `Lasso(; …)` convenience
  (`alpha = 1`). Like the Spglib/Sunny seam, the types live in the core package — named,
  validated, and dispatched on without the heavy dependency — while the actual solve
  lights up only under `using GLMNet`.
- **`solve_coefficients(::ElasticNet, X, y; groups)`** in the new
  `SCEFittingGLMNetExt` extension. GLMNet minimizes
  `(1/2n)·‖y − Xβ‖² + λ·[(1−α)/2·‖β‖₂² + α·‖β‖₁]` on the column-centered `(X, y)` the
  fit hands it, with `intercept = false` (so `j0` is still recovered analytically) and
  column `standardize` (the penalty acts per-column at `λ·std`, returning β on the
  original scale). `lambda = nothing` selects the penalty by K-fold cross-validation over
  GLMNet's automatic λ path — `select = :lambda_min` (lowest CV error) or `:lambda_1se`
  (sparsest within one SE); a numeric `lambda` fits at exactly that penalty.
- **Configuration-grouped, reproducible CV.** `fit` now passes per-row `groups` labels
  to `solve_coefficients`; for an energy+torque co-fit a configuration's energy row and
  all its torque-component rows share a label, so CV folds never split one configuration
  (which would leak within-configuration structure and bias λ selection). Folds are
  assigned deterministically by a seeded `hash` ranking — reproducible without taking on
  a `Random` dependency in the extension.
- The `AbstractEstimator` contract gains an optional `groups` keyword on
  `solve_coefficients`; `OLS`/`Ridge` accept and ignore it.
- Validated in a separate `test/glmnet/` environment (heavy Fortran-backed dependency,
  mirroring `test/sunny/` and `test/oracle/`): tiny-λ ≈ OLS, the analytic-`j0` centering
  invariant, CV support recovery and sparsity, `:lambda_1se` shrinkage, seeded
  reproducibility, ElasticNet ≠ Lasso, and an energy+torque grouped-CV co-fit. Core-only
  construction/validation and the deferred-backend error are covered in the main suite
  (`test/unit/test_fit.jl`).

### Added — Sunny.jl export (supercell + primitive-cell routes)

- **Conversion core** (`sce/sunny.jl`): turns a fitted `SCEModel` into the
  Sunny-representable channels — `ls=[1,1]` 2-body → a 3×3 exchange matrix
  (`M = (3/4π)·folded`, carrying Heisenberg / Dzyaloshinskii–Moriya / symmetric-Γ
  in one matrix) and `ls=[2]` 1-body → a traceless-symmetric single-ion tensor.
  `ls=[0…]` channels fold into `j0`; every other SALC (3-body+, higher `l`) is
  **skipped and reported**, since Sunny cannot represent it. The directed cluster
  members fold into one matrix per undirected supercell bond
  (`_sunny_supercell_terms`), and the whole conversion is gated **without Sunny** by
  reconstructing the energy (`_reconstruct_energy ≈ predict_energy − j0`).
- **`to_sunny(model; spins, g, mode, placement)`** in the `SCEFittingSunnyExt`
  extension (loaded by `using Sunny`): builds a real `Sunny.System` on the training
  supercell (P1, inhomogeneous), placing `set_exchange_at!` per bond (rescaled
  `J = M/(SₐS_b)`) and `set_onsite_coupling_at!` per atom. The system's classical
  energy reproduces `predict_energy − j0` exactly; exchange is independent of the
  spin length and Sunny mode, the single-ion term carries the
  classical/`:dipole`-quantum rescaling. Skipped channels are surfaced via `@warn`.
- **Primitive-cell unfold** (`placement = :primitive`/`:auto`, `_sunny_primitive`):
  recovers the chemical primitive cell from the space group's pure translations, groups
  supercell atoms into sublattices, and folds the supercell bonds onto one Sunny bond
  per **primitive** bond `(i, j, n)` (no multiplicity — Sunny's periodic replication
  restores it), for *unfolded* spin-wave dispersion. A `clean` flag detects when the
  model does not live on the primitive cell (interaction range reaching the supercell
  boundary), and the export falls back to the exact supercell route.
- Conversion math is a **core dependency-free** layer; only the `Sunny.System`
  assembly lives in the extension. Validated by `test/unit/test_sunny.jl` (Sunny-free:
  the `Z₁`/`Z₂` contractions, classification, energy reconstruction, the primitive fold
  and its clean detection, skip reporting) and the separate `test/sunny/` environment
  (the real `Sunny.System` energy vs the SCE energy across modes/spins for both routes
  — the primitive system reshaped back to the supercell reproduces the SCE energy,
  confirming the unfolded bonds and offsets — plus the skip warning and per-species
  spins).

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
- **VASP adapter** (`io/vasp.jl`, `module SCEFitting.VASP`): `read_poscar` /
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

- **Persistence** (`sce/persist.jl`): `SCEFitting.save(path, x)` and
  `SCEFitting.load(SCEBasis | SCEModel, path)` serialize a self-contained,
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
  `SpglibBackend` whose method lives in `ext/SCEFittingSpglibExt`.
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
