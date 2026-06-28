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
`τ_a = e_a × ∂E/∂e_a` (the analytic derivative of the same energy surface). An
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
  reused-atom clusters dropped); `pair_cutoff = Inf` is the whole WS cell. `AllImages`
  (every image, `R`-distinguished) is **only** for the future generalized-Bloch /
  spin-spiral path where `e^{iq·R}` resolves the images. For `cutoff < min_d dᵢ / 2` the
  two coincide. Do **not** "fix" a `> L/2` cutoff by folding aliases into a shorter shell
  (double-counts) — that regime is simply unresolvable from one supercell.
- **Energy units**: `Jφ` carry the DFT input unit (eV); `j0` is separate.
- **Torque**: `τ_a = e_a × ∂E/∂e_a`, design-matrix entry `(4π)^(N/2)·(e_a × ∂Φ/∂e_a)`
  (same scale and `μ`-mapping as the energy kernel). The torque design matrix `X_T`
  has no `j0` column. Its rows are flattened config-major, then atom-major, then
  `xyz`. The co-fit whitens the (centered) energy block by `√((1−w)/n_E)` and the
  torque block by `√(w/n_T)`; `j0` stays an energy-only quantity.

## Coupled ("linked") code sites — change one, check all

- `basis/Harmonics.jl` (`Zlm`, `grad_Zlm`) ↔ the on-sphere central-difference and
  closed-form agreement tests (`test/unit/test_harmonics.jl`). Normalization / sign
  drift silently biases `X`.
- **Energy kernel `evaluate` ↔ gradient kernel `accumulate_grad!`** (`basis/salc.jl`):
  identical `μ = idx[i] − ls[i] − 1` mapping, `ls`, `folded`, and `(4π)^(N/2)` scale.
  The gate is the finite-difference self-consistency `predict_torque ≈ e × ∇E_FD`
  (`test/unit/test_torque.jl`): the torque must be the exact derivative of the energy
  surface. Change one kernel, re-check the other.
- **Image selection ↔ neighbor list ↔ cluster edges** (`geometry/neighborlist.jl`,
  `clusters/enumerate.jl`, `sce/model.jl`): `SCEBasis` threads one `images` value to
  **both** `build_neighbor_list` and `candidate_clusters`/`build_clusters`; they must
  agree. `MinimumImage` keeps minimum-image pairs (no `i==j`) and admits a clique edge
  only at its atom-pair minimum-image distance with all atoms distinct; `AllImages` keeps
  every in-cutoff image and admits edges within the radial cutoff. The tie/cutoff
  tolerance is relative (`_SAME_DIST_RTOL`) on both sides so a degenerate WS-boundary
  shell is never split. The minimum-image search box is adaptive — change it and re-check
  the skewed-cell test. `images` is **not** persisted (the full SALC basis is stored and
  reloaded verbatim), so only `read_input`/`SCEBasis` carry it.
- **SALC construction ↔ the ground-truth invariance test** (`test/unit/test_salc.jl`,
  `test/unit/test_nbody.jl`, `test/oracle/runtests.jl`): every SALC must satisfy
  `Φ(g·e) = Φ(e)` (non-collinear spins, all `Lf`, **all body orders**) and
  `Φ(−e) = Φ(e)`. The combined-space projection action (`_project_and_fold`) and the
  member transport (`_transport_term`) must use the *same* rotation direction and
  axis-relabel convention (`invperm(perm)`); the eigenvalue-exactly-0/1 idempotency
  assertion and the invariance test are the gates. At `N ≥ 3` also confirm SALCs are
  linearly independent (design-matrix rank = #SALC).
- **Design-matrix columns are identified by `SALCKey`** (`SALCBasis.keys`, sorted),
  not by construction order. The key must stay **injective**: `block` runs across all
  canonical `l`-orderings that share one sorted `ls` label (a proper-subgroup site
  stabilizer splits a degenerate multiset into several ordering orbits — see
  `test/unit/test_nbody.jl`). `SCEModel` re-pairs `jphi` to a basis **by key** on any
  reload (positionally paired only within a session); `fingerprint = hash(sorted keys)`
  guards against cross-basis confusion. This reload is realized by persistence
  (`MagestyRebuild.load(SCEModel, …)` rebuilds `jphi` in basis-key order from the
  per-`SALCKey` coefficients).
- **Persistence schema ↔ the serialized structs** (`sce/persist.jl`): `_to_doc` /
  `_from_doc` mirror the fields of `Crystal` / `Lattice` / `SpaceGroup` / `Interaction`
  / `SALCKey` / `SALCTerm` / `SALCMember` / `SALC` / `SCEBasis` / `SCEModel`. Add or
  rename a field on any of these and both halves (and `test/unit/test_persist.jl`'s
  round-trip) must follow; the space group is rebuilt via `_assemble_spacegroup` from
  the stored fractional ops, and the `UInt64` fingerprint is stored as a string and
  **recomputed** on load (never trusted — `hash` is Julia-version dependent). The TOML
  input reader (`sce/input.jl`) mirrors only the *setup* structs (crystal + interaction
  + symmetry), not the SALCs.
- **`coeftable` columns ↔ `SALCKey` fields** (`sce/coeftable.jl`): each result row is
  read straight off a `SALCKey` (`body` / `orbit_id` / `ls`→comma string / `Lf` /
  `block`) plus `jphi`; the `J` column pairs with `basis.salcs.keys` **positionally**
  (same order as the design matrix). Add or rename a `SALCKey` field → update the row
  builder, the `Tables.Schema`, and `test/unit/test_coeftable.jl`.
- **VASP torque target ↔ the model torque convention** (`io/dftsource.jl`, `io/vasp.jl`):
  the training torque from a constrained-noncollinear OSZICAR is `τ_a = −m_a × B_a`
  (`B` = constraining field), which must stay the *same* physical quantity, sign, and
  `3×n_atoms` config/atom/`xyz` layout as the model's `predict_torque = e_a × ∂E/∂e_a`
  (the design-matrix convention) — flip one side and the co-fit silently biases. Moments
  and field are rotated from the `SAXIS` frame by `Rz(α)·Ry(β)`. **DFT-code I/O is confined
  to `AbstractDFTSource` adapters** (namespaced submodules like `VASP`); the SCE pipeline
  consumes only `SpinDatum`/`SCEDataset` and stays DFT-code-agnostic. The readers are
  cross-checked against Magesty's parsers in the oracle.
- **Sunny export conversion ↔ the energy reconstruction** (`sce/sunny.jl`,
  `ext/MagestyRebuildSunnyExt.jl`): `_l1_pair_matrix` / `_l2_onsite_matrix` must satisfy
  `eₐ'·M·e_b = Σ folded·Z·Z` (the gate is the `Z₁`/`Z₂` contraction test); the per-bond
  matrix is `jϕ·(4π)^(N/2)·M` and the two directed members `(a,b,R)`/`(b,a,−R)` fold into
  one matrix on the canonical `a≤b` bond (reverse transposed). The whole chain is checked
  **without Sunny** by `_reconstruct_energy ≈ predict_energy − j0`; the extension then sets
  `J = M/(SₐS_b)` so the `Sunny.System` energy matches. Only `ls=[1,1]`/`ls=[2]` are
  representable — every other channel must be **reported as skipped**, never silently
  dropped. Change a harmonic normalization or the `(4π)^(N/2)` scale → both the matrix
  formulas and the energy gate move together.
- `solve_coefficients(est, X, y)` receives a **column-centered** `X` (⇒ the solver
  adds no intercept; `j0` is recovered analytically in `fit`). Every estimator —
  in-tree or in an extension — must honor this.

## Tests

| Command | Purpose |
|---|---|
| `julia --project -e 'using Pkg; Pkg.test()'` | unit + Aqua (default) |
| `TEST_MODE=all julia --project -e 'using Pkg; Pkg.test()'` | unit + Aqua + JET |
| `TEST_MODE=jet julia --project -e 'using Pkg; Pkg.test()'` | JET type-stability |
| `julia --project=test/oracle test/oracle/runtests.jl` | from-scratch numerics vs pinned Magesty |

The core suite (`runtests.jl`) dispatches on the `TEST_MODE` env var
(`default`/`all`/`unit`/`aqua`/`jet`) and never depends on Magesty. The oracle
suite is a separate environment that `dev`s a pinned Magesty.jl.

## Git (this rebuild)

Local git (`add` / `commit` / branch) is pre-authorized for this exploratory
rebuild — no per-action confirmation. Only remote operations (`push`,
registration) require confirmation. Branch off the default branch before
committing.

## References

- `STYLE_GUIDE.md` — package-specific style deltas.
- `SPEC.md` — realized architecture, types, public API.
- `docs/design-notes.md` — why the rebuild diverges from Magesty (the refinements).
- `CHANGELOG.md` — what landed in the v0 slice.
- `examples/heisenberg_chain.jl` — runnable end-to-end (recovers `J`);
  `examples/kagome_threebody.jl` — 3-body / multi-term SALCs, energy+torque co-fit.
- `references/` — supporting literature (notes tracked, PDFs local-only).
