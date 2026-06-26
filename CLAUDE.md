# CLAUDE.md

> Shared baseline (numerical-correctness priority, JP-conversation / EN-repo
> policy, Conventional Commits, Julia style, shared subagents) is inherited from
> `~/Packages/CLAUDE.md`. Only package-specific rules live here.

## Project goal

Clean, extensible, Julia-native rebuild of `Magesty.jl`: fit a spin-cluster
expansion (SCE) `E = j0 + Σφ Jφ Φφ({e_a})` — with per-atom torque
`τ_a = −m_a (e_a × B_a)` co-fit — to noncollinear DFT data. The numerical core
(tesseral harmonics, Clebsch–Gordan coupling, symmetry-adapted basis, design
matrices, regression) is **reimplemented from scratch**; `Magesty.jl` is used only
as a *pinned numerical oracle* in `test/oracle/`. Priority: numerical correctness
and reproducibility over stylistic concerns. See `SPEC.md` for architecture.

## Numerical / physics conventions

Easy to break silently — confirm before touching the algorithm.

- **Spin directions are unit vectors**; `directions[:, a]` (norm 1) is the
  direction, the moment magnitude `magmom[a]` is stored separately.
- **Spin layout `3 × n_atoms`** (rows x, y, z; columns atoms). Transposing breaks
  the pipeline.
- **Real (tesseral) spherical harmonics `Zₗₘ`** (Drautz, PRB 102 024104), per-site
  factor `(4π)^(−1/2)`; the design matrix carries `(4π)^(N/2)` so an N-body term
  cancels to O(1).
- **Time reversal**: even-`Σl_s` channel only; the projector averages
  `½(D(g) + (−1)^{Σl_s} D(g)ᵀ)`.
- **Lattice / reciprocal**: `Lattice.vectors` columns are `aᵢ`; `reciprocal =
  inv(vectors)` rows are `bᵢ` with `aᵢ·bⱼ = δᵢⱼ` (no 2π). Interplanar spacing
  `dᵢ = 1/‖row_i(reciprocal)‖`. Fractional coords are wrapped to `[0,1)` on
  periodic axes (neighbor-list precondition).
- **Energy units**: `Jφ` carry the DFT input unit (eV); `j0` is separate.

## Coupled ("linked") code sites — change one, check all

- `basis/Harmonics.jl` (`Zₗₘ`, `∂ᵢZlm`) ↔ the on-sphere central-difference and
  closed-form agreement tests. Normalization / sign drift silently biases `X`.
- **Design-matrix columns are addressed by `SALCKey`, asserted at every join**
  (basis ↔ dataset ↔ model). A column's meaning is its canonical key, never
  construction order. (This replaces Magesty's implicit "SALC order is
  load-bearing" coupling — there is no `column == SALC iteration index` rule.)
- `assemble_weighted_problem` invariants (column-centered ⇒ solver adds no
  intercept; row-scaled by `torque_weight` ⇒ solver does not re-weight) ↔ every
  `solve_coefficients` method.

## Tests

| Command | Purpose |
|---|---|
| `julia --project -e 'using Pkg; Pkg.test()'` | unit + Aqua (default) |
| `TEST_MODE=all julia --project -e 'using Pkg; Pkg.test()'` | unit + Aqua + JET |
| `TEST_MODE=jet  julia --project -e 'using Pkg; Pkg.test()'` | JET type-stability |
| `TEST_MODE=oracle` (under `test/oracle/`) | from-scratch numerics vs pinned Magesty |

`runtests.jl` dispatches on the `TEST_MODE` env var.

## Git (this rebuild)

Local git (`add` / `commit` / branch) is pre-authorized for this exploratory
rebuild — no per-action confirmation. Only remote operations (`push`,
registration) require confirmation. Branch off the default branch before
committing.

## References

- `STYLE_GUIDE.md` — package-specific style deltas.
- `SPEC.md` — architecture, types, public API.
- `references/` — supporting literature (notes tracked, PDFs local-only).
