# Style guide

> The shared Julia style lives in `~/Packages/CLAUDE.md` ("Julia style"): official
> Julia + DFTK guides, naming, `for i = 1:n` vs `for x in xs`, explicit named
> tuples, ≤ 92 cols, hot-path `SVector`/`MVector`/`@views`/`@inbounds`. Only
> package-specific additions are listed here.

## Argument order

Follows the shared `(data) → spec → options → verbosity` shape, with one
carve-out:

- **Strategy-dispatched verbs put the strategy first** (so multiple dispatch is
  the public surface): `solve_coefficients(estimator, X, y)`,
  `read_configs(source)`, `analyze_symmetry(backend, crystal; tol)`.
- **Evaluation verbs keep predictor-first**: `predict_energy(model, data)`,
  `r2_energy(model, data)`.

## Package-specific naming

- `SALCKey` — the canonical structural identity of a basis column. Column order is
  `sort` over the keys; design-matrix columns are addressed by key, never by
  construction order.
- Hot-path harmonics come in a validated `Zₗₘ` and an unchecked, buffered
  `Zₗₘ_unsafe`.
- Rank parameter `R = N + 1` on `CoupledBasis{R}` (N = body order); `SALC` itself is
  not parametric — its terms store `folded::Array{Float64}` with runtime rank.
  **Never promote `l` / `l_max` to a type parameter** — they are runtime fields /
  loop bounds (avoids specialization blow-up).

## Inner constructors

Types whose outer constructor validates or normalizes inputs (e.g. `Crystal`
wraps fractional coordinates) define the logic in an **inner** constructor, so an
exactly field-typed call cannot bypass it to the default constructor.
