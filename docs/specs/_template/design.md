# Design: <title>

Status: draft (YYYY-MM-DD)

## Summary

<!-- 1-3 paragraphs describing the chosen design. List alternatives only
     when you need to explain why they were rejected. -->

## Module layout

<!-- Affected files / modules / types and their responsibilities. -->

| Target | Change |
|---|---|
| `src/<layer>/<file>.jl` | <!-- e.g., change the signature of `bar` --> |

## API

<!-- Public / internal API additions, changes, deletions. Note type
     annotations, docstring style, and which keywords are deliberately
     defaultless. -->

```julia
# Example
function new_fn(x::Foo; opt::Bool = false)::Bar
```

## Types and conventions

<!-- Impact on physics conventions, units, numerical conventions. List any
     new invariants. If this diverges from SLCE.jl, the ledger row. -->

## Impact on coupled sites

<!-- Which "Coupled code sites" in CLAUDE.md does this touch? -->

- [ ] `Harmonics` / `AngularMomentum` / SALC projection / normalization tests:
- [ ] `SALCKey` order / design columns / TOML persistence / `coeftable`:
- [ ] Energy kernel ↔ torque kernel (`sce/model.jl`):
- [ ] Decor engine ↔ pointed moment basis ↔ `momentfit.jl` doors:
- [ ] Readers' `zero_moment_atol` ↔ dataset doors:
- [ ] Upstream divergence ledger (SLCE.jl):
- [ ] `.claude/agents/` references:
- [ ] `SPEC.md` / `docs/src/api.md` updates:

## Test strategy

<!-- New tests, modified tests, benchmarks. For every new gate name the
     ORACLE (closed form / hand calculation / invariant / independent
     implementation / labeled pin). -->

## Risks and open items

<!-- Anything that may change numerical results, unresolved decisions,
     deferred alternatives. -->
