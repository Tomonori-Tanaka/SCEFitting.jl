# Tasklist: cross-orbit alias groups

Status: landed (2026-09-11)

## Milestones

### M1 — detection + ties in the basis

- [x] (2026-09-11) `src/basis/aliases.jl`: `AliasGroup`, `_ColumnTies`, `_column_ties`, `alias_groups`, info message
- [x] (2026-09-11) `include("basis/aliases.jl")` + `public alias_groups, AliasGroup, n_columns` in `src/SCEFitting.jl`
- [x] (2026-09-11) `SCEBasis.ties` + `alias_rtol` keyword + `n_columns`
- [x] (2026-09-11) gates 1, 2, 4, 5 of the test strategy; `test/unit/test_aliases.jl` registered in `test/runtests.jl`

### M2 — column space through the fit

- [x] (2026-09-11) `_fold_columns` in `design.jl`; `SCEPredictor(f)` expansion; `coef`/`SCEFit` docstrings
- [x] (2026-09-11) `salc_groups` / `group_costs` / `penalty_metric` / `select_fit` in column space
- [x] (2026-09-11) gates 3, 7

### M3 — disclosure

- [x] (2026-09-11) `SCEPredictor.split`; `coeftable` columns; read-out warnings; persist v7 + legacy read
- [x] (2026-09-11) gate 6

### M4 — docs and ledger

- [x] (2026-09-11) `theory/resolvability.md`, `SPEC.md`, `api.md`, `CHANGELOG.md`, `CLAUDE.md` (coupled site + ledger)
- [x] (2026-09-11) `build_salc_basis` docstring / `_reduce_orbit_salcs` comment wording
- [x] (2026-09-11) real-system check on Nd2Fe14B (l02 = 10 groups, 169/169 tied columns full rank; l06 = 30; l044_c4.0 = 20) recorded in CHANGELOG

## Exit checklist

- [x] (2026-09-11) `make test-unit` passes (4 threads, 39233); `make test-sunny` (26); Aqua/JET left to CI.
- [x] (2026-09-11) `make test-pin` passes unchanged (104 at 4 and at 1 threads).
- [x] (2026-09-11) `make docs` builds (strict; a pre-existing stale `@ref` in `build_neighbor_list` fixed on the way).
- [x] (2026-09-11) `test/unit/test_aliases.jl`: hand geometry, closed-form SALC value, hand bond-sum energy, bond-weighted mean (Cm); alias-free bases bitwise unchanged.
- [x] (2026-09-11) `SPEC.md`, `docs/src/api.md`, guides, README updated.
- [x] (2026-09-11) `bench/BENCH_LOG.md` entry (identity fold; `_column_ties` ≤ 4.4 % of the build on the only aliasing fixture).
- [x] (2026-09-11) Tier 2 panel run; all majors applied (Cm per-bond gate, `:unequal_norm`, `alias_rtol` cap, inner constructors, `_ColumnTies` naming, shared MGS helper, `_column_ties` accessor, doc sites). Deferred minors: tied-path `_fold_columns` peak memory, `:legacy` naming, per-load `@info`, `[interaction].alias_rtol` in the setup file.
- [ ] ~~If module names or Makefile targets changed: `.claude/agents/` swept.~~
- [x] (2026-09-11) Divergence-ledger row (SLCE freezes the shared orbits; this package ties and discloses).
- [x] (2026-09-11) `CHANGELOG.md` `[Unreleased]` updated.
- [x] (2026-09-11) `Status:` line and the `docs/specs/README.md` row updated in sync.
- [x] (2026-09-11) Implementation commit hash appended below.

### Landing

- Implementation commit: `6e3b1a1` (2026-09-11).
