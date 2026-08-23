# Tasklist: `[moment]` section in the TOML setup file

Status: landed (2026-08-24) — Q1–Q4 resolved as recommended; commit hash in tasklist.md

This file holds coarse-grained, commit-sized milestones. Day-to-day tracking
goes through `TaskCreate` in-session.

## Milestones

### M1 — reader + constructor (`feat(io): [moment] section and MomentBasis(path)`)

- [x] `src/SCEFitting.jl`: move `include("io/input.jl")` after `basis/momentbasis.jl`
      (and the I/O section comment with it).
- [x] `src/io/input.jl`: schema docstring, `_moment_from_input`,
      `_species_list_from_input`, `read_setup` の `moment` フィールド、
      `MomentBasis(path; backend, tol, tie_tol)`、`soc` / unknown-key refusal.
- [x] `src/basis/momentbasis.jl`: `MomentSpec` docstring に TOML 経路の 1 行。
- [x] `test/unit/test_input.jl`: design.md §Test strategy 1–5.
- [x] `make test-all` 緑 (38485)。

### M2 — docs (`docs: [moment] TOML schema`)

- [x] `docs/src/guide/io.md` スキーマ小節、`docs/src/guide/moment.md` TOML 経路
      （`docs/src/api.md` は bare `MomentBasis` が自動で拾うので編集しない）。
- [x] `SPEC.md` I/O 節（`tie_tol` + `moment`）、`CHANGELOG.md` `[Unreleased]`、
      `CLAUDE.md` の upstream ledger 行 + coupled-site bullet、`docs/specs/README.md` の行。
- [x] `make docs` strict 緑 (after listing `torque_weight_per_site` in api.md, a pre-existing omission from `d11d367`)。

### M3 — acceptance on bcc Fe (package 外、コミットしない)

- [x] `~/jijs/scefit/for_package/fe_bcc/3x3x3cubic/scefit/l024_c4.1/input.toml` に
      `[moment]` を追加、`fit.jl` を `MomentBasis("input.toml")` に簡略化、
      83 SALC / 同一 RMSE を確認（`m_p5.0_s4.1` 等の variant も `input.toml` だけの差に）。
      → 実施済み: 10 variant すべて変換、`l024_c4.1` 再実行で fit_summary.txt がヘッダ 2 行以外同一。

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~ items
that do not apply.

- [x] `make test-all` passes (4 threads).
- [x] `make test-pin` passes (85), no recapture.
- [x] `make docs` builds (strict).
- [ ] ~~If results changed: regression or validation test added, oracle
      independent of the implementation.~~ (no numerical change; the new gates
      use hand-written specs / documented defaults as oracles)
- [x] If public API changed: `SPEC.md` updated (`api.md` picks the new method up via the bare binding).
- [ ] ~~If a hot path was touched: before / after recorded in
      `bench/BENCH_LOG.md`.~~
- [x] Tier 2 review panel run — reduced to API + maintainability axes (Q4, agreed):
      API approve (3 minor, applied), maintainability 2 major (Bool <: Integer leak,
      `soc` test match string) + 7 minor, all applied; Tier 1 `code-reviewer` final pass:
      0 major, 2 minor (empty species array door added).
- [ ] ~~If module names or Makefile targets changed: `.claude/agents/` swept.~~
- [x] If this diverges from SLCE.jl: divergence ledger row in `CLAUDE.md`.
- [x] `CHANGELOG.md` `[Unreleased]` updated.
- [x] `Status:` line in this file and the table in `docs/specs/README.md`
      updated in sync.
- [ ] Implementation commit hash appended below.
