# Tasklist: <title>

Status: draft (YYYY-MM-DD)

This file holds coarse-grained, commit-sized milestones. Day-to-day tracking
goes through `TaskCreate` in-session.

## Milestones

### M1 — <name>

- [ ] <!-- coarse step -->

### M2 — <name>

- [ ] <!-- ... -->

## Exit checklist

Run through every item once implementation lands. ~~Strike through~~ items
that do not apply.

- [ ] `make test-all` passes (4 threads).
- [ ] `make test-pin` passes, or pins recaptured with the reason in
      `test/pin/PIN.md`.
- [ ] `make docs` builds (strict).
- [ ] If results changed: regression or validation test added, oracle
      independent of the implementation.
- [ ] If public API changed: `SPEC.md` and `docs/src/api.md` updated.
- [ ] If a hot path was touched: before / after recorded in
      `bench/BENCH_LOG.md`.
- [ ] Tier 2 review panel run (numerical / maintainability / performance /
      API axes) and findings resolved.
- [ ] If module names or Makefile targets changed: `.claude/agents/` swept.
- [ ] If this diverges from SLCE.jl: divergence ledger row in `CLAUDE.md`.
- [ ] `CHANGELOG.md` `[Unreleased]` updated.
- [ ] `Status:` line in this file and the table in `docs/specs/README.md`
      updated in sync.
- [ ] Implementation commit hash appended below.
