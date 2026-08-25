# `docs/specs/` index

Folders for mid-sized or larger development units. Operating rules live in
the "Managing development units" section of [CLAUDE.md](../../CLAUDE.md).
When starting a new spec, copy from [`_template/`](_template/).

## Status table

| Spec | Status | One-line summary |
|---|---|---|
| [260824-moment-toml-input](260824-moment-toml-input/) | landed (2026-08-24) | Optional `[moment]` TOML section → `MomentSpec`; `MomentBasis(path)` |
| [260824-penalty-metric](260824-penalty-metric/) | landed (2026-08-25) | Basis-intrinsic penalty metric (both channels), unpenalized μ₀, λ selection for `MomentFit` |
| [260824-pointed-nbody-general](260824-pointed-nbody-general/) | landed (2026-08-25) | General body order in the pointed enumeration; `nbody` cap raised to 4, `cutoff_star` per star order |
| [260825-moment-lsum-per-body](260825-moment-lsum-per-body/) | landed (2026-08-25) | Per-body-order `lsum` in the pointed moment spec (`[interaction].lsum`'s spelling) |

(Work before this index was introduced — the v0 rebuild, the SLCE.jl
carve-out, and the pointed-moment backport — is recorded in `CHANGELOG.md`
and in the design record referenced there; it is not back-filled here.)

This table and the `Status:` line in each `tasklist.md` are duplicated
intentionally; update both when a spec lands. Reconcile in bulk when drift is
found.

## Completion criteria

Each `tasklist.md` ends with a shared exit checklist (see
[`_template/tasklist.md`](_template/tasklist.md)). In particular,
**`.claude/agents/` is easy to forget** — whenever module names or Makefile
targets change, sweep the agent files as well.
