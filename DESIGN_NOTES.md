# Design notes index

Index of design discussions, investigations, and on-hold ideas. Follow the
links for detail. Active development units live under `docs/specs/`;
historical benchmarks live in `bench/BENCH_LOG.md` and `git log`.

Operating rules: [`docs/design-notes/README.md`](docs/design-notes/README.md).

## Design proposals

| Topic | Status | Last update |
|---|---|---|

Completed proposals are folded into their corresponding spec under
`docs/specs/` and removed from this index.

## Investigations and standing rationale

- [Why the rebuild diverges from Magesty.jl](docs/design-notes.md) — the
  deliberate refinements (generalized neighbor list, Wigner–Seitz
  resolvability, Wigner-D from the package's own `Zₗₘ`, `SALCKey` addressing,
  …); standing reference, kept at its original path.
- [Upstream divergence ledger (SLCE.jl ↔ this package)](CLAUDE.md#upstream-divergence-ledger-slcejl--this-package)
  — the spellings and polarities deliberately kept different from the joint
  spin–lattice sibling; consult before any backport in either direction.

## Performance backlog

- `bench/BENCH_LOG.md` names the two gate scripts and their regression
  thresholds (2026-08-21 baseline entry); no open backlog items.
