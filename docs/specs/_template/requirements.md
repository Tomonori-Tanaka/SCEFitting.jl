# Requirements: <title>

Status: draft (YYYY-MM-DD)

## Goal

<!-- What we want to achieve. 1-3 sentences. -->

## Background

<!-- Why now. Constraints, requester, prior investigations, the upstream
     SLCE.jl state if this is a backport. -->

## Scope

Includes:

- <!-- e.g., API change in `fitting/` -->

Excludes:

- <!-- e.g., refactor of `clusters/` (separate spec) -->

## Invariants

<!-- Things that must NOT change: physics conventions, numerical conventions,
     public API, persistence format, etc. -->

- Spin layout stays `3 × n_atoms` with unit columns; `magmoms` separate.
- Real tesseral `Zₗₘ`; the `(4π)^(N/2)` design scale applied exactly once.
- `SALCKey` column addressing and ordering; existing TOML models reload and
  predict identically.
- Torque sign `τ = −e × ∂E/∂e`; co-fit whitening.
- Minimum-image / Wigner–Seitz resolvability (no alias folding).
- Threaded ≡ serial bitwise gates keep passing.
- ...

## Completion criteria

- [ ] <!-- e.g., `make test-all` passes -->
- [ ] <!-- e.g., new gate with an implementation-independent oracle -->
- [ ] <!-- e.g., new API reflected in `docs/src/api.md` and `SPEC.md` -->

## References

- Related issues / PRs:
- Related specs / design notes / upstream SLCE.jl commits:
