# Contributing to SCEFitting.jl

Thanks for your interest in contributing. SCEFitting.jl fits spin-cluster
expansion (SCE) models — energies, torques, and adiabatic site moments — to
noncollinear DFT data. Numerical correctness, reproducibility, and physical
consistency take precedence over stylistic refactoring.

## Reporting issues

- **Bugs**: open a GitHub issue using the *Bug report* template. Please include
  a minimal reproducer (a small `input.toml` or `Crystal` + `BasisSpec`,
  expected vs observed numbers).
- **Feature requests**: use the *Feature request* template. Describe the
  scientific use case so we can judge scope.
- **Security issues**: see [SECURITY.md](SECURITY.md). Do not file public
  issues for vulnerabilities.

## Development workflow

1. Fork the repository and create a topic branch from `main` (`fix/<slug>`,
   `feat/<slug>`, `refactor/<slug>`, `docs/<slug>`, …).
2. Make changes. For non-trivial work (multiple design decisions, days of
   effort, behavioral changes) we use spec folders under `docs/specs/`; a
   template is at [docs/specs/_template/](docs/specs/_template/).
3. Add or update tests. Numerical changes must come with a regression or
   validation test whose expected value is **independent of the
   implementation** (closed form, hand calculation, invariant, independent
   implementation); captured-output pins are allowed only when labeled as
   change detectors (see `test/pin/PIN.md`).
4. Run the local checks before opening a PR:
   ```bash
   make test-all      # unit + Aqua + JET (4 threads — required)
   make test-pin      # change detectors over the SALC chain
   make docs          # strict Documenter build (executes every example)
   ```
   `make test-ci` runs the full CI matrix (adds the Sunny / GLMNet extension
   suites and the runnable examples); `make ci-local` reproduces CI from a
   cold environment.
5. Update documentation as needed:
   - User-facing changes → `docs/src/` (guides, tutorials)
   - New public API → `SPEC.md` and `docs/src/api.md`
   - Examples (`examples/`) that exercise the changed code path
   - `CHANGELOG.md` `[Unreleased]`

## Commit messages

We follow [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>(<scope>): <subject>
```

- Types: `feat` / `fix` / `docs` / `test` / `refactor` / `perf` / `chore` / `style`
- Subject: imperative mood, lowercase, no trailing period
- Breaking changes: include `BREAKING CHANGE: ...` in the commit body
- Backports from the sibling SLCE.jl cite the upstream commit SHA

## Style

- US English in source, comments, docstrings, commit messages, and PRs.
- See [STYLE_GUIDE.md](STYLE_GUIDE.md) for the package-specific rules
  (argument order, `SALCKey`, inner constructors) on top of the shared Julia
  style.
- Hot-path guidance (StaticArrays, threading, bounds checks) lives in
  [CLAUDE.md](CLAUDE.md) under "Performance guidelines"; benchmark fixtures
  and the regression rule in [bench/README.md](bench/README.md) and
  [bench/BENCH_LOG.md](bench/BENCH_LOG.md).

## Physics conventions

These are easy to break silently — confirm before touching the algorithm:

- Spin directions are unit vectors; layout `3 × n_atoms`; moment magnitudes
  are a separate quantity.
- Real tesseral spherical harmonics `Zₗₘ` with the `(4π)^(N/2)` design scale
  applied exactly once.
- Only even-`Σl` (time-reversal-even) labels; improper rotations through the
  per-site `(−1)^l` parity.
- Torque `τ_a = −e_a × ∂E/∂e_a` (Landau–Lifshitz sign).
- Minimum-image / Wigner–Seitz periodic resolvability — aliases are never
  folded into a shorter shell.
- Energy unit follows the DFT input (eV); `j0` is stored separately from `Jφ`.

The full list is in [CLAUDE.md](CLAUDE.md) ("Numerical / physics conventions")
and the theory pages of the [documentation](https://tomonori-tanaka.github.io/SCEFitting.jl/dev/).

## Pull requests

- One logical change per PR. Small follow-ups are preferred over giant PRs.
- Fill in the PR template; check the boxes that apply.
- CI must pass: tests on both operating systems, the extension suites, the
  examples, the downstream SCEMonteCarlo.jl suite, the pin tier, and the docs
  build.

## License

By contributing, you agree that your contributions will be licensed under the
[MIT License](LICENSE).
