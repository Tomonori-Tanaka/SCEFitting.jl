<!-- Thanks for contributing to SCEFitting.jl. -->

## Summary

<!-- 1-3 sentences: what does this PR do and why. -->

## Type of change

- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New feature (non-breaking change that adds functionality)
- [ ] Breaking change (fix or feature that would change existing behavior)
- [ ] Docs / tests / chore only

## Numerical correctness

If this PR changes any numerical result:

- [ ] I have explained *why* the result changes.
- [ ] I have added a regression or validation test whose expected value is
      independent of the implementation (or a labeled change-detector pin).
- [ ] I have updated `docs/` / `examples/` if user-facing behavior changes.

If this PR touches physics conventions (signs, units, normalization, SALC
projection, Clebsch–Gordan, spherical harmonics, resolvability):

- [ ] I have consulted the theory pages of the
      [documentation](https://tomonori-tanaka.github.io/SCEFitting.jl/dev/).
- [ ] I have updated all coupled sites listed under "Coupled code sites" in
      [CLAUDE.md](../CLAUDE.md).

## Checks

- [ ] `make test-all` passes locally (4 threads).
- [ ] `make test-pin` passes, or the pins were recaptured with the reason recorded
      in `test/pin/PIN.md`.
- [ ] `make docs` builds (strict).
- [ ] Commit messages follow
      [Conventional Commits](https://www.conventionalcommits.org/).
- [ ] Public API changes are reflected in `SPEC.md` and `docs/src/api.md`.
- [ ] `CHANGELOG.md` `[Unreleased]` updated.

## Related issues

<!-- Closes #123, refs #456. -->
