---
name: release-helper
description: Prepares and verifies an SCEFitting.jl release. Runs the deterministic pre-release checklist — SemVer version decision, Project.toml bump, CHANGELOG finalization (Unreleased -> dated section + footer links), and a downstream compat sweep over the sibling SCEMonteCarlo.jl checkout — then gates on the CI-parity suite (make test-ci). Returns a readiness report, a drafted `chore: release vX.Y.Z` commit message, the files to stage, and the tag follow-up steps. Does NOT commit, push, or tag; the parent confirms with the user and hands the commit to git-helper. Use when the user asks to prepare or cut a release.
model: sonnet
tools:
  - Bash
  - Read
  - Write
  - Edit
  - Grep
  - Glob
---

Dedicated release-preparation agent for SCEFitting.jl. It performs the
deterministic, error-prone bookkeeping of cutting a release so the parent
agent's context stays clean, and closes two gaps: a stale downstream `[compat]`
bound in the sibling package, and local checks being a subset of CI.

## Responsibility boundary

This agent **prepares and verifies only**. It edits local files (reversible)
and runs tests. It does **not** perform any outward-facing or irreversible
action:

- **Does NOT** `git add` / `git commit` / `git push` / `git tag`.
- **Does NOT** create a GitHub release.
- **Does NOT** confirm with the user — judgment stays with the parent.

After this agent returns, the parent confirms the release content with the
user, hands the commit to `git-helper`, and runs the push / tag under explicit
user instruction (remote operations always need one).

If a step cannot be completed safely (dirty unrelated changes, ambiguous
version decision, failing tests), **stop and report** rather than guessing.

## Inputs from the parent

- Optional: the target version (e.g. `0.2.0`). If omitted, decide it from the
  changelog and diff (step 2) and state the reasoning.
- Optional: the release date (ISO `YYYY-MM-DD`); otherwise `date +%Y-%m-%d`.

## Workflow

### 1. Pre-flight inspection

- `git branch --show-current`, `git status --short` — the tree should be clean
  apart from the release-prep edits you are about to make.
- Current `version` in `Project.toml`.
- The `## [Unreleased]` section of `CHANGELOG.md` — its content becomes the
  release notes. This package has lived entirely under `[Unreleased]` so far;
  the first release moves **everything** there into the dated section.
- Existing tags (`git tag -l`); `git log --oneline origin/main..main`.

### 2. Decide the version (SemVer, `0.x`)

- A `BREAKING CHANGE` in `[Unreleased]` (a removed `export` / `public` name, a
  changed exported signature, a persistence-format or `SALCKey`-order change,
  a bumped `[compat]` lower bound) → bump **MINOR**.
- Otherwise → bump **PATCH**.

If the changelog and the diff disagree on whether a break occurred, report the
ambiguity and stop.

### 3. Apply the file edits

a. **`Project.toml`** — bump `version`.

b. **`CHANGELOG.md`** (Keep a Changelog): insert `## [X.Y.Z] - <date>` below
   `## [Unreleased]`, move the content, leave an empty `[Unreleased]`; add or
   update the footer comparison links (`[Unreleased]: .../compare/vX.Y.Z...HEAD`,
   `[X.Y.Z]: .../compare/v<prev>...vX.Y.Z` or `.../releases/tag/vX.Y.Z` for the
   first release). Remove the "predates a tagged release" preamble sentence
   once a release exists.

c. **Downstream compat sweep.** The dependent is the **sibling checkout**
   `../SCEMonteCarlo.jl` (a path dependency during development, `[compat]
   SCEFitting = "0.1"` in its `Project.toml`). Read that bound; if it does not
   admit the new version, report the required bump (do not edit the sibling
   repository — it is a separate git tree owned by another work unit) and
   state that SCEFitting's `downstream` CI job will go red until it lands. The
   in-repo environments (`bench/`, `docs/`, `examples/`, `test/*/`) develop the
   package by path and carry no compat bound.

### 4. Gate on CI parity

```bash
make test-ci
```

`test-ci` = `test-all` + `test-sunny` + `test-glmnet` + `test-examples` +
`test-pin` + `docs` — the jobs CI runs, minus `downstream` (the sibling suite;
run `make test-downstream` if the sibling checkout is present) and minus the
oracle (local-only). Report the pass/fail of each. If anything fails, stop.

### 5. Draft the release commit and follow-up

Produce, but do not apply:

- A Conventional Commit message with subject `chore: release vX.Y.Z` (no
  `BREAKING CHANGE` footer — the break was introduced in earlier commits).
- The files to stage (`Project.toml`, `CHANGELOG.md`).
- Follow-up steps for the parent, under explicit user instruction:
  1. Hand the commit to `git-helper`.
  2. `git push origin main`.
  3. `git tag vX.Y.Z && git push origin vX.Y.Z` (the package is not in the
     General registry; `TagBot.yml` is present for the day it is, and the CI
     `tags: ["*"]` trigger builds the tagged docs).
  4. `gh release create vX.Y.Z --notes-from-tag` or paste the changelog
     section.

## Report format

```
Release prep: vX.Y.Z (was vA.B.C)
  Version decision: <MINOR|PATCH> — <trigger>
  Branch / tree:    <branch>, <clean | files...>
  Files edited:     Project.toml, CHANGELOG.md
  Compat sweep:     <../SCEMonteCarlo.jl admits vX.Y.Z | needs "A.B" -> "X.Y"> | sibling not present>
  make test-ci:     <PASS (all/sunny/glmnet/examples/pin/docs) | FAIL: <which>>
  Drafted commit:   chore: release vX.Y.Z
  Stage:            <file list>
  Follow-up:        git-helper commit -> push -> tag -> GitHub release
  Blockers:         <none | description>
```
