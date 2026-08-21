# Real-data parity tier (S6)

Runs the moment channel on the two real constrained-DFT sets the pointed-moment
design record was developed on, and compares SCEFitting against its upstream
SLCE.jl **column by column**. This is the package's "named system" check
(FeGe B20 2×2×2, FeRh B2 4×4×4), not a unit test.

**Not a persistent gate.** There is no CI job: the tier needs the sibling
`~/Packages/SLCE.jl` checkout and data that live outside the repository
(`~/jijs`, the design record's `step2_assets`). Run it by hand after anything
that touches `basis/momentbasis.jl`, `fitting/momentfit.jl`, the decor engine
or the extxyz / EMBSET readers:

```bash
julia --project=test/parity -t 4 test/parity/runtests.jl
```

`Manifest.toml` here is gitignored. Data locations are read from the
environment with defaults to the author's layout; a block whose data are absent
is **skipped loudly** (`@warn` + `@test_skip`), never silently passed:

| variable | default | used by |
|---|---|---|
| `SCE_PARITY_FEGE_DIR` | `~/Packages/_brain_storming/adiabatic-moment-sce/step2_assets/fege` | `tau01/03/05.extxyz` |
| `SCE_PARITY_FEGE_POSCAR` | `~/jijs/magesty/fege/2x2x2/lc_exp/ge_include/202601/mfa-tau01/300k5/nelect_0/input/POSCAR` | the 64-atom supercell |
| `SCE_PARITY_FERH_DIR` | `~/jijs/magesty/ferh/4x4x4/afm/tau_0to0.6/350k3` | `input/rh_random/{POSCAR,sample-NNN.INCAR}`, `output/rh_random/lambda100/EMBSET{,_mint}` |

## What it checks, and what each expected value is

Three kinds of statement, labeled as such in `runtests.jl`:

- **Reference numbers (acceptance)** — the FeGe `full` / `lsum2` in-sample σ,
  config-level 5-fold CV σ and held-out τ0.3 σ, and the FeRh gated σ, at
  **relative 1 %** against the values SLCE.jl produced on the same protocol
  (design record §3.3). They are an agreement statement between two
  implementations on real data, not physics oracles; the protocol (gate
  2.2e-4 μB, `Random.seed!(1)` + `randperm` folds) is the upstream script's.
- **Census integers (change detectors)** — kept-row counts and the FeRh mode-1
  antiparallel census (`n_anti` = 3833/7744 on the Rh orbit, 0/7744 on Fe).
  Exact equality; these pin what the archived data contain under the current
  mode rule and must be recaptured, with the reason recorded, if either changes.
- **Column parity (absolute)** — the pointed design columns of both packages on
  the same configurations, matched by `SALCKey`, compared with **no free scale**
  (`‖a − b‖ ≤ 1e-10 ‖b‖` per column and elementwise `atol = 1e-12·max|b|`),
  plus the dataset targets, gates and masks bitwise, and held-out predictions
  at 1e-10. Both readers are also compared on one extxyz file (bitwise fields).

The SLCE.jl revision actually compared against is printed at the top of the run
(a `[sources]` path cannot pin one) and belongs in any result you quote.
