# Benchmark log — SCEFitting.jl

A running record of performance numbers for the core hotspots, kept so refactors can
be checked for regression after the fact. Scripts live in [this directory](.).

> **History before 2026-08-11 lives in `SLCE.jl/bench/BENCH_LOG.md`.** This package
> was revived on 2026-08-11 by cutting `698a841` of the package that is now
> SLCE.jl, so every optimisation entry dated before that (the 128-atom SALC build,
> the design-matrix halving, the `dnPl` cache) was measured on the shared ancestor
> and applies here too. Rather than copy them, this log starts at the carve-out.

**How to log.** When you touch a hot path (`basis/salcbasis.jl`, `clusters/`, the
`sce/model.jl` design kernels, `fitting/estimators.jl`, `basis/Harmonics.jl`), run the
relevant script in this directory before and after and append an entry below:

- a one-line context header (date, branch/commit, machine, Julia version, threads,
  fixture size `n`/`lmax`/`m`),
- a **Before / After** table (median time, MiB, allocs),
- a short note (wall-time vs alloc trade-off, downstream impact, follow-ups).

Entries are append-only history — keep them after merging.

---

## Penalty metric: construction cost and λ-path neutrality — 2026-08-24

**Context**: 2026-08-24 · `main` · local macOS (darwin 24.6, aarch64) · julia 1.12.7 ·
**threads = 4** · `bench/bench_solver.jl` (new sections) · fixture: bcc Fe 2×2×2
(16 atoms), `nbody = 2`, `cutoff = 4.1`, `lmax = 2`, Spglib backend → **18 SALCs**.

Not a before/after: `penalty_metric` is new, and the question it has to answer is
whether attaching it by default in `Ridge(basis; ...)` / `GroupAdaptiveRidge(basis; ...)`
is affordable, and whether it slows the λ path down.

### Reference-ensemble construction

| `nconfig` | `torque_weight` | time | allocated |
|---|---|---|---|
| 500 | 0.0 | 33.2 ms | 48.1 MiB |
| 500 | 1.0 | 102.9 ms | 86.2 MiB |
| 2000 | 0.0 | 130.9 ms | 192.5 MiB |
| 2000 | 1.0 | 439.6 ms | 344.9 MiB |

Linear in `nconfig`, as it must be. The torque form costs ~3.4× the energy form (it adds
a gradient evaluation and an `n_atoms` loop per configuration). At the default
`nconfig = 2000` this is a **0.1–0.4 s one-off** on this basis; it scales with the column
count, so budget proportionally for a production basis (108 columns at 3×3×3 bcc Fe).

`@allocated` is CUMULATIVE allocation, dominated by the evaluation kernel's per-call
churn — the same churn `_design_energy` pays — not by anything the function holds. The
accumulation that matters for memory (never materializing the
`nconfig · 3 · n_atoms × p` torque design, 280 MiB at 3×3×3 bcc Fe) is invisible on a
fixture this small; it is a structural property of the loop, not a measured one here.

### λ path

| | time |
|---|---|
| `select_fit`, uniform penalty | 0.7 ms |
| `select_fit`, basis metric | 0.7 ms |

Ratio 0.95 — no measurable difference, which is the expected result: the metric is one
extra multiply per column per IRLS step and changes no loop structure. The cost of the
metric is entirely in its construction, and a caller who wants to avoid it can pass
`metric = nothing` or reuse a vector across estimators.

**Follow-up**: the moment channel has no benchmark script at all
(`bench/` covers clusters, design matrices, SALC build, the solver, and two end-to-end
fixtures), so `penalty_metric(::MomentBasis)` is unmeasured. Worth adding alongside a
`bench_moment.jl`.

---

## Baseline for the pointed-moment backport — 2026-08-21 (phase 0, S-bench)

**Why this entry exists.** The backport plan needs a performance baseline that
predates the work, plus an explicit rule for what counts as a regression. This
entry is that baseline; the rule is below it. Nothing here is a claim about
correctness.

**Context**: 2026-08-21 · branch `main` · local macOS (darwin 24.6, aarch64) ·
julia 1.12.6 · **threads = 1** (the scripts' default) · all six scripts at their
**stress defaults** (no positional arguments).

### The regression gates — and why exactly these

| script | role | why |
|---|---|---|
| `bench_salcbasis` | **GATE** | the projector / fold hot path — the code the backport touches |
| `bench_design_matrix` | **GATE** | `evaluate_salc` per (row × column) — the other half of the same chain |
| `bench_clusters` | context | upstream of the basis; the backport does not touch orbit enumeration |
| `bench_end_to_end` | context | a sum of the two gates plus the solve; moves for reasons the gates already explain |
| `bench_nd2fe14b` | context | realistic many-orbit / few-ops shape — a sanity read, too coarse to gate |
| `bench_solver` | context | pure BLAS; independent of everything the backport changes |

**How to judge.** Run each gate script **five times** and take the median of the
reported medians.

- **Allocation count is the primary gate, at zero tolerance.** Measured
  2026-08-21: three consecutive runs of each gate reported allocation counts and
  MiB that were *bit-for-bit identical* — the quantity is deterministic, so any
  change at all is a real change and must be explained in the entry that causes
  it.
- **Wall time is the secondary gate, at 5%.** Measured run-to-run spread over
  three runs was at most **1.4%** (salcbasis) and under 1.1% (design matrices),
  so 5% carries roughly 3.5x headroom over this machine's noise. A threshold at
  the noise floor would be a coin flip, not a gate.
- A gate that trips is not automatically a veto — it is a required line in the
  commit's log entry, saying what was traded for what.

### Measured — SCEFitting.jl @ `116c4fd`, SLCE.jl @ `a132928`

Same machine, same session, same fixtures. Both packages emit **the same
`n_salcs`** on every fixture below, so these are like-for-like.

| script / target | fixture | SCEFitting med | SLCE med | SCEFitting allocs | SLCE allocs |
|---|---|---|---|---|---|
| `build_salc_basis` | bcc Fe 4x4x4 (128 at.), lmax 3, cutoff 6.0 | **2803 ms** | **1998 ms** | 70,872,871 | 60,778,374 |
| `_design_energy` | bcc Fe 4x4x4, 100 cfg, lmax 2 | 466 ms | 474 ms | 11,299,224 | 11,299,408 |
| `_design_torque` | bcc Fe 4x4x4, 100 cfg, lmax 2 | 991 ms | 1026 ms | 7,527,192 | 7,527,376 |
| `build_neighbor_list` | bcc Fe 4x4x4, nbody 3, cutoff 6.0 | 3.66 ms | 3.71 ms | 24 | 24 |
| `build_clusters` | bcc Fe 4x4x4, nbody 3, cutoff 6.0 | 450 ms | 483 ms | 10,390,312 | 10,390,312 |
| basis build (end-to-end) | bcc Fe 4x4x4, lmax 2 | 934 ms | 846 ms | 18,103,347 | 16,419,513 |
| dataset (energy) | bcc Fe 4x4x4, 100 cfg | **465 ms** | **802 ms** | 11,300,233 | 17,759,959 |
| dataset (energy+torque) | bcc Fe 4x4x4, 100 cfg | **1452 ms** | **1845 ms** | 18,827,426 | 25,287,343 |
| `fit` (OLS, energy) | 46 columns | 0.20 ms | 0.15 ms | 320 | 326 |
| `fit` (OLS, co-fit w=0.5) | 46 columns, 38400 torque rows | 25.6 ms | 29.1 ms | 341 | 351 |
| `build_salc_basis` | Nd2Fe14B 68 at., nbody 3, cutoff 4.0 | 386 ms | 347 ms | 13,408,478 | 12,434,646 |
| `_design_energy` | Nd2Fe14B, 103 cfg | 120 ms | 121 ms | 1,723,810 | 1,725,398 |
| `_design_torque` | Nd2Fe14B, 103 cfg | 275 ms | 299 ms | 1,099,600 | 1,101,188 |
| `solve_coefficients` OLS | M=10000, P=800 | 2236 ms | 2117 ms | — | — |
| `solve_coefficients` Ridge | M=10000, P=800 | 126.6 ms | 126.7 ms | — | — |

**Two differences are large enough to matter to the backport route choice:**

1. **`build_salc_basis` is ~29% faster in SLCE** (1998 vs 2803 ms; 60.8M vs 70.9M
   allocations) at identical output (`n_salcs = 129` both, no drops on this
   fixture), and ~10% faster on Nd2Fe14B. Whatever SLCE gained since the carve-out
   was never backported. Under a facade route this comes for free; under a
   verbatim-port route it has to be ported or written off.
2. **Dataset construction is ~1.7x SLOWER in SLCE** (802 vs 465 ms energy-only).
   SLCE's dataset does strictly more work — the resolvability classification and
   the ASR reparametrisation, which SCEFitting does not have. This is a real cost
   of the facade route, and it is paid once per dataset, not per fit.

### Repairs made while taking the baseline

- **`bench/bench_solver.jl` was broken in BOTH packages**: it calls
  `solve_coefficients`, which is `public` but *not* exported, so `using SLCE` /
  `using SCEFitting` alone left it undefined and the script died on line 13.
  Fixed with an explicit `using <Pkg>: solve_coefficients`. The script had been
  dead since the export/public split; nothing runs `bench/` in CI, so nothing
  caught it.

---

## Slot-based SALC term layout (D2) — 2026-08-21

**Context**: branch `feat/pointed-moment` · baseline `6967364` (pre-D2) measured on
the same machine in the same session as the after numbers · local macOS (darwin
24.6, aarch64) · julia 1.12.6 · threads = 1 · stress defaults.

`SALCTerm` moved from a per-site `ls::Vector{Int}` to a per-axis
`slots::Vector{Slot}` (a `Slot` is `(site::Int, factor::SiteFactor)`, so four
words per axis where an `l` was one).

### GATE — `bench_salcbasis` (bcc Fe 4x4x4, 128 atoms, lmax 3, cutoff 6.0)

| | before (`6967364`) | after (D2) | delta |
|---|---|---|---|
| median | 2806 ms | **2914 ms** | **+3.8 %** |
| allocs | 70,873,129 | **73,263,913** | **+2,390,784 (+3.4 %)** |
| memory | 4256.84 MiB | 4389.00 MiB | +3.1 % |
| `n_salcs` | 129 | 129 | — |

**The allocation gate trips (zero tolerance), deliberately.** What was traded:
the term label is now channel-aware, which is the whole point of the slice, and
it costs one `Vector{Slot}` per raw term at construction plus one more in
`_canonicalize_members`' remap, where the old layout allocated one
`Vector{Int}`. Wall time stays inside the 5 % secondary gate.

Note the scale: the finished basis holds 137,472 terms over 109,568 members, so
the +2.39 M allocations are dominated by the *pre-reduction* term population,
not by the surviving one.

### GATE — `bench_design_matrix` (bcc Fe 4x4x4, 100 cfg, lmax 2)

| | before (`116c4fd`) | after (D2) | delta |
|---|---|---|---|
| `_design_energy` median | 466 ms | 469 ms | +0.6 % |
| `_design_energy` allocs | 11,299,224 | **11,299,224** | **0** |
| `_design_torque` median | 991 ms | 1030 ms | +3.9 % |
| `_design_torque` allocs | 7,527,192 | **7,527,192** | **0** |

**The evaluation half is allocation-neutral to the byte.** Reading `l` and the
site through a slot instead of a parallel `ls` vector costs nothing in the
kernels — the indirection is resolved at compile time.

### One micro-optimisation tried and rejected

`_function_vector` (the orbit-reduction key builder) allocates one
`Vector{Int}` per term for `_term_spin_ls(t)`. Keying the dictionary on
`t.slots` directly removes that allocation — measured 73,263,913 -> 72,988,981 —
but costs **+44 % wall time** (2914 -> 4206 ms), because `hash(::Slot)` goes
through `objectid` and that dictionary is probed once per nonzero tensor entry.
Reverted; the comment in `_function_vector` records the measurement so nobody
re-tries it.

---

## Mixed-channel (decor) SALC engine (D3) — 2026-08-21

**Context**: branch `feat/pointed-moment` · baseline = the D2 entry above,
measured on the same machine · local macOS (darwin 24.6, aarch64) · julia 1.12.6 ·
threads = 1 · stress defaults.

The production pure-spin path is untouched by the new engine. Two shared pieces
did move: `SALCScratch` gained an `rl` buffer (the `SolidHarmonics` batch
workspace the displacement axes need) and `AngularMomentum.build_real_bases`
gained a `keep` path predicate, called once per coupling path with a default
that always accepts.

### GATE — `bench_salcbasis` (bcc Fe 4x4x4, 128 atoms, lmax 3, cutoff 6.0)

| | D2 | D3 | delta |
|---|---|---|---|
| median | 2914 ms | 2941 ms | +0.9 % |
| allocs | 73,263,913 | **73,263,913** | **0** |
| `n_salcs` | 129 | 129 | — |

**Allocation-neutral to the byte.** The `keep` predicate costs nothing: with the
default closure the call inlines away, and the basis build never constructs a
`SALCScratch`.

### GATE — `bench_design_matrix` (bcc Fe 4x4x4, 100 cfg, lmax 2)

| | D2 | D3 | delta |
|---|---|---|---|
| `_design_energy` median | 469 ms | 478 ms | +1.9 % |
| `_design_energy` allocs | 11,299,224 | **11,299,316** | **+92** |
| `_design_torque` median | 1030 ms | 1025 ms | −0.5 % |
| `_design_torque` allocs | 7,527,192 | **7,527,284** | **+92** |

**The allocation gate trips (zero tolerance) by exactly 92 on each, and the 92
is fully attributed**: `_design_energy` / `_design_torque` construct one
task-local `SALCScratch` per column and the bench has 46 SALCs, so the new `rl`
field is built 46 times. Adding one `Vector{Float64}(undef, 4)` field to a
struct constructed 46 times costs exactly 92 allocations (measured directly on a
two-field vs three-field toy struct: 140 -> 232). It is per *scratch*, not per
term or per entry — 8e-6 of the total — and it is what buys a single evaluation
workspace shared by both channels instead of a second one allocated per term.
Wall time stays inside the 5 % secondary gate on both.
