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

## Pointed body order: the N = 3 → 4 cost curve, stage by stage — 2026-08-25

**Context**: 2026-08-25 · `main` · local macOS (darwin 24.6, aarch64) · julia 1.12.7 ·
**threads = 4** · `bench/bench_moment.jl` (new; `make bench-moment`) · fixture: bcc Fe
3×3×3 (54 atoms, 2592 space-group ops), `lmax_mark = lmax_env = 2`, `lsum = 4`,
`cutoff_pair = 4.1`, `cutoff_star = [4.1, 2.5]` (3-body to 3NN, 4-body to 1NN).

Not a before/after either: body order 4 did not exist before this spec, and **body
order 3 is unchanged bit-for-bit** (the pin tier holds that). The N = 3 column is the
regression baseline for later work; the N = 4 column is what raising the door costs.

Read the two columns carefully: rows **(a)** and **(b)** are that body order ALONE
(the stage is called once per order), while **(c)**–**(f)** are the whole `nbody = N`
build or the whole basis — so (c) at N = 4 contains the 3-body work of the N = 3 row.
"(a)+(b) is 2 % of (c)" below uses that reading.

| Stage | N = 3 (med) | N = 4 (med) | note |
|---|---|---|---|
| (a) star members | 12 ms / 38 MiB | 4 ms / 25 MiB | 76,788 vs 72,576 members |
| (b) orbit reduction | 26 ms / 45 MiB | 14 ms / 25 MiB | 14 vs 3 orbits |
| (c) whole `MomentBasis` | **2.15 s / 6.5 GiB** | **4.36 s / 13.3 GiB** | 83 vs 95 SALCs |
| (d) `moment_resolvability` | 1.54 s / 1.9 GiB | 1.88 s / 2.6 GiB | uncached (see below) |
| (e) `_design_moment`, 8 cfg | 63 ms / 61 MiB | 64 ms / 72 MiB | |
| (e) `penalty_metric`, 256 cfg | 1.47 s / 1.9 GiB | 1.81 s / 2.2 GiB | |
| (f) TTFX (cold process, start → exit) | 9.5 s | 13.0 s | child at 4 threads |

### Before / after: `Threads.@threads :greedy` (review panel, same session)

Three loops were on the default schedule, which cuts the index range into one
**contiguous** chunk per thread — and in all three the expensive iterations are
contiguous **at the end** (orbits and columns are emitted in ascending body order, and
both projection and column evaluation grow steeply with it). So every high-body item
landed in the last chunk and ran serially while the other threads idled. `:dynamic`
would not have helped: it is the same chunking, differing only in thread affinity.
`:greedy` hands out iterations individually. Each task writes only its own slot and the
output is re-sorted by key, so the result is **bitwise identical at any schedule** —
`make test-pin` is unchanged (104, no recapture) and the threaded ≡ serial gates stay
green.

| Loop | before | after |
|---|---|---|
| `MomentBasis` orbit loop → stage (c), N = 4 | 5.73 s | **4.36 s** (−24 %) |
| `MomentBasis` orbit loop → stage (c), N = 3 | 2.07 s | 2.15 s (+4 %, within this stage's noise — 14 orbits over 4 threads is barely skewed, and the stage allocates 6.5 GiB) |
| `penalty_metric(::SCEBasis)` column loop, N = 4 basis | 2.00 s | **1.81 s** (−10 %) |
| `_design_moment` column loop (was `:dynamic`) | 73 ms | 64 ms (−12 %) |

The N = 4 build is where it pays: 3 four-body orbits against 4 threads is the worst
case for contiguous chunking. A basis with more, more uniform orbits will see less.

**What the split says.** The wall is the **SALC projection**, not the enumeration:
(a) + (b) is under 2 % of (c) at either order. Raising the door to 4 costs **+2.2 s of
projection for 12 extra columns** — the Reynolds projector's `eigen(Symmetric(P))` runs
on a per-block carrier of dimension `assignments × paths × (2L_f+1)`, and the 4-body
label `(0; 1, 1, 2)` has three assignments where the 3-body ones have one. Design-matrix
evaluation, by contrast, barely moves (+12 %): once the basis exists, a 4-body column
costs about what a 3-body one does.

**Per-order radii are doing the work.** At 2.5 Å the 4-body star sees the 8 first
neighbours, so `C(8,3)·4! = 1344` members per marked atom — *fewer* members than the
3-body sector at 4.1 Å, where `C(z,2)·3!` runs over three shells. Held at a single
4.1 Å radius the same 4-body sector is 2,774,736 members and 115,614 P1 orbits (M0
measurement, 2026-08-24), which is why `cutoff_star` became per star order: the honest
4-body probe is unreachable without it.

**Measurement note.** `moment_resolvability(mb)` caches its default-`rtol` result on
the basis, so timing that call reports the cache (0.000 ms). Stage (d) passes
`rtol = 1e-10` — the same value, explicitly — to take the same work uncached. Anyone
re-running this must keep that, or the row silently becomes a cache-hit benchmark.

**The first TTFX take (12.1 s / 17.9 s) was discarded, not carried forward.** It was
invalid three ways, all in the child-process harness: `setenv` REPLACED the child's
environment rather than extending it, so the Makefile's `JULIA_NUM_THREADS = 4` was
dropped and `Base.julia_cmd()` does not carry `-t` — the child built single-threaded
inside a 4-thread table; and the clock was started inside the child, after `using
SCEFitting` and the Spglib extension load, so the window excluded exactly the package
load TTFX is about. The harness now uses `addenv`, passes `-t` explicitly, and times
the child from the parent. The 9.5 / 13.0 s in the table above are the corrected take,
and they are not comparable to the discarded pair.

**Follow-ups (not done here).** `_eval_term_mixed` is not `Val(D)`-specialized, and the
`N!` re-anchoring expands every translation class instead of carrying a multiplicity
weight; both were deliberately left alone because bit-identity at `N ≤ 3` came first.

## Penalty metric: construction cost and λ-path neutrality — 2026-08-24

**Context**: 2026-08-24 · `main` · local macOS (darwin 24.6, aarch64) · julia 1.12.7 ·
**threads = 4** · `bench/bench_solver.jl` (new sections) · fixture: bcc Fe 3×3×3
(54 atoms), `nbody = 2`, `cutoff = 4.1`, `lmax = 2`, Spglib backend → **21 SALCs**;
the λ path on 400 configurations over 25 λ.

Not a before/after: `penalty_metric` is new, and the questions it has to answer are
whether attaching it by default in `Ridge(basis; ...)` / `GroupAdaptiveRidge(basis; ...)`
is affordable, and whether it slows the λ path down.

### Reference-ensemble construction

| `nconfig` | `torque_weight` | time | allocated |
|---|---|---|---|
| 2048 | 0.0 | 0.56 s | 814 MiB |
| 2048 | 1.0 | 1.25 s | 651 MiB |
| 8192 | 0.0 | 2.36 s | 3256 MiB |
| 8192 | 1.0 | 5.29 s | 2603 MiB |

Linear in `nconfig`, as it must be, and linear in the column count. The torque form
costs ~2.2× the energy form (it adds a gradient evaluation and an `n_atoms` loop per
configuration). `@allocated` is CUMULATIVE allocation, dominated by the evaluation
kernel's per-call churn — the same churn `_design_energy` pays — not by anything the
function holds; the loop never materializes the `nconfig · 3 · n_atoms × p` torque
design, which would be 280 MiB of live memory here.

**This table set the default.** Paired with the accuracy measurement below, 2048 is
where cost and noise meet: 8192 would be a 5 s constructor on a 21-column basis and
proportionally worse on a production one, for a factor-2 reduction in a noise that is
already an order of magnitude below the systematic factor the metric removes.

### Accuracy of the estimate (relative standard error of `mⱼ` over 8 seeds)

Measured separately, on bcc Fe 2×2×2 with `lmax = 2` at `nbody = 2` (18 columns) and
`nbody = 3` (82 columns), as the across-seed spread divided by the across-seed mean:

| `nconfig` | 2-body median / max | 3-body median / max |
|---|---|---|
| 256 | 8.5 % / 16.0 % | 9.4 % / 16.4 % |
| 512 | 6.0 % / 9.7 % | 6.3 % / 11.2 % |
| 2048 | 3.3 % / 5.3 % | 3.3 % / 6.0 % |
| 8192 | 1.6 % / 2.8 % | 1.6 % / 2.8 % |

`1/√nconfig` throughout, and **body order barely moves it** — the concern that
higher-body columns have heavier-tailed `Φ²` and would need a larger ensemble is not
borne out at these orders.

### λ path

| | time |
|---|---|
| `select_fit`, uniform penalty | 2.5 ms |
| `select_fit`, basis metric | 2.5 ms |

Ratio 1.01 — no difference, which is the expected result: the metric is one extra
multiply per column per IRLS step and changes no loop structure. The cost of the metric
is entirely in its construction, which is why `SCEFitting.with_lambda` exists (move one
estimator along the path rather than rebuilding it per point).

**Follow-up**: the moment channel has no benchmark script at all (`bench/` covers
clusters, design matrices, SALC build, the solver, and two end-to-end fixtures), so
`penalty_metric(::MomentBasis)` — the chunked path, whose chunk size is a byte budget —
is unmeasured. Worth adding alongside a `bench_moment.jl`.

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
