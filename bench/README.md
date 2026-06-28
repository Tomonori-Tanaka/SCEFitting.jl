# SCEFitting.jl benchmarks

Performance benchmarks for the core hotspots. This is **separate from
`test/oracle/`** — the oracle checks *numerical* agreement with Magesty.jl; these
scripts measure *speed and allocations*. Each script is standalone and runs in the
`bench/` environment (its own `Project.toml`/`Manifest.toml`, which develops the
package from `..`).

## Setup (once)

```bash
julia --project=bench -e 'using Pkg; Pkg.develop(path="."); \
    Pkg.add(["BenchmarkTools", "Spglib", "Statistics", "Printf", "Random", "LinearAlgebra"]); \
    Pkg.instantiate()'
```

## Run

```bash
julia --project=bench bench/bench_salcbasis.jl       [n] [lmax]      # SALC projection (the hotspot)
julia --project=bench bench/bench_clusters.jl        [n] [nbody]     # neighbour list + orbit reduction
julia --project=bench bench/bench_design_matrix.jl   [n] [m] [lmax]  # energy + torque design matrices
julia --project=bench bench/bench_solver.jl                          # OLS vs Ridge solve, size sweep
julia --project=bench bench/bench_end_to_end.jl      [n] [m] [lmax]  # SCEBasis build + fit
```

Positional arguments are optional (defaults in each script header). `n` is the bcc
Fe supercell size: `n = 2 → 16` atoms, `n = 3 → 54`, `n = 4 → 128`. Start at `n = 2`
to smoke-test, then scale up (the SALC build at `n = 4`, `lmax = 2` is the heavy
case — minutes, not seconds).

## Fixtures and harness

`fixtures.jl` (shared, `include`d by every script) provides:

- `bcc_fe(n; a = 2.87)` — an `n×n×n` bcc Fe supercell (`Crystal`),
- `rand_configs(crystal, m; seed)` — `m` random unit-column spin configs,
- `interaction(; nbody, cutoff, lmax, isotropy)` — the matching `Interaction`,
- `bench_one(label, f; ntrials)` — warm-up + `@timed`/`@allocations` (median/min ms,
  alloc count, MiB), for hot paths too heavy for `@benchmark` sampling,
- `bench_header(title)`, `argn(i, default)`.

## Recording results

Append a before/after entry to [`../.claude/bench_log.md`](../.claude/bench_log.md)
when you touch a hot path (`basis/salcbasis.jl`, `clusters/`, `sce/model.jl` design
kernels, `fitting/estimators.jl`, `basis/Harmonics.jl`). Note machine, Julia version,
thread count, and the fixture size (`n`, `lmax`, `m`).
