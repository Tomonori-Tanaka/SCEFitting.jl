# End-to-end — the two coarse wall-time indicators a user actually feels:
# building the basis (`SCEBasis` = symmetry + clusters + SALC), and fitting
# (energy-only and energy+torque co-fit).
#
#   julia --project=bench bench/bench_end_to_end.jl [n] [m] [lmax] [cutoff]

using SCEFitting
import Spglib
include(joinpath(@__DIR__, "fixtures.jl"))

n      = argn(1, 4)
m      = argn(2, 100)
lmax   = argn(3, 2)
cutoff = argf(4, 6.0)

bench_header("end-to-end — bcc Fe $(n)×$(n)×$(n), $m configs, lmax=$lmax, " *
             "cutoff=$cutoff")
cr   = bcc_fe(n)
spec = basis_spec(; nbody = 2, cutoff = cutoff, lmax = lmax)
println("atoms = $(n_atoms(cr))")

bench_one("SCEBasis (sym+clusters+SALC)",
          () -> SCEBasis(cr, spec; backend = SpglibBackend()))

b    = SCEBasis(cr, spec; backend = SpglibBackend())
cfgs = rand_configs(cr, m)
rng  = MersenneTwister(2)
tqs  = [randn(rng, 3, n_atoms(cr)) for _ = 1:m]  # random targets — timing only

bench_one("SCEDataset (energy)", () -> SCEDataset(b, cfgs, randn(rng, m)))
bench_one("SCEDataset (energy+torque)", () -> SCEDataset(b, cfgs, randn(rng, m), tqs))

ds  = SCEDataset(b, cfgs, randn(rng, m))
dst = SCEDataset(b, cfgs, randn(rng, m), tqs)
println("n_salcs = $(n_salcs(b))   energy rows = $m   " *
        "torque rows = $(length(dst.y_T))")

bench_one("fit (OLS, energy)", () -> fit(SCEFit, ds, OLS()))
bench_one("fit (OLS, co-fit w=0.5)",
          () -> fit(SCEFit, dst, OLS(); torque_weight = 0.5))
