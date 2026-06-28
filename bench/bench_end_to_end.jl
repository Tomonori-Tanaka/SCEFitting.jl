# End-to-end — the two coarse wall-time indicators a user actually feels:
# building the basis (`SCEBasis` = symmetry + clusters + SALC), and fitting.
#
#   julia --project=bench bench/bench_end_to_end.jl [n] [m] [lmax]

using SCEFitting
import Spglib
include(joinpath(@__DIR__, "fixtures.jl"))

n    = argn(1, 2)
m    = argn(2, 30)
lmax = argn(3, 2)

bench_header("end-to-end — bcc Fe $(n)×$(n)×$(n), $m configs, lmax=$lmax")
cr  = bcc_fe(n)
itr = interaction(; nbody = 2, cutoff = 2.6, lmax = lmax)
println("atoms = $(n_atoms(cr))")

bench_one("SCEBasis (sym+clusters+SALC)",
          () -> SCEBasis(cr, itr; backend = SpglibBackend()))

b    = SCEBasis(cr, itr; backend = SpglibBackend())
cfgs = rand_configs(cr, m)
ds   = SCEDataset(b, cfgs, randn(m))            # random targets — timing only
println("n_salcs = $(n_salcs(b))")

bench_one("fit (OLS)", () -> fit(SCEFit, ds, OLS()))
