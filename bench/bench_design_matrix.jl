# Design-matrix assembly — energy (`evaluate_salc`) and torque (`accumulate_grad!`)
# kernels over all SALC columns × configurations.
#
#   julia --project=bench bench/bench_design_matrix.jl [n] [m] [lmax]
#
# `n` supercell size, `m` number of spin configurations (default 30), `lmax` angular
# cutoff. The torque block is the heavier of the two (gradient + cross product per
# atom, 3·n_atoms rows per config).

using SCEFitting
import Spglib
include(joinpath(@__DIR__, "fixtures.jl"))

n    = argn(1, 2)
m    = argn(2, 30)
lmax = argn(3, 2)

bench_header("design matrices — bcc Fe $(n)×$(n)×$(n), $m configs, lmax=$lmax")
cr   = bcc_fe(n)
itr  = interaction(; nbody = 2, cutoff = 2.6, lmax = lmax)
b    = SCEBasis(cr, itr; backend = SpglibBackend())
cfgs = rand_configs(cr, m)
println("atoms = $(n_atoms(cr))   n_salcs = $(n_salcs(b))   configs = $m")

# `_design_energy` / `_design_torque` are internal kernels (not exported).
bench_one("_design_energy", () -> SCEFitting._design_energy(b, cfgs))
bench_one("_design_torque", () -> SCEFitting._design_torque(b, cfgs))
