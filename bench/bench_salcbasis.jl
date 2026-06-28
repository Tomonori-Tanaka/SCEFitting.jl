# SALC basis construction — the dominant hotspot (symmetry-adapted projection).
#
#   julia --project=bench bench/bench_salcbasis.jl [n] [lmax]
#
# `n` is the bcc Fe supercell size (default 2 → 16 atoms; 4 → 128 atoms), `lmax`
# the per-species angular cutoff (default 2, i.e. up to single-ion anisotropy).

using SCEFitting
import Spglib                                   # activate the SpglibBackend extension
include(joinpath(@__DIR__, "fixtures.jl"))

n    = argn(1, 2)
lmax = argn(2, 2)

bench_header("build_salc_basis — bcc Fe $(n)×$(n)×$(n), lmax=$lmax")
cr  = bcc_fe(n)
itr = interaction(; nbody = 2, cutoff = 2.6, lmax = lmax)
println("atoms = $(n_atoms(cr))")

# Stage the inputs the SALC build consumes (these are cheap relative to it).
sg = analyze_symmetry(SpglibBackend(), cr)
nl = build_neighbor_list(cr, itr.pair_cutoff, MinimumImage())
cs = build_clusters(cr, nl, sg; nbody = itr.nbody, selection = MinimumImage())

bench_one("build_salc_basis", () -> build_salc_basis(cr, sg, cs;
              lmax_by_species = itr.lmax, isotropy = itr.isotropy))

b = build_salc_basis(cr, sg, cs; lmax_by_species = itr.lmax, isotropy = itr.isotropy)
println("→ n_salcs = $(length(b))")
