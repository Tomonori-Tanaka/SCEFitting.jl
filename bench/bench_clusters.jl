# Cluster enumeration — neighbour list + symmetry-orbit reduction.
#
#   julia --project=bench bench/bench_clusters.jl [n] [nbody]
#
# The N≥3 clique check (all C(N,2) edges simultaneously at minimum image) is the
# costly part, so try `nbody = 3` on a larger `n`.

using SCEFitting
import Spglib
include(joinpath(@__DIR__, "fixtures.jl"))

n     = argn(1, 2)
nbody = argn(2, 2)

bench_header("cluster enumeration — bcc Fe $(n)×$(n)×$(n), nbody=$nbody")
cr  = bcc_fe(n)
itr = interaction(; nbody = nbody, cutoff = 2.6)
println("atoms = $(n_atoms(cr))")

sg = analyze_symmetry(SpglibBackend(), cr)

bench_one("build_neighbor_list",
          () -> build_neighbor_list(cr, itr.pair_cutoff, MinimumImage()))
nl = build_neighbor_list(cr, itr.pair_cutoff, MinimumImage())
println("→ neighbour pairs = $(length(nl.pairs))")

bench_one("build_clusters",
          () -> build_clusters(cr, nl, sg; nbody = nbody, selection = MinimumImage()))
cs = build_clusters(cr, nl, sg; nbody = nbody, selection = MinimumImage())
println("→ orbits by body = $(sort(collect(k => length(v) for (k, v) in cs.by_body)))")
