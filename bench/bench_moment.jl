# Pointed moment basis — the body-order cost curve, stage by stage.
#
#   julia --project=bench bench/bench_moment.jl [n] [cutoff_pair] [star3] [star4]
#
# `n` is the bcc Fe supercell size (default 3 → 54 atoms), `cutoff_pair` the 2-body
# mark–environment radius in Å (default 4.1 = 3NN), `star3` / `star4` the star radii of
# body orders 3 and 4 (defaults 4.1 and 2.5 — the 4-body probe is only affordable on
# the first shell, which is what per-order radii exist for).
#
# The five stages are timed SEPARATELY because they scale differently with N: member
# generation is `C(z, N−1)·N!`, orbit reduction is a fold over those members, the SALC
# projection is an `eigen` on a per-block carrier space, resolvability is a symbolic
# rank, and the design matrix is a per-column sweep. Reporting one build time hides
# which of them is the wall.
#
# TTFX (first call from a cold process) is measured in a child process, because the
# parent has already compiled everything by the time it gets there.

using SCEFitting
import Spglib                                   # activate the SpglibBackend extension
include(joinpath(@__DIR__, "fixtures.jl"))

using SCEFitting: MomentSpec, MomentBasis, moment_resolvability, penalty_metric,
                  _design_moment, _pointed_star_candidates, _orbits_from_members,
                  _star_cutoff, _star_cutoff_envelope, n_salcs, n_ops,
                  AbstractSymmetryBackend

# The symmetry analysis is not what this script measures, and re-running Spglib inside
# every timed `MomentBasis` call would fold a constant into stage (c). Hand the finished
# space group back through the backend seam instead.
struct FixedSG <: AbstractSymmetryBackend
    sg::Any
end
SCEFitting.analyze_symmetry(b::FixedSG, ::Crystal; tol::Real = 1e-5) = b.sg

n     = argn(1, 3)
cpair = argf(2, 4.1)
star3 = argf(3, 4.1)
star4 = argf(4, 2.5)

_spec(N) = MomentSpec(; lmax_env = [2], sampled = [true], lmax_mark = 2, nbody = N,
                      cutoff_pair = cpair, lsum = 4, isotropy = true,
                      cutoff_star = N == 3 ? star3 : [star3, star4])

# ---- TTFX child: one cold MomentBasis build, nothing else ---------------------------
if get(ENV, "SCEFIT_BENCH_MOMENT_TTFX", "") != ""
    N = parse(Int, ENV["SCEFIT_BENCH_MOMENT_TTFX"])
    t0 = time_ns()
    cr = bcc_fe(n)
    sg = analyze_symmetry(SpglibBackend(), cr)
    MomentBasis(cr, _spec(N); backend = FixedSG(sg))
    @printf("ttfx N=%d  %8.2f s (process start → first MomentBasis)\n",
            N, (time_ns() - t0) / 1e9)
    exit(0)
end

bench_header("MomentBasis — bcc Fe $(n)×$(n)×$(n), cutoff_pair=$cpair, " *
             "cutoff_star=[$star3, $star4]")
cr = bcc_fe(n)
sg = analyze_symmetry(SpglibBackend(), cr)
println("atoms = $(n_atoms(cr))   space-group ops = $(n_ops(sg))")
cfgs = rand_configs(cr, 8; seed = 20260825)

for N in (3, 4)
    println("\n── body order $N " * "─"^54)
    spec = _spec(N)
    # (a) member generation, at THIS order's radius (the shared neighbour list is built
    #     once at the envelope of the per-order radii, exactly as the ctor does it).
    nl = build_neighbor_list(cr, _star_cutoff_envelope(spec), MinimumImage(); tol = 1e-8)
    bench_one("(a) star members", () -> _pointed_star_candidates(cr, nl, spec, N))
    stars = _pointed_star_candidates(cr, nl, spec, N)
    @printf("    members = %d   (radius %.2f Å)\n",
            length(stars), _star_cutoff(spec, N)[1, 1])

    # (b) orbit reduction over those members
    bench_one("(b) orbit reduction", () -> _orbits_from_members(cr, sg, stars, N))
    orbs = _orbits_from_members(cr, sg, stars, N)
    @printf("    orbits  = %d\n", length(orbs))

    # (c) the whole build. The projection is what is left after (a) and (b), and for
    #     the orders that matter it is nearly all of it — printed as a share so a
    #     regression in either half is visible without a second harness.
    tot = bench_one("(c) MomentBasis (a+b+projection)",
                    () -> MomentBasis(cr, spec; backend = FixedSG(sg)))
    mb = MomentBasis(cr, spec; backend = FixedSG(sg))
    @printf("    n_salcs = %d   (body %d: %d)\n", n_salcs(mb), N,
            count(k -> k.body == N, mb.salc_basis.keys))

    # (d) the periodic-resolvability gate — run once per MomentDataset construction,
    #     so its cost is paid on every dataset, not only at build time. The default-rtol
    #     result is CACHED on the basis, so timing that call would report the cache;
    #     pass the default value explicitly to take the same work uncached.
    bench_one("(d) moment_resolvability", () -> moment_resolvability(mb; rtol = 1e-10))
    res = moment_resolvability(mb)
    @printf("    rank = %d / %d kept\n", res.rank, length(res.kept))

    # (e) the two per-fit sweeps over the finished basis
    bench_one("(e) _design_moment (8 cfg)", () -> _design_moment(mb, cfgs, cfgs))
    bench_one("(e) penalty_metric (256 cfg)",
              () -> penalty_metric(mb; nconfig = 256, seed = 1))
end

# ---- (f) TTFX, one cold process per body order --------------------------------------
println("\n── (f) TTFX " * "─"^54)
for N in (3, 4)
    run(setenv(`$(Base.julia_cmd()) --project=$(@__DIR__) $(@__FILE__) $n $cpair $star3 $star4`,
               "SCEFIT_BENCH_MOMENT_TTFX" => string(N)))
end
