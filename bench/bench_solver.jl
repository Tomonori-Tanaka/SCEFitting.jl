# Regression solve — `solve_coefficients` for OLS vs Ridge across matrix sizes.
#
#   julia --project=bench bench/bench_solver.jl
#
# Synthetic, column-centred design matrices (the solver contract: `X` is already
# centred and adds no intercept). Cheap enough for `@belapsed` sampling.

using SCEFitting
#  is public but unexported (the estimators drive it for you),
# so it has to be named explicitly --  alone does not bring it in.
using SCEFitting: solve_coefficients
import Spglib                                   # activates SpglibBackend
include(joinpath(@__DIR__, "fixtures.jl"))      # bench_header, @belapsed, mean, MersenneTwister

bench_header("solve_coefficients — OLS vs Ridge")

for (M, P) in ((2_000, 200), (5_000, 500), (10_000, 800))
    rng = MersenneTwister(0)
    X = randn(rng, M, P)
    X .-= mean(X; dims = 1)                     # column-centred (solver contract)
    y = randn(rng, M)
    t_ols   = @belapsed solve_coefficients(OLS(), $X, $y) samples = 3 evals = 1
    t_ridge = @belapsed solve_coefficients(Ridge(lambda = 1e-3), $X, $y) samples = 3 evals = 1
    @printf("M=%-6d P=%-5d   OLS=%8.3f ms   Ridge=%8.3f ms   ratio=%.2f\n",
            M, P, 1e3 * t_ols, 1e3 * t_ridge, t_ridge / t_ols)
end


# --- the penalty metric -----------------------------------------------------------
#
# `penalty_metric` runs in a CONSTRUCTOR (`GroupAdaptiveRidge(basis; ...)`), so its cost
# is user-visible before any fitting happens. It threads over columns and ACCUMULATES
# the torque block instead of materializing it, which is what keeps a large basis off
# the `nconfig · 3 · n_atoms × p` cliff (280 MiB at 3x3x3 bcc Fe / 108 columns / 2000
# configs — the fixture below is far too small to show that, so read the numbers here
# as the per-column-per-configuration cost, not as evidence for the accumulation).
# `@allocated` reports CUMULATIVE allocation, dominated by the evaluation kernel's churn
# rather than by anything this function holds.

bench_header("penalty_metric — reference ensemble (bcc Fe 2x2x2, lmax = 2)")

let
    cr = bcc_fe(2)
    b = SCEBasis(cr, basis_spec(; nbody = 2, cutoff = 4.1, lmax = 2);
                 backend = SpglibBackend())
    println("atoms = $(n_atoms(cr))   n_salcs = $(n_salcs(b))")
    for K in (500, 2000), w in (0.0, 1.0)
        t = @belapsed penalty_metric($b; torque_weight = $w, nconfig = $K) samples = 3 evals = 1
        mem = @allocated penalty_metric(b; torque_weight = w, nconfig = K)
        @printf("K=%-5d w=%.1f   %9.1f ms   %8.1f MiB\n", K, w, 1e3 * t, mem / 2^20)
    end
end

# --- a λ path with and without the metric -----------------------------------------
#
# The metric changes no inner loop (one extra multiply per column per IRLS step), so the
# path cost should be unchanged; what it adds is the one-off construction above.

bench_header("select_fit — λ path, uniform vs basis metric")

let
    cr = bcc_fe(2)
    b = SCEBasis(cr, basis_spec(; nbody = 2, cutoff = 4.1, lmax = 2);
                 backend = SpglibBackend())
    rng = MersenneTwister(7)
    cfgs = rand_configs(cr, 60)
    ds = SCEDataset(b, cfgs, randn(rng, 60))
    lams = 10.0 .^ range(-1, -7; length = 12)
    lw = SCEFitting.cost_weights(b; theta = 1.0)
    est_u = GroupAdaptiveRidge(lw.labels, lw.weights; lambda = 1.0)
    est_m = GroupAdaptiveRidge(b; lambda = 1.0, theta = 1.0, metric_nconfig = 500)
    t_u = @belapsed select_fit($ds, $est_u; lambdas = $lams) samples = 3 evals = 1
    t_m = @belapsed select_fit($ds, $est_m; lambdas = $lams) samples = 3 evals = 1
    @printf("uniform=%8.1f ms   metric=%8.1f ms   ratio=%.2f\n",
            1e3 * t_u, 1e3 * t_m, t_m / t_u)
end
