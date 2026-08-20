using Test

# Construction internals are public but unexported (the `SCEBasis` constructor drives them
# for you); the suite reaches them by qualification. Brought into the test module here so
# every included unit file sees them without each repeating the import.
using SCEFitting: build_neighbor_list, build_clusters, build_salc_basis, evaluate_salc,
    salcs, SALC, SALCKey, SALCBasis, NeighborPair, NeighborList, ClusterMember,
    ClusterOrbit, ClusterSet, analyze_symmetry, n_ops, SymOp, SpaceGroup,
    interplanar_spacing, solve_coefficients

const TEST_MODE = get(ENV, "TEST_MODE", "default")

# A misspelled TEST_MODE used to select NO branch below and exit green with a total of
# zero assertions — CI sets this value from YAML, so a typo there silently deletes the
# entire suite. [Backported from SLCE.jl a1ac9af.]
const _TEST_MODES = ("default", "all", "unit", "aqua", "jet")
TEST_MODE in _TEST_MODES || error(
    "TEST_MODE=\"$TEST_MODE\" is not one of $(join(_TEST_MODES, ", ")) — refusing to " *
    "run zero tests and report success")

# Several gates compare a threaded result against a SERIAL reference; at one thread
# `Threads.@threads` *is* the serial reference, so those testsets pass while asserting
# nothing (the threaded-vs-serial design builds in test_threading.jl, the deterministic
# basis builds in test_salc.jl / test_nbody.jl). This used to be a `@warn` that scrolled
# past 28k assertions, plus a line in CI.yml. Refuse instead.
# [Backported from SLCE.jl a1ac9af.]
if Threads.nthreads() == 1 && get(ENV, "SCEFITTING_ALLOW_SINGLE_THREAD", "0") != "1"
    error("run the suite with `-t N` for N > 1 (CI pins 4): the threaded-vs-serial " *
          "gates are vacuous at one thread. Set SCEFITTING_ALLOW_SINGLE_THREAD=1 to " *
          "override.")
end

@testset "SCEFitting.jl" begin
    if TEST_MODE in ("default", "all", "unit")
        include("unit/testutils.jl")     # shared helpers (rand_unit / randcfg /
                                         # canonical-member gates)
        include("unit/test_geometry.jl")
        include("unit/test_harmonics.jl")
        include("unit/test_angmom.jl")
        include("unit/test_coupledbasis.jl")
        include("unit/test_symmetry.jl")
        include("unit/test_clusters.jl")
        include("unit/test_imageselection.jl")
        include("unit/test_truncation.jl")
        include("unit/test_ws_nbody.jl")
        include("unit/test_salc.jl")
        include("unit/test_normalization.jl")
        include("unit/test_fit.jl")
        include("unit/test_selection.jl")
        include("unit/test_validation.jl")
        include("unit/test_torque.jl")
        include("unit/test_nbody.jl")
        include("unit/test_persist.jl")
        include("unit/test_input.jl")
        include("unit/test_dftsource.jl")
        include("unit/test_dataset.jl")
        include("unit/test_embset.jl")
        include("unit/test_coeftable.jl")
        include("unit/test_sunny.jl")
        include("unit/test_introspect.jl")
        include("unit/test_threading.jl")
    end
    if TEST_MODE in ("default", "all", "aqua")
        include("aqua.jl")
    end
    if TEST_MODE in ("all", "jet")
        include("jet.jl")
    end
end
