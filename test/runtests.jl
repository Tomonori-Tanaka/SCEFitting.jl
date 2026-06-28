using Test

const TEST_MODE = get(ENV, "TEST_MODE", "default")

@testset "MagestyRebuild.jl" begin
    if TEST_MODE in ("default", "all", "unit")
        include("unit/test_geometry.jl")
        include("unit/test_harmonics.jl")
        include("unit/test_angmom.jl")
        include("unit/test_coupledbasis.jl")
        include("unit/test_symmetry.jl")
        include("unit/test_clusters.jl")
        include("unit/test_imageselection.jl")
        include("unit/test_ws_nbody.jl")
        include("unit/test_salc.jl")
        include("unit/test_fit.jl")
        include("unit/test_torque.jl")
        include("unit/test_nbody.jl")
        include("unit/test_persist.jl")
        include("unit/test_input.jl")
        include("unit/test_coeftable.jl")
        include("unit/test_sunny.jl")
        include("unit/test_vaspio.jl")
        include("unit/test_site_engine.jl")
    end
    if TEST_MODE in ("default", "all", "aqua")
        include("aqua.jl")
    end
    if TEST_MODE in ("all", "jet")
        include("jet.jl")
    end
end
