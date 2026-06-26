using Test

const TEST_MODE = get(ENV, "TEST_MODE", "default")

@testset "MagestyRebuild.jl" begin
    if TEST_MODE in ("default", "all", "unit")
        include("unit/test_geometry.jl")
    end
    if TEST_MODE in ("default", "all", "aqua")
        include("aqua.jl")
    end
    if TEST_MODE in ("all", "jet")
        include("jet.jl")
    end
end
