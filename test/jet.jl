using MagestyRebuild
using JET

@testset "JET" begin
    JET.test_package(MagestyRebuild; target_modules = (MagestyRebuild,))
end
