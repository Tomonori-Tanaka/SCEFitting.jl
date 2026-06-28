using SCEFitting
using JET

@testset "JET" begin
    JET.test_package(SCEFitting; target_modules = (SCEFitting,))
end
