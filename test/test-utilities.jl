"Test singleton type."
@define_singleton TestSingleton <: Real

@testset "define singleton" begin
    
    @test TestSingleton <: Real
    @test isa(TESTSINGLETON, TestSingleton)
    @test repr(@doc(TESTSINGLETON)) == repr(@doc(TestSingleton)) ==
        repr(doc"Test singleton type.")

end
