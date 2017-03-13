@testset "Univariate integral" begin
    f, dom = integral_substitution(InvOddsRatio(), x->exp(-x^2), 0..Inf)
    @test hquadrature(f, dom.left, dom.right)[1] ≈ √π/2
end
