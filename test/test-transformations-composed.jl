@testset "composed transformation calculations" begin
    a = Affine(1.0,2.0)
    test_univariate(Logit() ∘ a)
    test_univariate(a ∘ Logit())
    test_univariate(Logistic() ∘ a)
    test_univariate(a ∘ Logistic())
    test_univariate(OddsRatio() ∘ a)
    test_univariate(a ∘ OddsRatio())
end

@testset "composed transformation domains" begin
    @test_throws Exception domain(Logit() ∘ Affine(2.0, 5.0) ∘ Logistic())
end
