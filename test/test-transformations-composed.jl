@testset "Composed transformations" begin
    test_univariate(Logit() ∘ Affine(1.0,2.0))
end


c = Logit() ∘ Affine(1.0,2.0)
