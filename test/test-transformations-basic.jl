@testset "Univariate transformations" begin
    test_univariate(Logit())
    test_univariate(Logistic())
    test_univariate(Log())
    test_univariate(Exp())
    test_univariate(OddsRatio())
    test_univariate(InvOddsRatio())
    for _ in 1:100
        test_univariate(Power(abs(randn())), 10)
    end
    for _ in 1:10
        a = Affine(randn(), randn())
        test_univariate(a, 10)
    end
end
