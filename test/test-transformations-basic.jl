@testset "domain in image" begin
    left, right = extrema(domain_in_image(LOGIT, 0..10.0))
    @test left ≈ LOGISTIC(0.0)
    @test right ≈ LOGISTIC(10.0)
end

@testset "Univariate transformations" begin
    test_univariate(LOGIT, AD_exceptions = logit_exceptions())
    test_univariate(LOGISTIC, AD_exceptions = logistic_exceptions())
    test_univariate(LOG)
    test_univariate(EXP)
    test_univariate(ODDSRATIO)
    test_univariate(INVODDSRATIO)
    test_univariate(REALCIRCLE)
    test_univariate(INVREALCIRCLE)
    for _ in 1:100
        test_univariate(Power(randn()), N = 10)
    end
    for _ in 1:10
        test_univariate(Affine(randn(), randn()), N = 10)
    end
end
