@testset "domain in image" begin
    left, right = extrema(domain_in_image(LOGIT, 0..10.0))
    @test left ≈ LOGISTIC(0.0)
    @test right ≈ LOGISTIC(10.0)
end

@testset "Univariate transformations" begin
    test_univariate(LOGIT)
    test_univariate(LOGISTIC)
    test_univariate(LOG)
    test_univariate(EXP)
    test_univariate(ODDSRATIO)
    test_univariate(INVODDSRATIO)
    for _ in 1:100
        test_univariate(Power(abs(randn())), 10)
    end
    for _ in 1:10
        a = Affine(randn(), randn())
        test_univariate(a, 10)
    end
end

