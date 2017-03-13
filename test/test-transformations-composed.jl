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

@testset "bridge tests" begin
    t = bridge(0..∞, InvOddsRatio(), -1.0..1.0)
    dom = domain(t)
    img = image(t)
    @test dom == 0..∞
    @test img == -1.0..1.0
    xs = vcat([0], sort(collect(rand(dom) for _ in 1:10000)), [∞])
    ys = t.(xs)
    @test all(y ∈ img for y in ys)
    @test issorted(ys)
    ymin, ymax = extrema(ys)
    @test ymin == -1
    @test ymax == 1
end
