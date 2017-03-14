@testset "composed transformation calculations" begin
    a = Affine(1.0,2.0)
    test_univariate(LOGIT ∘ a)
    test_univariate(a ∘ LOGIT)
    test_univariate(LOGISTIC ∘ a)
    test_univariate(a ∘ LOGISTIC)
    test_univariate(ODDSRATIO ∘ a)
    test_univariate(a ∘ ODDSRATIO)
end

@testset "composed transformation domains" begin
    @test_throws Exception domain(LOGIT ∘ Affine(2.0, 5.0) ∘ LOGISTIC)
end

@testset "bridge default test" begin
    @test bridge(0..1, -1..1) == Affine(2, -1)
end

function bridge_complex_test(dom, img, mapping = nothing)
    t = if mapping == nothing
        bridge(dom, img)
    else
        bridge(dom, mapping, img)
    end
    @test domain(t) == dom
    @test image(t) == img
    left, right = extrema(dom)
    xs = vcat([left], sort(collect(rand(dom) for _ in 1:10000)), [right])
    ys = t.(xs)
    @test all(y ∈ img for y in ys)
    @test issorted(ys)
    ymin, ymax = extrema(ys)
    yleft, yright = extrema(img)
    @test ymin == yleft
    @test ymax == yright
end

@testset "bridge complex test" begin
    bridge_complex_test(0..∞, -1.0..1.0)
    bridge_complex_test(0..∞, ℝ)
    bridge_complex_test(-1.0..1.0, 0..∞)
    bridge_complex_test(ℝ, 0..∞)
end
