import Compat: ∘

@testset "composed transformation calculations" begin
    a = Affine(1.5,2.0)
    test_univariate(LOGIT ∘ a, AD_exceptions = logit_exceptions(inv(a)))
    test_univariate(a ∘ LOGIT, AD_exceptions = logit_exceptions())
    test_univariate(LOGISTIC ∘ a, AD_exceptions = logistic_exceptions())
    test_univariate(a ∘ LOGISTIC, AD_exceptions = logistic_exceptions())
    test_univariate(ODDSRATIO ∘ a)
    test_univariate(a ∘ ODDSRATIO)
    test_univariate(a ∘ a)
    isidentity(a ∘ inv(a))
    isidentity(LOGISTIC ∘ LOGIT)
end

@testset "composed transformation show" begin
    a = Affine(1,2)
    b = LOGISTIC
    @test sprint(show, a ∘ b) == sprint(show, a) * " ∘ " * sprint(show, b)
end

@testset "composed transformation domains" begin
    @test_throws Exception domain(LOGIT ∘ Affine(2.0, 5.0) ∘ LOGISTIC)
end

@testset "bridge default test" begin
    @test bridge(0..1, -1..1) == Affine(2, -1)
end

"""
Test for bijections between domain `dom` and image `img`.

When `RR`, domain and image of the actual transformation will be infinite.
"""
function bridge_complex_test(dom, img; mapping = nothing, RR = false)
    t = if mapping == nothing
        @inferred bridge(dom, img)
        bridge(dom, img)
    else
        # @inferred bridge(dom, mapping, img) BROKEN
        bridge(dom, mapping, img)
    end
    @test domain(t) == (RR ? ℝ : dom)
    @test image(t) == (RR ? ℝ : img)
    left, right = extrema(dom)
    xs = vcat([left], sort(collect(rand(dom) for _ in 1:10000)), [right])
    ys = t.(xs)
    @test all(y ∈ img for y in ys)
    @test issorted(ys, rev = !isincreasing(t))
    ymin, ymax = extrema(ys)
    yleft, yright = extrema(img)
    @test ymin == yleft
    @test ymax == yright
end

@testset "bridge complex test" begin
    bridge_complex_test(ℝ⁺, ℝ)
    bridge_complex_test(ℝ, ℝ⁺)
    bridge_complex_test(-1.0..1.0, ℝ⁺)
    bridge_complex_test(ℝ⁺, -1.0..1.0)
    bridge_complex_test(ℝ, 𝕀; mapping = REALCIRCLE)
    bridge_complex_test(ℝ, 𝕀; mapping = REALCIRCLE ∘ Multiply(4.0))
    bridge_complex_test(𝕀, ℝ; mapping = INVREALCIRCLE ∘ Multiply(4.0))
    bridge_complex_test(𝕀, 0..5.0; RR = true)
    bridge_complex_test(-∞..5, -∞..9; RR = true)
    bridge_complex_test(ℝ⁺, ℝ⁻; RR = true)
    bridge_complex_test(ℝ⁻, ℝ⁺; RR = true)
    bridge_complex_test(ℝ, 𝕀)
    bridge_complex_test(𝕀, ℝ)
end
