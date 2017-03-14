@testset "interval constructors" begin
    @test PositiveRay(1.0) == 1.0..∞
    @test NegativeRay(2.0) == -∞..2.0
    @test RealLine() == ℝ == -∞..∞
end

@testset "interval printing" begin
    @test sprint(show, 1.0..2.0) == "(1.0..2.0)"
    @test sprint(show, 1.0..∞) == "(1.0..∞)"
    @test sprint(show, -∞..2.0) == "(-∞..2.0)"
    @test sprint(show, ℝ) == "ℝ"
end

@testset "intervals basics" begin
    seg = 1.0..2.0
    posray = 0.0..∞
    negray = -∞..1.5
    # numbers in seg
    @test 1.0 ∈ seg
    @test 1.5 ∈ seg
    @test 2.0 ∈ seg
    @test 0.0 ∉ seg
    @test ∞ ∉ seg
    @test_throws MethodError "string" ∈ seg
    # methods of seg
    @test width(seg) == 1.0
    @test middle(seg) == 1.5
    @test linspace(seg, 10) == linspace(1.0, 2.0, 10)
    # numbers in posray
    @test 1.0 ∈ posray
    @test ∞ ∈ posray
    @test 0 ∈ posray
    @test -1 ∉ posray
    @test -∞ ∉ posray
    @test_throws MethodError "string" ∈ posray
    # numbers in negray
    @test -∞ ∈ negray
    @test 0 ∈ negray
    @test 1.5 ∈ negray
    @test 2 ∉ negray
    @test ∞ ∉ negray
    @test_throws MethodError "string" ∈ negray
    # numbers in the real line
    @test 1 ∈ ℝ
    @test ∞ ∈ ℝ
    @test -∞ ∈ ℝ
    @test_throws MethodError "string" ∈ ℝ
    # special intervals
    @test 𝕀== 0.0..1.0
    @test ℝ⁺ == 0.0..∞
    @test ℝ⁻ == -∞..0.0
end

@testset "intervals intersections" begin
    seg = 1.0..2.0
    posray = 0.0..∞
    negray = -∞..1.5
    # intersections with ℝ
    @test seg ∩ ℝ == seg
    @test ℝ ∩ seg == seg
    @test posray ∩ ℝ == posray
    @test ℝ ∩ posray == posray
    @test negray ∩ ℝ == negray
    @test ℝ ∩ negray == negray
    @test ℝ ∩ ℝ == ℝ
    # empty intersections
    @test_throws Exception 𝕀∩ seg
    @test_throws Exception posray ∩ ℝ⁻
    @test_throws Exception ℝ⁻ ∩ seg
    # non-empty intersections
    let seg2 = 1.5..3.0
        @test seg ∩ seg2 == seg2 ∩ seg == 1.5..2.00
    end
    @test seg ∩ posray == posray ∩ seg == seg
    @test seg ∩ negray == negray ∩ seg == 1.0..1.5
    @test negray ∩ posray == posray ∩ negray == 0.0..1.5
end
