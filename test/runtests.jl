using ContinuousTransformations
using Base.Test
import ForwardDiff: derivative
using InferenceUtilities

######################################################################
# test: utilities
######################################################################

"Test singleton type."
ContinuousTransformations.@define_singleton TestSingleton <: Real

@testset "define singleton" begin
    @test TestSingleton <: Real
    @test isa(TESTSINGLETON, TestSingleton)
    @test repr(@doc(TESTSINGLETON)) == repr(@doc(TestSingleton)) ==
        repr(doc"Test singleton type.")
end

######################################################################
# test: intervals
######################################################################

@testset "interval constructors" begin
    @test_throws ArgumentError Segment(NaN, NaN)
    @test_throws ArgumentError Segment(-Inf, Inf)
    @test_throws ArgumentError Segment(2, -1)
    @test_throws ArgumentError Segment(2, 2)
    @test_throws ArgumentError PositiveRay(-Inf)
    @test_throws ArgumentError PositiveRay(NaN)
    @test_throws ArgumentError NegativeRay(-Inf)
    @test_throws ArgumentError NegativeRay(NaN)
end

@testset "interval equality" begin
    @test PositiveRay(1.0) == PositiveRay(1)
    @test NegativeRay(2.0) == NegativeRay(2)
    @test RealLine() == ℝ
    @test isa(Segment(1,2.0), Segment{Float64})
    @test Segment(1, 2.0) == Segment(1, 2) == Segment(1.0, 2.0)
end

@testset "named interval constants" begin
    @test ℝ == RealLine()
    @test ℝ⁺ == PositiveRay(0.0)
    @test ℝ⁻ == NegativeRay(0.0)
end

@testset "interval isapprox" begin
    @test ℝ ≈ ℝ
    @test !(ℝ ≈ Segment(1,2))
    @test !(ℝ ≈ ℝ⁺)
    @test !(ℝ ≈ ℝ⁻)
    @test PositiveRay(1) ≈ PositiveRay(1+eps())
    @test NegativeRay(-2) ≈ NegativeRay(-2+eps())
    @test Segment(1, 2) ≈ Segment(1+eps(), 2+eps())
end

@testset "intervals ∈, width, extrema, finiteness" begin
    seg = Segment(1.0, 2.0)
    posray = PositiveRay(0.0)
    negray = NegativeRay(1.5)
    # numbers in seg
    @test 1.0 ∈ seg
    @test 1.5 ∈ seg
    @test 2.0 ∈ seg
    @test 0.0 ∉ seg
    @test Inf ∉ seg
    @test_throws MethodError "string" ∈ seg
    # methods of seg
    @test width(seg) == 1.0
    @test middle(seg) == 1.5
    @test linspace(seg, 10) == linspace(1.0, 2.0, 10)
    # numbers in posray
    @test 1.0 ∈ posray
    @test Inf ∈ posray
    @test 0 ∈ posray
    @test -1 ∉ posray
    @test -Inf ∉ posray
    @test_throws MethodError "string" ∈ posray
    # numbers in negray
    @test -Inf ∈ negray
    @test 0 ∈ negray
    @test 1.5 ∈ negray
    @test 2 ∉ negray
    @test Inf ∉ negray
    @test_throws MethodError "string" ∈ negray
    # numbers in the real line
    @test 1 ∈ ℝ
    @test Inf ∈ ℝ
    @test -Inf ∈ ℝ
    @test_throws MethodError "string" ∈ ℝ
    # finiteness
    @test isfinite(seg) && !isinf(seg)
    @test !isfinite(posray) && isinf(posray)
    @test !isfinite(negray) && isinf(negray)
    @test !isfinite(posray) && isinf(posray)
    @test !isfinite(ℝ) && isinf(ℝ)
end

@testset "intervals intersections" begin
    seg = Segment(1.0, 2.0)
    posray = PositiveRay(0.0)
    negray = NegativeRay(1.5)
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
    let seg2 = Segment(1.5, 3.0)
        @test seg ∩ seg2 == seg2 ∩ seg == Segment(1.5, 2.0)
    end
    @test seg ∩ posray == posray ∩ seg == seg
    @test seg ∩ negray == negray ∩ seg == Segment(1.0, 1.5)
    @test negray ∩ posray == posray ∩ negray == Segment(0.0, 1.5)
    @test posray ∩ PositiveRay(2) == PositiveRay(2)
    @test posray ∩ PositiveRay(-2) == posray
    @test negray ∩ NegativeRay(-7) == NegativeRay(-7)
    @test negray ∩ NegativeRay(7) == negray
end

@testset "interval printing" begin
    @test repr(ℝ) == "ℝ"
end

######################################################################
# test: univariate transformation
######################################################################

@testset "univariate transformation basics" begin
    @test_throws DomainError Affine(0, 1)
    @test_throws DomainError Affine(-1, 2.0)
    @test Affine(1, 2.0) == Affine(1.0, 2.0)
    a = Affine(1,2)
    @test domain(a) == ℝ
    @test image(a) == ℝ
    @test isincreasing(a)

    @test domain(NEGATION) == ℝ
    @test image(NEGATION) == ℝ
    @test !isincreasing(NEGATION)
    
    @test domain(LOGISTIC) == ℝ
    @test image(LOGISTIC) == Segment(0, 1)
    @test isincreasing(LOGISTIC)

    @test domain(REALCIRCLE) == ℝ
    @test image(REALCIRCLE) == Segment(-1, 1)
    @test isincreasing(REALCIRCLE)

    @test domain(EXP) == ℝ
    @test image(EXP) == ℝ⁺
    @test isincreasing(EXP)
end

"""
    test_univariate_scalar(f, x; [AD_exceptions::Dict])

Test univariate transformation `f` with `x`. Tests for:

1. type stability,

2. transformed value being in the image,

3. inverse,

4. log Jacobian determinant (using automatic differentiation).

`forwarddiff_exceptions` is a dictionary handling exceptions that ForwardDiff cannot cope with at the moment. See [this discussion](https://github.com/JuliaDiff/ForwardDiff.jl/issues/209).
"""
function test_univariate(t::UnivariateTransformation, x; AD_exceptions = Dict())
    @test length(t) == 1
    @test size(t) == ()
    @test @isinferred t(x)
    y = t(x)
    @test y ∈ image(t)
    @test @isinferred inverse(t, y)
    @test inverse(t, y) ≈ x
    @test @isinferred logjac(t, x)
    lj = logjac(t, x)
    deriv = get(AD_exceptions, x, derivative(t, x))
    @test lj ≈ log(abs(deriv)) rtol=1e-5
end

function test_univariate_random(t::UnivariateTransformation; N=500,
                                AD_exceptions = Dict())
    for _ in 1:N
        test_univariate(t, randn(), AD_exceptions = AD_exceptions)
    end
end

const logistic_AD_exceptions = Dict(-Inf => 0.0)

@testset "basic univariate transformations" begin
    test_univariate_random(Affine(1,2))
    test_univariate_random(NEGATION)
    test_univariate_random(LOGISTIC; AD_exceptions = logistic_AD_exceptions)
    test_univariate_random(REALCIRCLE)
    test_univariate_random(EXP)
end

random_segment() = Segment(sort(randn(2))...)

function test_affine_bridge(x, y)
    x1, x2 = extrema(x)
    y1, y2 = extrema(y)
    a = affine_bridge(x, y)
    if !isincreasing(a)
        y2, y1 = y1, y2
    end
    @test a(x1) ≈ y1
    @test a(x2) ≈ y2
    @test a(x) ≈ y
end

@testset "affine bridge" begin
    test_affine_bridge(ℝ, ℝ)
    for _ in 1:100
        test_affine_bridge(random_segment(), random_segment())
        test_affine_bridge(PositiveRay(randn()), PositiveRay(randn()))
        test_affine_bridge(NegativeRay(randn()), NegativeRay(randn()))
        test_affine_bridge(PositiveRay(randn()), NegativeRay(randn()))
        test_affine_bridge(NegativeRay(randn()), PositiveRay(randn()))
    end
end

@testset "interval negation" begin
    @test NEGATION(Segment(1.0, 2.0)) == Segment(-2.0, -1.0)
    @test NEGATION(PositiveRay(9.0)) == NegativeRay(-9.0)
    @test NEGATION(NegativeRay(7.0)) == PositiveRay(-7.0)
    @test NEGATION(ℝ) == ℝ
end

function test_transformation_to(y)
    @test @isinferred transformation_to(y)
    t = transformation_to(y)
    @test image(t) == y
    test_univariate_random(t)
end

@testset "transformations to an interval" begin
    test_transformation_to(Segment(1, 2))
    test_transformation_to(PositiveRay(9.0))
    test_transformation_to(NegativeRay(-7.0))
end

@testset "show" begin
    @test repr(EXP) == "x ↦ exp(x)"
    @test repr(REALCIRCLE) == "x ↦ realcircle(x)"
    @test repr(LOGISTIC) == "x ↦ logistic(x)"
    @test repr(NEGATION) == "x ↦ -x"
    @test repr(Affine(1,0)) == "x ↦ x"
    @test repr(Affine(2,0)) == "x ↦ 2⋅x"
    @test repr(Affine(2,3)) == "x ↦ 2⋅x + 3"
    @test repr(transformation_to(Segment(0,1))) == "x ↦ 0.5⋅realcircle(x) + 0.5"
end

"""
    rand_Inf!(x, [p])

Replace each element of `x` with Inf or -Inf (equal probability), total with IID probability `p`.
"""
function rand_Inf!(x, p = 0.02)
    for i in 1:length(x)
        a = rand()
        a < p && (x[i] = Inf)
        rand() < 0.5 && (x[i] *= -1)
    end
end

function test_array_transformation(t, dims; N = 500)
    at = ArrayTransformation(t, dims)
    @test image(at) == fill(image(t), dims)
    @test length(at) == prod(dims)
    @test size(at) == dims
    @test_throws DimensionMismatch at(ones(dims .+ 1))
    @test_throws DimensionMismatch logjac(at, ones(dims .+ 1))
    for _ in 1:N
        x = randn(dims)
        @test @isinferred logjac(at, x)
        @test logjac(at, x) == sum(logjac.(t, x))
        ## log jacobian may be meaningless at Inf, introduce Inf's after
        rand_Inf!(x)
        y = t.(x)
        @test @isinferred at(x)
        @test at(x) == y
        @test @isinferred inverse(at, y)
        @test inverse(at, y) == inverse.(t, y)
    end
end

@testset "array transformations" begin
    test_array_transformation(transformation_to.(Segment(1,2)), (2,3))
    test_array_transformation(transformation_to.(ℝ), (4,5))
    test_array_transformation(EXP, (3,7,2))
    test_array_transformation(REALCIRCLE, (3,2))
    @test_throws ArgumentError ArrayTransformation(EXP, -1, 2)
    @test_throws MethodError ArrayTransformation(EXP, "a fish")
    @test repr(ArrayTransformation(EXP, 2, 3)) == repr(EXP) * " for (2, 3) elements"
    @test repr(ArrayTransformation(EXP, 2)) == repr(EXP) * " for 2 elements"
end

@testset "transformation tuple univariate" begin
    ts = transformation_to.((PositiveRay(1.0), NegativeRay(1.0), NegativeRay(1.0),
                             ℝ, Segment(0.0,1.0)))
    tt = TransformationTuple(ts)
    @test length(tt) == sum(length, ts)
    @test image(tt) == image.(ts)
    @test repr(tt) == "TransformationTuple" * repr(ts)
    for i in 1:length(ts)
        @test tt[i] == ts[i]
    end
    x = randn(length(tt))
    @test @isinferred tt(x)
    y = tt(x)
    @test y == map((t,x) -> t(x), ts, tuple(x...))
    @test @isinferred logjac(tt, x)
    @test logjac(tt, x) == sum(map(logjac, ts, tuple(x...)))
    @test @isinferred inverse(tt, y)
    @test inverse(tt, y) == [map(inverse, ts, y)...]
end

@testset "transformation tuple mixed" begin
    ts = (EXP, ArrayTransformation(REALCIRCLE, 2))
    tt = TransformationTuple(ts)
    @test length(tt) == sum(length, ts) == 3
    @test image(tt) == image.(ts)
    @test repr(tt) == "TransformationTuple" * repr(ts)
    for i in 1:length(ts)
        @test tt[i] == ts[i]
    end
    x = randn(length(tt))
    @test @isinferred tt(x)
    y = tt(x)
    @test y == (ts[1](x[1]), ts[2](x[2:3]))
    @test @isinferred logjac(tt, x)
    @test logjac(tt, x) == logjac(ts[1], x[1]) + logjac(ts[2], x[2:3])
    @test @isinferred inverse(tt, y)
    @test inverse(tt, y) == vcat(inverse.(ts, y)...)
end

@testset "transformation tuple inference" begin
    t = TransformationTuple((transformation_to(Segment(0.0,10.0)),
                             ArrayTransformation(Affine(1,0), 2)))
    @test @isinferred t(ones(3))
    @test @isinferred logjac(t, ones(3))
    @test @isinferred inverse(t, (1.0, ones(2)))
end

@testset "transformation tuple by row" begin
    t = TransformationTuple((transformation_to(Segment(0.0,10.0)),
                             ArrayTransformation(Affine(1,0), 2)))
    x = randn(length(t))
    tx = t(x)
    N = 100
    y = repeat(x', inner = (N, 1))
    ty = map_by_row(t, y)
    @test length(ty) == length(tx)
    @test all(size.(ty, 1) .== N)
    @test ty[1] == fill(tx[1], N)
    @test ty[2] == repeat(tx[2]', inner = (N, 1))
end
