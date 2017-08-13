using ContinuousTransformations
using ArgCheck
using Base.Test
import ForwardDiff: derivative

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

# """
#     rand_in(x::AbstractInterval)

# Return a random float in an interval (for testing). Do not use
# directly when testing extrema, as these may happen with `0`
# probability.
# """
# rand_in(::RealLine, ::Type{Val{true}}) = randn()
# rand_in(ray::PositiveRay, ::Type{Val{true}}) = ray.left + abs(randn())
# rand_in(ray::NegativeRay, ::Type{Val{true}}) = ray.right - abs(randn())
# rand_in(seg::Segment, ::Type{Val{true}}) = seg.left + width(seg) * rand()

# function rand(x::AbstractInterval; left_prob = 0.1, right_prob = 0.1)
#     @argcheck left_prob + right_prob ≤ 1
#     z = rand()
#     left, right = extrema(x)
#     if z < left_prob
#         left
#     elseif z > 1-right_prob
#         right
#     else
#         rand(x, Val{true})
#     end
# end


# """
#     random_interval(interval)

# Return a function that generates a random sub-interval of the given `interval`. Useful for unit tests.
# """
# function random_interval(interval::AbstractInterval)
#     a = rand(interval; right_prob = 0)
#     b = rand(interval; left_prob = 0)
#     if a > b
#         a, b = b, a
#     end
#     Interval(a, b)
# end

# """
#     scalars_outside_interval(interval)

# Return a vector of values outside the `interval` (for testing).
# """
# scalars_outside_interval(::RealLine) = []
# scalars_outside_interval(ray::PositiveRay) = ray.left - [0.1, 1.0, 2.0, Inf]
# scalars_outside_interval(ray::NegativeRay) = ray.right + [0.1, 1.0, 2.0, Inf]
# function scalars_outside_interval(seg::Segment)
#     vcat(scalars_outside_interval(PositiveRay(seg.left)),
#          scalars_outside_interval(NegativeRay(seg.right)))
# end

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

######################################################################
# test: univariate transformation
######################################################################

@testset "univariate transformation basics" begin
    @test_throws DomainError Affine(0, 1)
    @test_throws DomainError Affine(-1, 2.0)
    @test Affine(1, 2.0) == Affine(1.0, 2.0)
    a = Affine(1,2)
    @test image(a) == ℝ
    @test isincreasing(a)

    @test image(NEGATION) == ℝ
    @test !isincreasing(NEGATION)
    
    @test image(LOGISTIC) == Segment(0, 1)
    @test isincreasing(LOGISTIC)

    @test image(REALCIRCLE) == Segment(-1, 1)
    @test isincreasing(EXP)

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
    @inferred t(x)
    y = t(x)
    @test y ∈ image(t)
    @test t(y, INV) ≈ x
    logjac = t(x, LOGJAC)
    deriv = get(AD_exceptions, x, derivative(t, x))
    @test logjac ≈ log(abs(deriv))
end

function test_univariate_random(t::UnivariateTransformation; N=500, AD_exceptions = Dict())
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



