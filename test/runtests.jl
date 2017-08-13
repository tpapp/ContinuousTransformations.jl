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
    @test PositiveRay(1.0) == PositiveRay(1.0)
    @test NegativeRay(2.0) == NegativeRay(2.0)
    @test RealLine() == ℝ
    @test isa(Segment(1,2.0), Segment{Float64})
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










# """
# Test univariate transformation `f` with `x`. Tests for:

# 1. type of the transformed value,
# 2. whether it is in the range,
# 3. inverse,
# 4. jacobian determinant and its log using automatic differentiation.

# `forwarddiff_exceptions` is a dictionary handling exceptions that ForwardDiff
# cannot cope with at the moment. See eg
# ## workaround for https://github.com/JuliaDiff/ForwardDiff.jl/issues/209
# """
# function test_univariate_scalar{T}(f::UnivariateTransformation, x::T;
#                                    AD_exceptions = Dict())
#     y = f(x)
#     @test y ∈ image(f)
#     @test inv(f)(y) ≈ x
#     (y2, deriv) = f(x, DERIV)
#     expected_deriv = get(AD_exceptions, x, derivative(f, x))
#     @test y == y2
#     @test deriv ≈ expected_deriv
#     (y3, jac) = f(x, JAC)
#     @test y == y3
#     @test jac ≈ abs(expected_deriv)
#     (y3, logjac) = f(x, LOGJAC)
#     @test y == y3
#     @test logjac ≈ log(abs(expected_deriv))
# end

# # some exceptions below
# logit_exceptions(t=Multiply(1.0)) = Dict(t(1.0) => Inf)

# logistic_exceptions() = Dict(-Inf => 0.0)

# """
# Test that univariate transformations map an interval the correct way.
# """
# function test_univariate_interval(f::UnivariateTransformation, x::AbstractInterval)
#     y = f(x)
#     left, right = extrema(x)
#     f_left, f_right = f(left), f(right)
#     if isincreasing(f)
#         @test f_left < f_right
#     else
#         @test f_right < f_left
#         f_right, f_left = f_left, f_right
#     end
#     y_left, y_right = extrema(y)
#     @test y_left ≈ f_left
#     @test y_right ≈ f_right
# end

# """
#     test_univariate(f; N=500, AD_exceptions = Dict())

# Test univariate transformation `f`, called with random numbers and intervals generated by `randform`, which should return numbers in the domain.
# """
# function test_univariate(f; N=500, AD_exceptions = Dict())
#     dom = domain(f)
#     for i in 1:500
#         test_univariate_scalar(f, rand(dom), AD_exceptions = AD_exceptions)
#     end
#     for i in 1:500
#         test_univariate_interval(f, random_interval(dom))
#     end
#     for x in scalars_outside_interval(dom)
#         @test_throws DomainError f(x)
#     end
# end
