using ContinuousTransformations
using ValidatedNumerics
using Base.Test
import ForwardDiff: derivative
using Cubature

"""
Test univariate transformation `f` with `x`. Tests for:

1. type of the transformed value,
2. whether it is in the range,
3. inverse,
4. jacobian determinant and its log using automatic differentiation.
"""
function test_univariate(f::UnivariateTransformation, x::Real)
    y = f(x)
    @test isa(y, typeof(x))
    @test y ∈ f(domain(f))
    @test inv(f)(y) ≈ x
    (y2, jac) = f(x, JAC)
    @test y == y2
    @test jac ≈ abs(derivative(f, x))
    (y3, logjac) = f(x, LOGJAC)
    @test y == y3
    @test logjac ≈ log(jac)
end

@testset "Logit" begin
    for i in 1:1000 test_univariate(Logit(), rand()) end
end

@testset "Logistic" begin
    for i in 1:1000 test_univariate(Logistic(), rand()) end
end

@testset "Log" begin
    for i in 1:1000 test_univariate(Log(), abs(randn())) end
end

@testset "Exp" begin
    for i in 1:1000 test_univariate(Exp(), randn()) end
end

@testset "OddsRatio" begin
    for i in 1:1000 test_univariate(OddsRatio(), rand()) end
end

@testset "InvOddsRatio" begin
    for i in 1:1000 test_univariate(InvOddsRatio(), abs(randn())) end
end

@testset "Power" begin
    for i in 1:1000 test_univariate(Power(abs(randn())), abs(randn())) end
end

@testset "UpperBound" begin
    for i in 1:1000 test_univariate(Affine(randn(), randn()), randn()) end
end

@testset "Integral" begin
    f, dom = integral_substitution(InvOddsRatio(), x->exp(-x^2), 0..Inf)
    @test hquadrature(f, dom.lo, dom.hi)[1] ≈ √π/2
end
