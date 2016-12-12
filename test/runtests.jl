using ContinuousTransformations
using Base.Test
import ForwardDiff: derivative

"""
Test transformation `c` with `x`. Tests for:

1. type of the transformed value,
2. range of the transformed value,
3. inverse,
4. log(abs(det(Jacobian))) using automatic differentiation.
"""
function test_transform(c::UnivariateTransformation, x::Real)
    y = transform(c, x)
    @test isa(y, Real)
    @test y ∈ c
    @test invert(c, y) ≈ x
    (y2, lj) = transform_logjac(c, x)
    @test y == y2
    @test lj ≈ log(abs(derivative(x->transform(c,x), x)))
end

@testset "LowerUpperBound" begin
    for i in 1:1000
        lower = randn()
        test_transform(LowerUpperBound(lower, lower+1+abs(randn())), randn())
        test_transform(UNIT_INTERVAL, randn())
    end
end

@testset "LowerBound" begin
    for i in 1:1000
        test_transform(LowerBound(randn()), randn())
        test_transform(POSITIVE_REAL, randn())
    end
end

@testset "UpperBound" begin
    for i in 1:1000
        test_transform(UpperBound(randn()), randn())
    end
end
