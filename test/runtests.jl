using ContinuousTransformations
using Base.Test
import ForwardDiff: derivative

function test_transform(c::UnivariateTransformation, x::Real)
    y = transform(c, x)
    @test isa(y, Real)
    @test y ∈ c
    @test invert(c, y) ≈ x
    (y2, lj) = transform_logjac(c, x)
    @test y == y2
    @test lj ≈ log(abs(derivative(x->transform(c,x), x)))
end

for i in 1:1000
    lower = randn()
    test_transform(LowerUpperBound(lower, lower+1+abs(randn())), randn())
end

for i in 1:1000
    test_transform(LowerBound(randn()), randn())
end

for i in 1:1000
    test_transform(UpperBound(randn()), randn())
end
