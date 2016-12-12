module ContinuousTransformations

using StatsFuns
import Base: in

export
    transform,
    transform_logjac,
    invert,
    UnivariateTransformation,
    LowerUpperBound,
    LowerBound,
    UpperBound

######################################################################
# utilities
######################################################################

"""
Transform a value (in ℝ or ℝⁿ).
"""
function transform end

"""
Transform a value (in ℝ or ℝⁿ), and return a tuple of the result and the 
log of the determinant of the Jacobian.
"""
function transform_logjac end

"""
Invert a transformation. `invert(c, transform(c, x)) == x` should hold.
"""
function invert end

abstract UnivariateTransformation

"""
Transform ℝ to the finite interval (lower,upper), using the logistic function.
"""
immutable LowerUpperBound{T} <: UnivariateTransformation
    lower::T
    upper::T
    function LowerUpperBound(lower, upper)
        (isfinite(lower) && isfinite(upper)) ||
            error(ArgumentError("Bounds need to be finite."))
        lower < upper ||
            error(ArgumentError("Bounds need to be in the right order"))
        new(lower, upper)
    end
end

width(lu::LowerUpperBound) = lu.upper-lu.lower

LowerUpperBound{T <: Real}(lower::T, upper::T) = LowerUpperBound{T}(lower, upper)

LowerUpperBound(lower::Real, upper::Real) =
    LowerUpperBound(promote(lower, upper)...)

in(x::Real, lu::LowerUpperBound) = lu.lower ≤ x ≤ lu.upper

transform(lu::LowerUpperBound, x::Real) = lu.lower+logistic(x)*width(lu)

transform_logjac(lu::LowerUpperBound, x::Real) =
    (transform(lu, x), -2*log1pexp(-x)-x+log(width(lu)))

invert(lu::LowerUpperBound, y) = logit((y-lu.lower)/width(lu))

"""
Transform ℝ to the interval (lower,∞), using an exponential transformation.
"""
immutable LowerBound{T} <: UnivariateTransformation
    lower::T
    function LowerBound(lower)
        isfinite(lower) || error(ArgumentError("Bound needs to be finite."))
        new(lower)
    end
end

LowerBound{T}(lower::T) = LowerBound{T}(lower)

in(x::Real, l::LowerBound) = l.lower ≤ x

transform(l::LowerBound, x::Real) = exp(x)+l.lower

transform_logjac(l::LowerBound, x::Real) = (transform(l, x), x)

invert(l::LowerBound, y::Real) = log(y-l.lower)

"""
Transform ℝ to the interval (-∞, upper), using an exponential transformation.
"""
immutable UpperBound{T} <: UnivariateTransformation
    upper::T
    function UpperBound(upper)
        isfinite(upper) || error(ArgumentError("Bound needs to be finite."))
        new(upper)
    end
end

UpperBound{T}(upper::T) = UpperBound{T}(upper)

in(x::Real, u::UpperBound) = x ≤ u.upper

transform(u::UpperBound, x::Real) = u.upper-exp(x)

transform_logjac(u::UpperBound, x::Real) = (transform(u, x), x)

invert(u::UpperBound, y::Real) = log(u.upper-y)

end # module
