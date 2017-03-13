export
    AbstractInterval, ∞,
    RealLine, ℝ,
    PositiveRay, ℝ⁺,
    NegativeRay, ℝ⁻,
    Segment, 𝕀, width,
    Interval, ..

abstract AbstractInterval

import Base: in, show, middle, linspace, intersect

"The real line [-∞,∞]."
immutable RealLine <: AbstractInterval
end 

show(io::IO, ::RealLine) = print("-∞..∞")

const ℝ = RealLine()

const ∞ = Inf

in(x::Real, ::RealLine) = true

"""
The interval [left,∞).
"""
@auto_hash_equals immutable PositiveRay{T <: Real} <: AbstractInterval
    left::T
    function PositiveRay(left)
        @argcheck isfinite(left) "Need finite endpoint."
        new(left)
    end
end

show(io::IO, ray::PositiveRay) = println("$(ray.left)..∞")

PositiveRay{T}(left::T) = PositiveRay{T}(left)

const ℝ⁺ = PositiveRay(0.0)

in(x::Real, ray::PositiveRay) = ray.left ≤ x

"""
The interval (-∞,right).
"""
@auto_hash_equals immutable NegativeRay{T <: Real} <: AbstractInterval
    right::T
    function NegativeRay(right)
        @argcheck isfinite(right) "Need finite endpoint."
        new(right)
    end
end

show(io::IO, ray::NegativeRay) = println("-∞..$(ray.right)")

NegativeRay{T}(right::T) = NegativeRay{T}(right)

const ℝ⁻ = NegativeRay(0.0)

in(x::Real, ray::NegativeRay) = x ≤ ray.right

"""
The interval [a,b], with a < b enforced.
"""
@auto_hash_equals immutable Segment{T <: Real} <: AbstractInterval
    left::T
    right::T
    function Segment(left, right)
        @argcheck isfinite(left) && isfinite(right) "Need finite endpoints."
        @argcheck left < right "Need strictly increasing endpoints."
        new(left, right)
    end
end

show(io::IO, s::Segment) = println("$(s.left)..$(s.right)")

Segment{T <: Real}(left::T, right::T) = Segment{T}(left, right)

Segment(left::Real, right::Real) = Segment(promote(left, right)...)

in(x::Real, s::Segment) = s.left ≤ x ≤ s.right

width(s::Segment) = s.right - s.left

middle(s::Segment) = middle(s.left, s.right)

linspace(s::Segment, n = 50) = linspace(s.left, s.right, n)

"Unit interval."
𝕀 = Segment(0.0, 1.0)

"""
Create a RealLine, Segment, or Ray, depending on the arguments.
"""
function Interval(left::Real, right::Real)
    if isfinite(left) && isfinite(right)
        Segment(left, right)
    elseif isfinite(left) && right == Inf
        PositiveRay(left)
    elseif left == -Inf && isfinite(right)
        NegativeRay(right)
    elseif left == -Inf && right == Inf
        ℝ
    else
        throw(ArgumentError("Can't interpret ($left, $right) as an interval."))
    end
end

@inline ..(left, right) = Interval(left, right)

"""
Return the image of an interval for a monotone map (increasing or
decreasing).
"""
function monotone_map_interval(f, x::AbstractInterval, increasing)
    left = isa(x, Union{Segment, PositiveRay}) ? x.left : -Inf
    right = isa(x, Union{Segment, NegativeRay}) ? x.right : Inf
    increasing ? Interval(f(left), f(right)) : Interval(f(right), f(left))
end

######################################################################
# intersections
######################################################################

## general fallback method. define specific methods with the following
## argument ordering Segment < PositiveRay < NegativeRay < RealLine < all
intersect(a::AbstractInterval, b::AbstractInterval) = intersect(b, a)
    
intersect(a::RealLine, b::AbstractInterval) = b

"Helper function for forming a segment when possible."
@inline function _maybe_segment(a, b)
    # NOTE Decided not to represent the empty interval, as it has no use in
    # the context of this package. Best to throw an error as soon as possible.
    a < b ? Segment(a, b) : error("intersection of intervals is empty")
end

intersect(a::Segment, b::Segment) =
    _maybe_segment(max(a.left, b.left), min(a.right, b.right))

intersect(a::Segment, b::PositiveRay) =
    _maybe_segment(max(b.left, a.left), a.right)

intersect(a::Segment, b::NegativeRay) =
    _maybe_segment(a.left, min(a.right, b.right))

intersect(a::PositiveRay, b::PositiveRay) = PositiveRay(max(a.left, b.left))

intersect(a::PositiveRay, b::NegativeRay) = _maybe_segment(a.left, b.right)

intersect(a::NegativeRay, b::NegativeRay) = NegativeRay(min(a.right, b.right))
