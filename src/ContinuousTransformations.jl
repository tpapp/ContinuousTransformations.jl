module ContinuousTransformations

using ArgCheck
using AutoHashEquals
using MacroTools
using Parameters
using StatsFuns

export
    AbstractInterval, RealLine, ℝ, PositiveRay, ℝ⁺, NegativeRay, ℝ⁻,
    Segment, width,
    UnivariateTransformation, image, isincreasing, INV, LOGJAC,
    Affine,
    Negation, NEGATION, Logistic, LOGISTIC, RealCircle, REALCIRCLE, Exp, EXP,
    affine_bridge, default_transformation, transformation_to
    
import Base:
    in, length, size, ∘,
    middle, linspace, intersect, extrema, isfinite, isinf, isapprox
    

######################################################################
# utilities
######################################################################

"""
    define_isapprox(T, fields...)

Define an `isapprox` method, comparing the given fields in type `T`.
"""
macro define_isapprox(T, fields...)
    body = foldl((a, b) -> :($(a) && $(b)),
                 :(isapprox(x.$(field), y.$(field); rtol=rtol, atol=atol))
                 for field in fields)
    quote
        function Base.isapprox(x::$T, y::$T; rtol::Real=√eps(), atol::Real=0)
            $(body)
        end
    end
end

"""
Define a singleton type with the given name and supertype (specified
as `name <: supertype`), and a constant which defaults to the name in
uppercase.
"""
macro define_singleton(name_and_supertype, constant = nothing)
    @capture name_and_supertype name_ <: supertype_
    if constant == nothing
        constant = Symbol(uppercase(string(name)))
    end
    quote
        Core.@__doc__ struct $(esc(name)) <: $(esc(supertype)) end
        Core.@__doc__ const $(esc(constant)) = $(esc(name))()
    end
end


######################################################################
# intervals
######################################################################

abstract type AbstractInterval end

isinf(x::AbstractInterval) = !isfinite(x)

isapprox(::AbstractInterval, ::AbstractInterval; rtol=√eps(), atol=0) = false

"""
    RealLine()

The (extended) real line [-∞,∞]. Use the constant ℝ.
"""
struct RealLine <: AbstractInterval end

isapprox(::RealLine, ::RealLine; rtol=√eps(), atol=0) = true

show(io::IO, ::RealLine) = print(io, "ℝ")

const ℝ = RealLine()

in(x::Real, ::RealLine) = true

extrema(::RealLine) = -Inf, Inf

isfinite(::RealLine) = false

"""
    PositiveRay(left)

The interval `[left, ∞]`.
"""
@auto_hash_equals struct PositiveRay{T <: Real} <: AbstractInterval
    left::T
    function PositiveRay{T}(left::T) where T
        @argcheck isfinite(left) "Need finite endpoint."
        new(left)
    end
end

PositiveRay{T}(left::T) = PositiveRay{T}(left)

in(x::Real, ray::PositiveRay) = ray.left ≤ x

extrema(ray::PositiveRay) = ray.left, Inf

isfinite(::PositiveRay) = false

const ℝ⁺ = PositiveRay(0.0)

@define_isapprox PositiveRay left

"""
    NegativeRay(right)

The interval `[-∞,right]`.
"""
@auto_hash_equals struct NegativeRay{T <: Real} <: AbstractInterval
    right::T
    function NegativeRay{T}(right::T) where T
        @argcheck isfinite(right) "Need finite endpoint."
        new(right)
    end
end

NegativeRay{T}(right::T) = NegativeRay{T}(right)

in(x::Real, ray::NegativeRay) = x ≤ ray.right

extrema(ray::NegativeRay) = -Inf, ray.right

isfinite(::NegativeRay) = false

const ℝ⁻ = NegativeRay(0.0)

@define_isapprox NegativeRay right

"""
    Segment(left, right)

The finite interval `[left, right]`, with ``-∞ < a < b < ∞`` enforced.
"""
@auto_hash_equals struct Segment{T <: Real} <: AbstractInterval
    left::T
    right::T
    function Segment{T}(left::T, right::T) where T
        @argcheck isfinite(left) && isfinite(right) "Need finite endpoints."
        @argcheck left < right "Need strictly increasing endpoints."
        new(left, right)
    end
end

Segment{T <: Real}(left::T, right::T) = Segment{T}(left, right)

Segment(left::Real, right::Real) = Segment(promote(left, right)...)

in(x::Real, s::Segment) = s.left ≤ x ≤ s.right

extrema(s::Segment) = s.left, s.right

isfinite(::Segment) = true

@define_isapprox Segment left right

width(s::Segment) = s.right - s.left

middle(s::Segment) = middle(s.left, s.right)

linspace(s::Segment, n = 50) = linspace(s.left, s.right, n)


######################################################################
# intersections
######################################################################

## general fallback method. define specific methods with the following
## argument ordering Segment < PositiveRay < NegativeRay < RealLine < all
intersect(a::AbstractInterval, b::AbstractInterval) = intersect(b, a)
    
intersect(a::RealLine, b::AbstractInterval) = b

"""
    _maybe_segment(a, b)

Helper function for forming a segment when possible. Internal, not exported.
"""
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

######################################################################
# abstract interface for transformations
######################################################################

"""
The log of the determinant of the Jacobian. Use as
```julia
transformation(x, LOGJAC)
```
"""
@define_singleton LogJac <: Any

"Inverse of the transformation. Use as
```julia
transformation(x, INV)
```"
@define_singleton Inv <: Any

"""
    image(transformation)

Return the image of the transformation.
"""
function image end


######################################################################
# univariate transformations
######################################################################

"""
Univariate monotone transformation, either *increasing* or *decreasing* on the whole domain (thus, a bijection).
"""
abstract type UnivariateTransformation <: Function end

length(::UnivariateTransformation) = 1

size(::UnivariateTransformation) = ()

"""
    isincreasing(transformation)

Return `true` (`false`), when the transformation is monotone increasing (decreasing).
"""
function isincreasing end

"""
    Affine(α, β)

Mapping ``ℝ → ℝ`` using ``x ↦ α⋅x + β``. `α > 0` is enforced, see `Negation`.
"""
@auto_hash_equals struct Affine{T <: Real} <: UnivariateTransformation
    α::T
    β::T
    function Affine{T}(α::T, β::T) where T
        @argcheck α > 0 DomainError()
        new(α, β)
    end
end

Affine(α::T, β::T) where T = Affine{T}(α, β)
Affine(α, β) = Affine(promote(α, β)...)

image(::Affine) = ℝ
(t::Affine)(x) = fma(x, t.α, t.β)
(t::Affine)(x, ::LogJac) = log(abs(t.α))
(t::Affine)(x, ::Inv) = fma(x, 1/t.α, -t.β/t.α)
isincreasing(t::Affine) = true

(t::Affine)(x::Segment) = Segment(t(x.left), t(x.right))
(t::Affine)(x::PositiveRay) = PositiveRay(t(x.left))
(t::Affine)(x::NegativeRay) = NegativeRay(t(x.right))
(t::Affine)(::RealLine) = ℝ

"""
    Negation.

Mapping ``ℝ → ℝ`` using ``x ↦ -x``.
"""
@define_singleton Negation <: UnivariateTransformation

image(::Negation) = ℝ
(t::Negation)(x) = -x
(t::Negation)(x, ::LogJac) = zero(x)
(t::Negation)(x, ::Inv) = -x
isincreasing(::Negation) = false

(::Negation)(x::Segment) = Segment(-x.right, -x.left)
(::Negation)(x::PositiveRay) = NegativeRay(-x.left)
(::Negation)(x::NegativeRay) = PositiveRay(-x.right)
(::Negation)(::RealLine) = ℝ

"""
    Logistic()

Mapping ``ℝ → (0,1)`` using ``x ↦ x/(1+x)``.
"""
@define_singleton Logistic <: UnivariateTransformation

image(::Logistic) = Segment(0, 1)
(t::Logistic)(x) = logistic(x)
(t::Logistic)(x, ::LogJac) = -(log1pexp(x)+log1pexp(-x))
(t::Logistic)(x, ::Inv) = logit(x)
isincreasing(::Logistic) = true

"""
    RealCircle()

Mapping ``ℝ → (-1,1)`` using ``x ↦ x/√(1+x^2)``.
"""
@define_singleton RealCircle <: UnivariateTransformation

image(::RealCircle) = Segment(-1, 1)
(t::RealCircle)(x) = isinf(x) ? sign(x) : x/√(1+x^2)
(t::RealCircle)(x, ::LogJac) = -1.5*log1psq(x)
(t::RealCircle)(x, ::Inv) = x/√(1-x^2)
isincreasing(::RealCircle) = true

"""
    Exp()

Mapping ``ℝ → ℝ⁺`` using ``x ↦ exp(x)``.
"""
@define_singleton Exp <: UnivariateTransformation

image(::Exp) = ℝ⁺
(t::Exp)(x) = exp(x)
(t::Exp)(x, ::LogJac) = x
(t::Exp)(x, ::Inv) = log(x)
isincreasing(::Exp) = true

######################################################################
# composed transformations
######################################################################

"""
    ComposedTransformation(f, g)

Compose two univariate transformations. Use the `∘` operator for construction.
"""
struct ComposedTransformation{Tf <: UnivariateTransformation,
                              Tg <: UnivariateTransformation} <:
                                  UnivariateTransformation
    f::Tf
    g::Tg
end

(c::ComposedTransformation)(x) = c.f(c.g(x))

function show(io::IO, c::ComposedTransformation)
    show(io, c.f)
    print(io, " ∘ ")
    show(io, c.g)
end

function (c::ComposedTransformation)(x, ::LogJac)
    y, log_g′x = c.g(x, LOGJAC)
    fy, log_f′y = c.f(y, LOGJAC)
    fy, log_f′y + log_g′x
end

image(t::ComposedTransformation) = t.f(image(t.g))

isincreasing(c::ComposedTransformation) = isincreasing(c.f) == isincreasing(c.g)

(t::ComposedTransformation)(x, ::Inv) = g(x, f(x, INV), INV)

∘(f::UnivariateTransformation, g::UnivariateTransformation) = ComposedTransformation(f, g)

∘(f::Affine, g::Affine) = Affine(f.α*g.α, fma(f.α, g.β, f.β))

######################################################################
# calculated transformations
######################################################################

function affine_bridge(x::Segment, y::Segment)
    α = width(y) / width(x)
    β = middle(y) - middle(x) * α
    Affine(α, β)
end

affine_bridge(::RealLine, ::RealLine) = Affine(1, 0)

affine_bridge(x::PositiveRay, y::PositiveRay) = Affine(1, y.left - x.left)
affine_bridge(x::NegativeRay, y::NegativeRay) = Affine(1, y.right - x.right)
affine_bridge(x::PositiveRay, y::NegativeRay) =
    Affine(1, y.right + x.left) ∘ Negation
affine_bridge(x::NegativeRay, y::PositiveRay) =
    Affine(1, y.left - x.right) ∘ Negation

default_transformation(::Segment) = REALCIRCLE
default_transformation(::PositiveRay) = EXP
default_transformation(::NegativeRay) = EXP
default_transformation(::RealLine) = Affine(1, 0)

"""
    transformation_to(y, [transformation])

Return a transformation that maps ℝ (or ℝⁿ when applicable) to `y`. The second argument may be used to specify a particular transformation, otherwise `default_transformation` is used.
"""
transformation_to(y::AbstractInterval,
                  transformation = default_transformation(y)) = 
                      affine_bridge(image(transformation), y) ∘ transformation

######################################################################
# vector transformations
######################################################################

struct VectorTransformation{T}
    transformations::T
end

image(t::VectorTransformation) = image.(t.transformations)
(t::VectorTransformation)(x) = map((t,x)->t(x), t.transformations, x)
(t::VectorTransformation)(x, ::LogJac) =
    sum(t(x, LOGJAC) for (t,x) in zip(t.transformations, x))
(t::VectorTransformation)(x, ::Inv) = map((t,x)->t(x, INV), t.transformations, x)
length(t::VectorTransformation) = sum(length(t) for t in t.transformations)
size(t::VectorTransformation) = (length(t), )

end # module
