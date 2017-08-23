module ContinuousTransformations

using ArgCheck
using AutoHashEquals
using MacroTools
using Parameters
using StatsFuns
using Lazy
using Unrolled

export
    AbstractInterval, RealLine, ℝ, PositiveRay, ℝ⁺, NegativeRay, ℝ⁻,
    Segment, width,
    UnivariateTransformation, image, isincreasing, inverse, logjac,
    Affine,
    Negation, NEGATION, Logistic, LOGISTIC, RealCircle, REALCIRCLE, Exp, EXP,
    affine_bridge, default_transformation, transformation_to,
    ArrayTransformation, TransformationTuple
    
import Base:
    in, length, size, ∘, show, getindex,
    middle, linspace, intersect, extrema, minimum, maximum, isfinite, isinf, isapprox
    

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

"""
    _fma(x, y, z)

Placeholder for `Base.fma` until https://github.com/JuliaDiff/ReverseDiff.jl/issues/86 is fixed.
"""
_fma(x, y, z) = x*y+z

######################################################################
# intervals
######################################################################

abstract type AbstractInterval end

isinf(x::AbstractInterval) = !isfinite(x)

isapprox(::AbstractInterval, ::AbstractInterval; rtol=√eps(), atol=0) = false

extrema(x::AbstractInterval) = minimum(x), maximum(x)

"""
    RealLine()

The (extended) real line [-∞,∞]. Use the constant ℝ.
"""
struct RealLine <: AbstractInterval end

isapprox(::RealLine, ::RealLine; rtol=√eps(), atol=0) = true

show(io::IO, ::RealLine) = print(io, "ℝ")

const ℝ = RealLine()

in(x::Real, ::RealLine) = true

minimum(::RealLine) = -Inf
maximum(::RealLine) = Inf

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
minimum(ray::PositiveRay) = ray.left
maximum(::PositiveRay{T}) where T = T(Inf)
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
minimum(::NegativeRay{T}) where T = -T(Inf)
maximum(ray::NegativeRay) = ray.right
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
minimum(s::Segment) = s.left
maximum(s::Segment) = s.right
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
    logjac(t, x)

The log of the determinant of the Jacobian of `t` at `x`.
```
"""
function logjac end

"""
    inverse(t, x)

Return ``t⁻¹(x)``.
"""
function inverse end

"""
    image(transformation)

Return the image of the transformation.
"""
function image end

"""
Continuous bijection ℝⁿ→ℝⁿ.
"""
abstract type ContinuousTransformation <: Function end

"""
    rhs_string(transformation, term)

Return the formula representing the hand side of the `transformation`, with `term` as the argument.
"""
function rhs_string end

transformation_string(t, x = "x") = x * " ↦ " * rhs_string(t, x)

show(io::IO, ::MIME"text/plain", t::ContinuousTransformation) =
    println(io, transformation_string(t, "x"))

show(io::IO, t::ContinuousTransformation) =
    print(io, transformation_string(t, "x"))

######################################################################
# univariate transformations
######################################################################

"""
Univariate monotone transformation, either *increasing* or *decreasing* on the whole domain (thus, a bijection).
"""
abstract type UnivariateTransformation <: ContinuousTransformation end

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
(t::Affine)(x) = _fma(x, t.α, t.β)
logjac(t::Affine, x) = log(abs(t.α))
inverse(t::Affine, x) = _fma(x, 1/t.α, -t.β/t.α)
isincreasing(t::Affine) = true

(t::Affine)(x::Segment) = Segment(t(x.left), t(x.right))
(t::Affine)(x::PositiveRay) = PositiveRay(t(x.left))
(t::Affine)(x::NegativeRay) = NegativeRay(t(x.right))
(t::Affine)(::RealLine) = ℝ

function rhs_string(t::Affine, term)
    @unpack α, β = t
    α == 1 || (term = "$(α)⋅" * term)
    β == 0 || (term = term * " + $(β)")
    term
end

"""
    Negation.

Mapping ``ℝ → ℝ`` using ``x ↦ -x``.
"""
@define_singleton Negation <: UnivariateTransformation

image(::Negation) = ℝ
(::Negation)(x) = -x
logjac(::Negation, x) = zero(x)
inverse(::Negation, x) = -x
isincreasing(::Negation) = false

(::Negation)(x::Segment) = Segment(-x.right, -x.left)
(::Negation)(x::PositiveRay) = NegativeRay(-x.left)
(::Negation)(x::NegativeRay) = PositiveRay(-x.right)
(::Negation)(::RealLine) = ℝ

rhs_string(::Negation, term) = "-" * term

"""
    Logistic()

Mapping ``ℝ → (0,1)`` using ``x ↦ x/(1+x)``.
"""
@define_singleton Logistic <: UnivariateTransformation

image(::Logistic) = Segment(0, 1)
(t::Logistic)(x) = logistic(x)
logjac(::Logistic, x) = -(log1pexp(x)+log1pexp(-x))
inverse(::Logistic, x) = logit(x)
isincreasing(::Logistic) = true

rhs_string(::Logistic, term) = "logistic($term)"

"""
    RealCircle()

Mapping ``ℝ → (-1,1)`` using ``x ↦ x/√(1+x^2)``.
"""
@define_singleton RealCircle <: UnivariateTransformation

image(::RealCircle) = Segment(-1, 1)
(t::RealCircle)(x) = isinf(x) ? sign(x) : x/√(1+x^2)
logjac(::RealCircle, x) = -1.5*log1psq(x)
inverse(::RealCircle, x) = x/√(1-x^2)
isincreasing(::RealCircle) = true

rhs_string(::RealCircle, term) = "realcircle($term)"

"""
    Exp()

Mapping ``ℝ → ℝ⁺`` using ``x ↦ exp(x)``.
"""
@define_singleton Exp <: UnivariateTransformation

image(::Exp) = ℝ⁺
(t::Exp)(x) = exp(x)
logjac(::Exp, x) = x
inverse(::Exp, x) = log(x)
isincreasing(::Exp) = true

rhs_string(::Exp, term) = "exp($term)"

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

rhs_string(t::ComposedTransformation, term) = rhs_string(t.f, rhs_string(t.g, term))

image(t::ComposedTransformation) = t.f(image(t.g))

function logjac(t::ComposedTransformation, x)
    @unpack f, g = t 
    y = g(x)
    log_g′x = logjac(g, x)
    log_f′y = logjac(f, y)
    log_f′y + log_g′x
end


inverse(t::ComposedTransformation, x) = inverse(t.g, inverse(t.f, x))

isincreasing(c::ComposedTransformation) = isincreasing(c.f) == isincreasing(c.g)


∘(f::UnivariateTransformation, g::UnivariateTransformation) =
    ComposedTransformation(f, g)

∘(f::Affine, g::Affine) = Affine(f.α*g.α, _fma(f.α, g.β, f.β))

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
    Affine(1, y.right + x.left) ∘ NEGATION
affine_bridge(x::NegativeRay, y::PositiveRay) =
    Affine(1, y.left + x.right) ∘ NEGATION

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
# array transformations
######################################################################

"""
    ArrayTransformation(transformation, dims)
    ArrayTransformation(transformation, dims...)

Apply transformation to a vector, returning an array of the given dimensions.
"""
struct ArrayTransformation{T <: UnivariateTransformation, D} <: ContinuousTransformation
    transformation::T
    function ArrayTransformation(transformation::T,
                                 dims::Tuple{Vararg{Int64, N}}) where {T,N}
        @argcheck all(dims .> 0) "Invalid dimensions."
        new{T,dims}(transformation)
    end
end

ArrayTransformation(transformation, dims::Int...) =
    ArrayTransformation(transformation, dims)

size(t::ArrayTransformation{T,D}) where {T,D} = D

length(t::ArrayTransformation) = prod(size(t))

function transformation_string(t::ArrayTransformation, term)
    s = size(t)
    s_print = length(s) == 1 ? s[1] : s
    "$(transformation_string(t.transformation, term)) for $(s_print) elements"
end

(t::ArrayTransformation)(x) = reshape((t.transformation).(x), size(t))

function logjac(t::ArrayTransformation, x)
    @argcheck length(x) == length(t) DimensionMismatch
    lj(x) = logjac(t.transformation, x)
    sum(lj, x)
end

inverse(t::ArrayTransformation, x) =
    reshape(inverse.(t.transformation, x), size(t))

image(t::ArrayTransformation) = fill(image(t.transformation), size(t))

@forward ArrayTransformation.t isincreasing

######################################################################
# transformation tuple
######################################################################

struct TransformationTuple{T <: Tuple{Vararg{ContinuousTransformation}}} <:
    ContinuousTransformation
    transformations::T
end

transformation_string(t::TransformationTuple, x) =
    "TransformationTuple" * repr(t.transformations)

length(t::TransformationTuple) = sum(length(t) for t in t.transformations)

image(t::TransformationTuple) = image.(t.transformations)

getindex(t::TransformationTuple, ix) = t.transformations[ix]

function next_indexes(acc, t)
    l = length(t)
    acc + l, acc + (1:l)
end

next_indexes(acc, t::UnivariateTransformation) = acc+1, acc+1

@unroll function transformation_indexes(ts)
    acc = 0
    result = ()
    @unroll for t in ts
        acc, ix = next_indexes(acc, t)
        result = (result..., ix)
    end
    result
end

function (t::TransformationTuple)(x)
    ts = t.transformations
    map((t,ix) -> t(x[ix]), ts, transformation_indexes(ts))
end

function logjac(t::TransformationTuple, x)
    ts = t.transformations
    sum(map((t,ix) -> logjac(t, x[ix]), ts, transformation_indexes(ts)))
end

inverse(t::TransformationTuple, y::Tuple) = vcat(map(inverse, t.transformations, y)...)

end # module
