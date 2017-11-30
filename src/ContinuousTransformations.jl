module ContinuousTransformations

using AutoHashEquals
using DocStringExtensions
using MacroTools
using Parameters
using StatsFuns
using Lazy
using Unrolled

export
    AbstractInterval, RealLine, ℝ, PositiveRay, ℝ⁺, NegativeRay, ℝ⁻,
    Segment, width,
    UnivariateTransformation, domain, image, isincreasing, inverse, logjac,
    Affine, IDENTITY,
    Negation, NEGATION, Logistic, LOGISTIC, RealCircle, REALCIRCLE, Exp, EXP,
    Logit, LOGIT, InvRealCircle, INVREALCIRCLE, Log, LOG,
    affine_bridge, default_transformation, bridge,
    ArrayTransformation, TransformationTuple, map_by_row,
    TransformLogLikelihood, get_transformation

import Base:
    in, length, size, ∘, show, getindex, middle, linspace, intersect, extrema,
    minimum, maximum, isfinite, isinf, isapprox

# disabled until https://github.com/JuliaDiff/ForwardDiff.jl/issues/261 is fixed
# const log1p = Base.Math.JuliaLibm.log1p # should be faster and more accurate?



# utilities


"""
    $SIGNATURES

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

Placeholder for `Base.fma` until
https://github.com/JuliaDiff/ReverseDiff.jl/issues/86 is fixed.
"""
_fma(x, y, z) = x*y+z



# abstract interface for transformations


"""
    logjac(t, x)

The log of the determinant of the Jacobian of `t` at `x`.
```
"""
function logjac end

"""
    inverse(t, x)

Return ``t⁻¹(x)``.

    inverse(t)

Return the transformation ``t⁻¹``.
"""
inverse(t, x) = inverse(t)(x)

"""
    domain(transformation)

Return the domain of the transformation.
"""
function domain end

"""
    image(transformation)

Return the image of the transformation.
"""
function image end

"""
Continuous bijection ``D ⊂ ℝⁿ→ I ⊂ ℝⁿ`` or ``D ⊂ ℝ → I ⊂ ℝ``.
"""
abstract type ContinuousTransformation <: Function end

"""
    rhs_string(transformation, term)

Return the formula representing the hand side of the `transformation`, with
`term` as the argument.
"""
function rhs_string end

transformation_string(t, x = "x") = x * " ↦ " * rhs_string(t, x)

show(io::IO, ::MIME"text/plain", t::ContinuousTransformation) =
    println(io, transformation_string(t, "x"))

show(io::IO, t::ContinuousTransformation) =
    print(io, transformation_string(t, "x"))

"""
    bridge(dom, img, [transformation])

Return a transformation that maps `dom` to `img`.

The `transformation` argument may be used to specify a particular transformation
family, otherwise `default_transformation` is used.
"""
function bridge end



# intervals

"""
Abstract supertype for all univariate intervals. It is not specified whether
they are open or closed.
"""
abstract type AbstractInterval end

isinf(x::AbstractInterval) = !isfinite(x)

isapprox(::AbstractInterval, ::AbstractInterval; rtol=√eps(), atol=0) = false

extrema(x::AbstractInterval) = minimum(x), maximum(x)

"""
    RealLine()

The real line. Use the constant [`ℝ`](@ref).
"""
struct RealLine <: AbstractInterval end

isapprox(::RealLine, ::RealLine; rtol=√eps(), atol=0) = true

show(io::IO, ::RealLine) = print(io, "ℝ")

"A constant for the real line."
const ℝ = RealLine()

in(x::Real, ::RealLine) = true

minimum(::RealLine) = -Inf
maximum(::RealLine) = Inf

isfinite(::RealLine) = false

"""
    PositiveRay(left)

The real numbers above `left`. See [`ℝ⁺`](@ref).
"""
@auto_hash_equals struct PositiveRay{T <: Real} <: AbstractInterval
    left::T
    function PositiveRay{T}(left::T) where T
        isfinite(left) || throw(ArgumentError("Need finite endpoint."))
        new(left)
    end
end

PositiveRay{T}(left::T) = PositiveRay{T}(left)

in(x::Real, ray::PositiveRay) = ray.left ≤ x
minimum(ray::PositiveRay) = ray.left
maximum(::PositiveRay{T}) where T = T(Inf)
isfinite(::PositiveRay) = false

"The positive real numbers."
const ℝ⁺ = PositiveRay(0.0)

@define_isapprox PositiveRay left

"""
    NegativeRay(right)

The real numbers below `right`. See [`ℝ⁻`](@ref).
"""
@auto_hash_equals struct NegativeRay{T <: Real} <: AbstractInterval
    right::T
    function NegativeRay{T}(right::T) where T
        isfinite(right) || throw(ArgumentError("Need finite endpoint."))
        new(right)
    end
end

NegativeRay{T}(right::T) = NegativeRay{T}(right)

in(x::Real, ray::NegativeRay) = x ≤ ray.right
minimum(::NegativeRay{T}) where T = -T(Inf)
maximum(ray::NegativeRay) = ray.right
isfinite(::NegativeRay) = false

"The negative real numbers."
const ℝ⁻ = NegativeRay(0.0)

@define_isapprox NegativeRay right

"""
    Segment(left, right)

The real numbers between `left` and `right`, with
``-∞ < \\text{left} < \\text{right} < ∞`` enforced.
"""
@auto_hash_equals struct Segment{T <: Real} <: AbstractInterval
    left::T
    right::T
    function Segment{T}(left::T, right::T) where T
        isfinite(left) && isfinite(right) ||
            throw(ArgumentError("Need finite endpoints."))
        left < right ||
            throw(ArgumentError("Need strictly increasing endpoints."))
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

"""
    $SIGNATURES

Width of a finite interval.
"""
width(s::Segment) = s.right - s.left

middle(s::Segment) = middle(s.left, s.right)

linspace(s::Segment, n = 50) = linspace(s.left, s.right, n)



# intersections

## general fallback method. define specific methods with the following
## argument ordering Segment < PositiveRay < NegativeRay < RealLine < all
intersect(a::AbstractInterval, b::AbstractInterval) = intersect(b, a)

intersect(a::RealLine, b::AbstractInterval) = b

"""
    $SIGNATURES

Helper function for forming a segment when possible. Internal, not exported.
"""
@inline function _maybe_segment(a, b)
    # NOTE Decided not to represent the empty interval, as it has no use in
    # the context of this package. Best to throw an error as soon as possible.
    a < b ? Segment(a, b) :
        throw(ArgumentError("intersection of intervals is empty"))
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
# univariate transformations
######################################################################

"""
Univariate monotone transformation, either *increasing* or *decreasing* on the
whole domain (thus, a bijection).
"""
abstract type UnivariateTransformation <: ContinuousTransformation end

length(::UnivariateTransformation) = 1

size(::UnivariateTransformation) = ()

"""
    isincreasing(transformation)

Return `true` (`false`), when the transformation is monotonically increasing
(decreasing).
"""
function isincreasing end

"""
Trait that is useful for domain and image calculations. See [`RRStable`](@ref).
"""
abstract type RRStability end

"""
Trait that indicates that a univariate transformation

1. maps ``ℝ`` to ``ℝ``,

2. supports mapping intervals, and

3. maps subtypes of [`AbstractInterval`](@ref) to the same type.
"""
struct RRStable <: RRStability end

"""
Trait that indicates that a univariate transformation is *not* [`RRStable`](@ref).
"""
struct NotRRStable <: RRStability end

"""
    $SIGNATURES

Return either the trait [`RRStable`](@ref) and `NotRRStable`.
"""
RR_stability(::UnivariateTransformation) = NotRRStable()

∘(::RRStability, ::RRStability) = NotRRStable()

∘(::RRStable, ::RRStable) = RRStable()

######################################################################
# Affine and Negation
#
# These play a special role since map subtypes of `AbstractInterval` to
# themselves, and thus are useful for composition.
######################################################################

"""
    Affine(α, β)

Mapping ``ℝ → ℝ`` using ``x ↦ α⋅x + β``.

``α > 0`` is enforced, see [`Negation`](@ref).
"""
@auto_hash_equals struct Affine{T <: Real} <: UnivariateTransformation
    α::T
    β::T
    function Affine{T}(α::T, β::T) where T
        α > 0 || throw(DomainError())
        new(α, β)
    end
end

Affine(α::T, β::T) where T = Affine{T}(α, β)
Affine(α, β) = Affine(promote(α, β)...)

"Identity (as an affine transformation)."
const IDENTITY = Affine(1, 0)

domain(::Affine) = ℝ
image(::Affine) = ℝ
(t::Affine)(x) = _fma(x, t.α, t.β)
logjac(t::Affine, x) = log(abs(t.α))
inverse(t::Affine) = Affine(1/t.α, -t.β/t.α)
isincreasing(t::Affine) = true

RR_stability(::Affine) = RRStable()

(t::Affine)(x::Segment) = Segment(t(x.left), t(x.right))
(t::Affine)(x::PositiveRay) = PositiveRay(t(x.left))
(t::Affine)(x::NegativeRay) = NegativeRay(t(x.right))
(t::Affine)(::RealLine) = ℝ

function rhs_string(t::Affine, term)
    @unpack α, β = t
    α == 1 || (term = "$(signif(α, 4))⋅" * term)
    β == 0 || (term = term * " + $(signif(β, 4))")
    term
end

"""
    Negation()

Mapping ``ℝ → ℝ`` using ``x ↦ -x``.
"""
@define_singleton Negation <: UnivariateTransformation

domain(::Negation) = ℝ
image(::Negation) = ℝ
(::Negation)(x) = -x
logjac(::Negation, x) = zero(x)
inverse(::Negation) = NEGATION
isincreasing(::Negation) = false

RR_stability(::Negation) = RRStable()

(::Negation)(x::Segment) = Segment(-x.right, -x.left)
(::Negation)(x::PositiveRay) = NegativeRay(-x.left)
(::Negation)(x::NegativeRay) = PositiveRay(-x.right)
(::Negation)(::RealLine) = ℝ

rhs_string(::Negation, term) = "-" * term

######################################################################
# transformations from ℝ to subsets
######################################################################

"""
    Logistic()

Mapping ``ℝ → (0,1)`` using ``x ↦ 1/(1+\\exp(-x))``.
"""
@define_singleton Logistic <: UnivariateTransformation

domain(::Logistic) = ℝ
image(::Logistic) = Segment(0, 1)
(t::Logistic)(x) = logistic(x)
logjac(::Logistic, x) = -(log1pexp(x)+log1pexp(-x))
inverse(::Logistic) = LOGIT
isincreasing(::Logistic) = true

rhs_string(::Logistic, term) = "logistic($term)"

"""
    RealCircle()

Mapping ``ℝ → (-1,1)`` using ``x ↦ x/√(1+x^2)``.
"""
@define_singleton RealCircle <: UnivariateTransformation

domain(::RealCircle) = ℝ
image(::RealCircle) = Segment(-1, 1)
(t::RealCircle)(x) = isinf(x) ? sign(x) : x/√(1+x^2)
logjac(::RealCircle, x) = -1.5*log1psq(x)
inverse(::RealCircle) = INVREALCIRCLE # x/√(1-x^2)
isincreasing(::RealCircle) = true

rhs_string(::RealCircle, term) = "realcircle($term)"

"""
    Exp()

Mapping ``ℝ → ℝ⁺`` using ``x ↦ \\exp(x)``.
"""
@define_singleton Exp <: UnivariateTransformation

domain(::Exp) = ℝ
image(::Exp) = ℝ⁺
(t::Exp)(x) = exp(x)
logjac(::Exp, x) = x
inverse(::Exp) = LOG
isincreasing(::Exp) = true

rhs_string(::Exp, term) = "exp($term)"

######################################################################
# transformations to ℝ from subsets
######################################################################

"""
    Logit()

Mapping ``(0,1) → ℝ`` using ``x ↦ \\log(x/(1-x))``.
"""
@define_singleton Logit <: UnivariateTransformation

domain(::Logit) = Segment(0, 1)
image(::Logit) = ℝ
(t::Logit)(x) = logit(x)
logjac(::Logit, x) = -log(x)-log1p(-x)
inverse(::Logit) = LOGISTIC
isincreasing(::Logit) = true

rhs_string(::Logit, term) = "logit($term)"

"""
    InvRealCircle()

Mapping ``(-1,1) → ℝ`` using ``x ↦ x/√(1-x^2)``.
"""
@define_singleton InvRealCircle <: UnivariateTransformation

domain(::InvRealCircle) = Segment(-1, 1)
image(::InvRealCircle) = ℝ
(t::InvRealCircle)(x) = x/√(1-x^2)
logjac(::InvRealCircle, x) = -1.5*log1p(-abs2(x))
inverse(::InvRealCircle) = REALCIRCLE
isincreasing(::InvRealCircle) = true

rhs_string(::InvRealCircle, term) = "realcircle⁻¹($term)"

"""
    Log()

Mapping ``ℝ → ℝ⁺`` using ``x ↦ \\exp(x)``.
"""
@define_singleton Log <: UnivariateTransformation

domain(::Log) = ℝ⁺
image(::Log) = ℝ
(t::Log)(x) = log(x)
logjac(::Log, x) = -log(x)
inverse(::Log) = EXP
isincreasing(::Log) = true

rhs_string(::Log, term) = "log($term)"

######################################################################
# composed transformations
######################################################################

"""
    ComposedTransformation(f, g)

Compose two univariate transformations, resulting in the mapping
``f∘g``, or ``x ↦ f(g(x))`.

Use the `∘` operator for construction.
"""
struct ComposedTransformation{Tf <: UnivariateTransformation,
                              Tg <: UnivariateTransformation} <:
                                  UnivariateTransformation
    f::Tf
    g::Tg
end

RR_stability(t::ComposedTransformation) = RR_stability(t.f) ∘ RR_stability(t.g)

(c::ComposedTransformation)(x) = c.f(c.g(x))

rhs_string(t::ComposedTransformation, term) =
    rhs_string(t.f, rhs_string(t.g, term))

# composed image and domain uses traits, effectively allowing calculation for
# composing with Affine and Negation

"""
    composed_image(f_RR_stability, g_RR_stability, f, g)
"""
composed_image(::RRStable, ::RRStability, f, g) = f(image(g))

composed_image(::RRStability, ::RRStable, f, g) = image(f)

composed_image(::RRStable, ::RRStable, f, g) = ℝ

"""
    composed_domain(f_RR_stability, g_RR_stability, f, g)
"""
composed_domain(::RRStable, ::RRStability, f, g) = domain(g)

composed_domain(::RRStability, ::RRStable, f, g) = inverse(g)(domain(f))

composed_domain(::RRStable, ::RRStable, f, g) = ℝ

image(t::ComposedTransformation) =
    composed_image(RR_stability(t.f), RR_stability(t.g), t.f, t.g)

domain(t::ComposedTransformation) =
    composed_domain(RR_stability(t.f), RR_stability(t.g), t.f, t.g)

function logjac(t::ComposedTransformation, x)
    @unpack f, g = t
    y = g(x)
    log_g′x = logjac(g, x)
    log_f′y = logjac(f, y)
    log_f′y + log_g′x
end

inverse(t::ComposedTransformation, x) = inverse(t.g, inverse(t.f, x))

inverse(t::ComposedTransformation) =
    ComposedTransformation(inverse(t.g), inverse(t.f))

isincreasing(c::ComposedTransformation) = isincreasing(c.f) == isincreasing(c.g)

∘(f::UnivariateTransformation, g::UnivariateTransformation) =
    ComposedTransformation(f, g)

∘(f::Affine, g::Affine) = Affine(f.α*g.α, _fma(f.α, g.β, f.β))

######################################################################
# calculated transformations
######################################################################

"""
    affine_bridge(interval1, interval1)

Return an affine transformation between two intervals of the same type.
"""
function affine_bridge(x::Segment, y::Segment)
    α = width(y) / width(x)
    β = middle(y) - middle(x) * α
    Affine(α, β)
end

affine_bridge(::RealLine, ::RealLine) = IDENTITY

affine_bridge(x::PositiveRay, y::PositiveRay) = Affine(1, y.left - x.left)
affine_bridge(x::NegativeRay, y::NegativeRay) = Affine(1, y.right - x.right)
affine_bridge(x::PositiveRay, y::NegativeRay) =
    Affine(1, y.right + x.left) ∘ NEGATION
affine_bridge(x::NegativeRay, y::PositiveRay) =
    Affine(1, y.left + x.right) ∘ NEGATION

"""
    default_transformation(dom, img)

Return a transformation from `dom` that can be mapped to `img` using
`affine_bridge`.
"""
default_transformation(::RealLine, ::Segment) = REALCIRCLE
default_transformation(::RealLine, ::PositiveRay) = EXP
default_transformation(::RealLine, ::NegativeRay) = EXP
default_transformation(::RealLine, ::RealLine) = IDENTITY

default_transformation(::Segment, ::RealLine) = INVREALCIRCLE
default_transformation(::PositiveRay, ::RealLine) = LOG
default_transformation(::NegativeRay, ::RealLine) = LOG

bridge(dom::RealLine, img, t = default_transformation(dom, img)) =
    affine_bridge(image(t), img) ∘ t

bridge(dom, img::RealLine, t = default_transformation(dom, img)) =
    t ∘ affine_bridge(dom, domain(t))

bridge(::RealLine, ::RealLine) = IDENTITY

######################################################################
# array transformations
######################################################################

"""
    ArrayTransformation(transformation, dims)
    ArrayTransformation(transformation, dims...)

Apply transformation to a vector, returning an array of the given dimensions.
"""
struct ArrayTransformation{T <: UnivariateTransformation,
                           D} <: ContinuousTransformation
    transformation::T
    function ArrayTransformation(transformation::T,
                                 dims::Tuple{Vararg{Int64, N}}) where {T,N}
        all(dims .> 0) || throw(ArgumentError("Invalid dimensions."))
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
    length(x) == length(t) || throw(DimensionMismatch())
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

"""
    TransformationTuple(tuple_of_transformations)
    TransformationTuple(transformations...)

A tuple of transformations. Given a vector of matching length, each takes as
many reals as needed, and returns the result as a tuple.
"""
@auto_hash_equals struct
    TransformationTuple{T <: Tuple{Vararg{ContinuousTransformation}}} <:
        ContinuousTransformation
    transformations::T
end

TransformationTuple(ts::ContinuousTransformation...) = TransformationTuple(ts)

function transformation_string(t::TransformationTuple, x)
    s = IOBuffer()
    print(s, "TransformationTuple")
    for (tt, ix) in zip(t.transformations,
                        transformation_indexes(t.transformations))
        print(s, "\n    ")
        print(s, transformation_string(tt, "$(x)[$(ix)]"))
    end
    String(take!(s))
end

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

inverse(t::TransformationTuple, y::Tuple) =
    vcat(map(inverse, t.transformations, y)...)

"""
    map_by_row(t, x)

Apply `t` to the rows of `x`, returning the results as a tuple of vectors or
matrices.

Useful for transforming (posterior) draws from a sample in ℝⁿ.
"""
function map_by_row(t::TransformationTuple, x::AbstractMatrix)
    ts = t.transformations
    ixs = transformation_indexes(ts)
    _t(t, ix::Int) = t.(@view(x[:, ix]))
    _t(t, ix::Range) = mapslices(t, @view(x[:, ix]), 2)
    map(_t, ts, ixs)
end

######################################################################
# loglikelihood wrapper
######################################################################

"""
    TransformLogLikelihood(ℓ, t::Union{Tuple, TransformationTuple})
    TransformLogLikelihood(ℓ, transformations...)

Return a callable that

1. transforms its vector argument using a transformation tuple (or a tuple of
transformations, converted as required) to a tuple of values,

2. calls `ℓ` with these, which should return a scalar,

3. returns the result above corrected by the log Jacobians.

Useful when `ℓ` is a log-likelihood function with a restricted domain, and `t`
is used to trasform to this domain from ``ℝⁿ``.
"""
struct TransformLogLikelihood{L, T <: TransformationTuple}
    loglikelihood::L
    transformation::T
end

TransformLogLikelihood(L, T::Tuple) =
    TransformLogLikelihood(L, TransformationTuple(T))

TransformLogLikelihood(L, ts::ContinuousTransformation...) =
    TransformLogLikelihood(L, TransformationTuple(ts))

@forward TransformLogLikelihood.transformation length

(f::TransformLogLikelihood)(x) =
    f.loglikelihood(f.transformation(x)...) + logjac(f.transformation, x)

get_transformation(d::TransformLogLikelihood) = d.transformation

function show(io::IO, f::TransformLogLikelihood)
    print(io, "TransformLogLikelihood of length $(length(f)), with ")
    print(io, f.transformation)
end

end
