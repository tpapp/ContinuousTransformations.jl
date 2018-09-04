module ContinuousTransformations

import Base.rand

using ArgCheck
using AutoHashEquals
import Distributions:
    ContinuousUnivariateDistribution, ContinuousMultivariateDistribution, logpdf
using DocStringExtensions
using LinearAlgebra
using MacroTools
using Parameters
using StatsFuns
using Lazy
# using Unrolled

export
    # intervals
    AbstractInterval, RealLine, ℝ, PositiveRay, ℝ⁺, NegativeRay, ℝ⁻,
    Segment, width,
    # general API
    domain, image, isincreasing, transform, inverse, logjac,
    transform_and_logjac, @increment_logjac, ContinuousTransformation,
    # univariate transformations
    UnivariateTransformation,
    Affine, IDENTITY,
    Negation, NEGATION, Logistic, LOGISTIC, RealCircle, REALCIRCLE, Exp, EXP,
    Logit, LOGIT, InvRealCircle, INVREALCIRCLE, Log, LOG,
    affine_bridge, default_transformation, bridge,
    # multivariate transformations NOTE API may change
    UnitVector, CorrelationCholeskyFactor, lkj_correlation_cholesky_logpdf,
    # grouped transformations
    GroupedTransformation, get_transformation, ArrayTransformation,
    TransformationTuple,
    # wrapped transformations
    TransformationWrapper, TransformLogLikelihood, get_loglikelihood,
    TransformDistribution, get_distribution, logpdf_in_domain, logpdf_in_image,
    # utilities
    ungrouping_map, all_finite

import Base: in, length, size, ∘, show, getindex, intersect, extrema, minimum,
    maximum, isfinite, isinf, isapprox

import Statistics: middle


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
    $SIGNATURES

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

"""
    $SIGNATURES

Define [`transform`](@ref) and [`logjac`](@ref) using
[`transform_and_logjac`](@ref).
"""
macro define_from_transform_and_logjac(T)
    S = esc(T)
    quote
        (t::$S)(x) = transform_and_logjac(t, x)[1]
        ContinuousTransformations.logjac(t::$S, x) = transform_and_logjac(t, x)[2]
    end
end


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
$TYPEDEF

Continuous bijection ``D ⊂ ℝ^n→ I ⊂ ℝ^n`` or ``D ⊂ ℝ → I ⊂ ℝ``.
"""
abstract type ContinuousTransformation <: Function end

"""
    transform(t, x)

Return ``t(x)``. NOTE: this equivalent to the callable syntax `t(x)`, which is
also available.
"""
transform(t::ContinuousTransformation, x) = t(x)

"""
    y, logjac = $SIGNATURES

Return the transformed value and the determinant of the log Jacobian.

Equivalent to `t(x), logjac(t, x)`, but may be faster because of reused
components.
"""
function transform_and_logjac(t::ContinuousTransformation, x)
    transform(t, x), logjac(t, x)
end

"""
    $SIGNATURES

Return `transformation(x)` while incrementing `logjac` with the log Jacobian
determinant.
"""
macro increment_logjac(logjac, transformation, x)
    quote
        y, lj = transform_and_logjac($(esc(transformation)), $(esc(x)))
        $(esc(logjac)) += lj
        y
    end
end

"""
    $SIGNATURES

Return `transformation(x)` while incrementing `logjac` with the log Jacobian
determinant.
"""
macro increment_logjac(logjac, transformation_x)
    @capture transformation_x transformation_(x_)
    quote
        @increment_logjac $(esc(logjac)) $(esc(transformation)) $(esc(x))
    end
end

"""
    rhs_string(transformation, term)

Return the formula representing the hand side of the `transformation`, with
`term` as the argument.
"""
function rhs_string end

transformation_string(t, x = "x") = x * " ↦ " * rhs_string(t, x)

# this is necessary because `display` calls this method, which is defined for
# Function and would match
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
$TYPEDEF

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
        @argcheck isfinite(left) "Need finite endpoint."
        new(left)
    end
end

PositiveRay(left::T) where {T} = PositiveRay{T}(left)

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
        @argcheck isfinite(right) "Need finite endpoint."
        new(right)
    end
end

NegativeRay(right::T) where {T} = NegativeRay{T}(right)

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
        @argcheck isfinite(left) && isfinite(right) "Need finite endpoints."
        @argcheck left < right "Need strictly increasing endpoints."
        new(left, right)
    end
end

Segment(left::T, right::T) where {T <: Real} = Segment{T}(left, right)

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


# univariate transformations

"""
$TYPEDEF

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
$TYPEDEF

Trait that is useful for domain and image calculations. See [`RRStable`](@ref).
"""
abstract type RRStability end

"""
$TYPEDEF

Trait that indicates that a univariate transformation

1. maps ``ℝ`` to ``ℝ``,

2. supports mapping intervals, and

3. maps subtypes of [`AbstractInterval`](@ref) to the same type.
"""
struct RRStable <: RRStability end

"""
$TYPEDEF

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


# Affine and Negation
#
# These play a special role since map subtypes of `AbstractInterval` to
# themselves, and thus are useful for composition.

"""
    Affine(α, β)

Mapping ``ℝ → ℝ`` using ``x ↦ α⋅x + β``.

``α > 0`` is enforced, see [`Negation`](@ref).
"""
@auto_hash_equals struct Affine{T <: Real} <: UnivariateTransformation
    α::T
    β::T
    function Affine{T}(α::T, β::T) where T
        @argcheck α > 0 DomainError
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
    α == 1 || (term = "$(round(α; sigdigits = 4))⋅" * term)
    β == 0 || (term = term * " + $(round(β, sigdigits = 4))")
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


# transformations from ℝ to subsets

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


# transformations to ℝ from subsets

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


# composed transformations

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


# calculated transformations

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
default_transformation(::RealLine, ::Segment) = LOGISTIC
default_transformation(::RealLine, ::PositiveRay) = EXP
default_transformation(::RealLine, ::NegativeRay) = EXP

default_transformation(::Segment, ::RealLine) = LOGIT
default_transformation(::PositiveRay, ::RealLine) = LOG
default_transformation(::NegativeRay, ::RealLine) = LOG

bridge(dom::RealLine, img, t = default_transformation(dom, img)) =
    affine_bridge(image(t), img) ∘ t

bridge(dom, img::RealLine, t = default_transformation(dom, img)) =
    t ∘ affine_bridge(dom, domain(t))

bridge(::RealLine, ::RealLine) = IDENTITY


# Special constants

"""
Hyperbolic tangent transformation.

An affine stretch of LOGISTIC to (-1, 1).
"""
const TANH = bridge(ℝ, Segment(-1.0, 1.0), LOGISTIC)


# UnitVector

"""
    UnitVector(n)

Transform `n-1` real numbers to a unit vector of length `n`, under the Euclidean
norm.
"""
struct UnitVector <: ContinuousTransformation
    n::Int
    function UnitVector(n::Int)
        @argcheck n ≥ 1 "Dimension should be positive."
        new(n)
    end
end

length(t::UnitVector) = t.n - 1

rhs_string(t::UnitVector, term) = "UnitVector($(t.n))($term)"

function transform_and_logjac!(ys::AbstractVector{T},
                               t::UnitVector, xs::AbstractVector) where T
    @unpack n = t
    @argcheck length(xs) == n - 1
    @argcheck length(ys) == n
    remainder = one(T)
    logjac = zero(T)
    for (i, x) in enumerate(xs)
        z = @increment_logjac logjac TANH(x)
        ys[i] = z * √remainder
        logjac += log(remainder) / 2
        remainder -= abs2(z) * remainder
    end
    ys[end] = √remainder
    ys, logjac
end

transform_and_logjac(t::UnitVector, zs::AbstractVector{T}) where T =
    transform_and_logjac!(Vector{T}(undef, t.n), t, zs)

function inverse!(xs::AbstractVector{T},
                  t::UnitVector, ys::AbstractVector) where T
    @unpack n = t
    @argcheck length(xs) == n - 1
    @argcheck length(ys) == n
    remainder = one(T)
    for i in 1:(n - 1)
        y = ys[i]
        z = y / √remainder
        xs[i] = inverse(TANH, z)
        remainder -= abs2(y)
    end
    xs
end

inverse(t::UnitVector, ys::AbstractVector{T}) where T =
    inverse!(Vector{T}(undef, t.n - 1), t, ys)

@define_from_transform_and_logjac UnitVector


# Cholesky factor of a correlation matrix

"""
    CorrelationCholeskyFactor(n)

Cholesky factor of a correlation matrix of size `n`.
"""
struct CorrelationCholeskyFactor <: ContinuousTransformation
    n::Int
    function CorrelationCholeskyFactor(n)
        @argcheck n ≥ 1 "Dimension should be positive."
        new(n)
    end
end

"""
    $SIGNATURES

Number of elements in a triangle of an ``n×n`` square matrix, including the
diagonal.
"""
triangle_length(n) = n * (n-1) ÷ 2

length(t::CorrelationCholeskyFactor) = triangle_length(t.n)

rhs_string(t::CorrelationCholeskyFactor, term) =
    "CorrelationCholeskyFactor($(t.n))($term)"

function transform_and_logjac(t::CorrelationCholeskyFactor,
                              x::AbstractVector{T}) where T
    @unpack n = t
    @argcheck length(x) == length(t)
    logjac = zero(T)
    cumulative_index = 0
    L = zeros(T, n, n)
    for i in 1:n
        k = i - 1
        _, lj = transform_and_logjac!(@view(L[i, 1:i]), UnitVector(i),
                                      @view(x[cumulative_index .+ (1:k)]))
        logjac += lj
        cumulative_index += k
    end
    LowerTriangular(L), logjac
end

@define_from_transform_and_logjac CorrelationCholeskyFactor

function inverse(t::CorrelationCholeskyFactor,
                 L::LowerTriangular{T}) where T
    @unpack n = t
    @argcheck size(L) == (n, n)
    x = Vector{T}(undef, triangle_length(n))
    cumulative_index = 0
    for i in 1:n
        k = i - 1
        inverse!(@view(x[cumulative_index .+ (1:k)]), UnitVector(i),
                 @view(L[i, 1:i]))
        cumulative_index += k
    end
    x
end


# helper function for Cholesky decompositions

"""
    $SIGNATURES

Log PDF of the LKJ distribution (Lewandowski et al 2009) for correlation
matrices.

A correlation matrix ``Ω=LL'`` has the density ``|Ω|^(η-1)``. However, it is
usually not necessary to construct ``Ω``, so this function is parametrized in
terms of a Cholesky decomposition `L`.

Note that this function **does not check if `L` yields a valid correlation
matrix**. Use [`CorrelationCholeskyFactor`](@ref) for generating valid
arguments.

Valid values are ``η > 0``. When ``η > 1``, the distribution is unimodal at the
identity, while ``0 < η < 1`` has a trough. ``η = 2`` is recommended as a vague
prior.

When ``η = 1``, the density is uniform. Note however that the function should
still be invoked, because of the Jacobian correction of the transformation.
"""
function lkj_correlation_cholesky_logpdf(L::Union{LowerTriangular,
                                                  UpperTriangular}, η)
    @argcheck η > 0
    z = diag(L)
    n = size(L, 1)
    sum(log.(z) .* ((n:-1:1) .+ 2*(η-1))) + log(2) * n
end


# grouped transformations

"""
$TYPEDEF

Abstract type for grouped transformations.

A grouped transformation takes a vector, and transforms contiguous blocks of
elements to some output type, determined by the specific transformation type.

All subtypes support

- `length`: return the length of the vector that can be used as an argument

- callable object for the transformation

- `logjac`, and `inverse`,

- `domain` and `image`, which may have specific interpretation for their result
  types depending on the concrete subtype.
"""
abstract type GroupedTransformation <: ContinuousTransformation end

"""
    $SIGNATURES

Return the transformation from a wrapper object.
"""
get_transformation(d) = d.transformation


# array transformations

"""
    ArrayTransformation(transformation, dims)
    ArrayTransformation(transformation, dims...)

Apply transformation to a vector, returning an array of the given dimensions.

[`domain`](@ref), [`image`](@ref), and [`isincreasing`](@ref) return the
corresponding values for the underlying transformation.
"""
struct ArrayTransformation{T <: UnivariateTransformation,
                           D} <: GroupedTransformation
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
    @argcheck length(x) == length(t) DimensionMismatch()
    lj(x) = logjac(t.transformation, x)
    sum(lj, x)
end

inverse(t::ArrayTransformation, x) =
    reshape(inverse.(t.transformation, x), size(t))

@forward ArrayTransformation.transformation isincreasing, image, domain


# transformation tuple

"""
    TransformationTuple(transformations::Tuple)
    TransformationTuple(transformations...)

A tuple of [`ContinuousTransformation`](@ref)s. Given a vector of matching
`length`, each takes as many reals as needed, and returns the result as a tuple.
"""
@auto_hash_equals struct
    TransformationTuple{T <: Tuple{Vararg{ContinuousTransformation}}} <:
        GroupedTransformation
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

domain(t::TransformationTuple) = domain.(t.transformations)

image(t::TransformationTuple) = image.(t.transformations)

getindex(t::TransformationTuple, ix) = t.transformations[ix]

function next_indexes(acc, t)
    l = length(t)
    acc + l, acc .+ (1:l)
end

next_indexes(acc, t::UnivariateTransformation) = acc+1, acc+1

#=@unroll=# function transformation_indexes(ts)
    acc = 0
    result = ()
    #=@unroll=# for t in ts
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


# wrapper

"""
$TYPEDEF

Wrap a transformation to achieve some specialized functionality.

Supports `length`, [`get_transformation`](@ref), and other methods depending on
the subtype.
"""
abstract type TransformationWrapper <: Function end


# loglikelihood wrapper

"""
    TransformLogLikelihood(ℓ, transformation::Union{Tuple, GroupedTransformation})

    TransformLogLikelihood(ℓ, transformations...)

Return a callable that

1. transforms its vector argument using a grouped transformation to a set of
   values,

2. calls `ℓ` (which should return a scalar) with this tuple.

3. returns the result above corrected by the log Jacobians.

Useful when `ℓ` is a log-likelihood function with a restricted domain, and
`transformations` is used to trasform to this domain from ``ℝ^n``.

See also [`get_transformation`](@ref), [`get_distribution`](@ref),
`Distributions.logpdf`, and [`logpdf_in_domain`](@ref).
"""
struct TransformLogLikelihood{L, T <: GroupedTransformation} <:
        TransformationWrapper
    loglikelihood::L
    transformation::T
end

TransformLogLikelihood(L, T::Tuple) =
    TransformLogLikelihood(L, TransformationTuple(T))

TransformLogLikelihood(L, ts::ContinuousTransformation...) =
    TransformLogLikelihood(L, TransformationTuple(ts))

"""
    $SIGNATURES

Return the log likelihood function.
"""
get_loglikelihood(t::TransformLogLikelihood) = t.loglikelihood

@forward TransformLogLikelihood.transformation length

(f::TransformLogLikelihood)(x) =
    f.loglikelihood(f.transformation(x)) + logjac(f.transformation, x)

function show(io::IO, f::TransformLogLikelihood)
    print(io, "TransformLogLikelihood of length $(length(f)), with ")
    print(io, f.transformation)
end


# transform distributions

"""
    TransformDistribution(distribution, transformation)

Given a `transformation` and a `distribution`, create a transformed distribution
object that has the distribution of `transformation(x)` with `x ∼ distribution`.

The transformation object is callable with the same syntax as
`transformation`. It also supports methods `rand`, `length`.

See also [`logpdf_in_domain`](@ref) and [`logpdf_in_image`](@ref).
"""
struct TransformDistribution{D <: Union{ContinuousUnivariateDistribution,
                                        ContinuousMultivariateDistribution},
                             T <: Union{UnivariateTransformation,
                                        GroupedTransformation}} <:
                                            TransformationWrapper
    distribution::D
    transformation::T
    function TransformDistribution(distribution::D, transformation::T) where
        {D <: ContinuousUnivariateDistribution, T <: UnivariateTransformation}
        new{D, T}(distribution, transformation)
    end
    function TransformDistribution(distribution::D, transformation::T) where
        {D <: ContinuousMultivariateDistribution, T <: GroupedTransformation}
        @argcheck(length(distribution) == length(transformation),
                  "Incompatible lengths.")
        new{D,T}(distribution, transformation)
    end
end

get_transformation(t::TransformDistribution) = t.transformation

"""
    $SIGNATURES

Return the wrapped distribution.
"""
get_distribution(t::TransformDistribution) = t.distribution

(t::TransformDistribution)(x) = t.transformation(x)

rand(t::TransformDistribution) = t(rand(t.distribution))

@forward TransformDistribution.distribution length

"""
    $SIGNATURES

The log pdf for a transformed distribution at `t(x)` in image, calculated in the
domain without performing the transformation.

The log pdf is adjusted with the log determinant of the Jacobian, ie the
following holds:

```julia logpdf_in_image(t, t(x)) == logpdf_in_domain(t, x) ```

See [`logpdf_in_image`](@ref).

!!! note

    Typical usage of this function would be drawing some random `x`s from the
    contained distribution (possibly also used for some other purpose), and
    obtaining the log pdfs at `t(y)` with the same values.
"""
function logpdf_in_domain(t::TransformDistribution, x)
    # NOTE adjusting by -logjac because of the derivative of the inverse rule
    logpdf(t.distribution, x) - logjac(t.transformation, x)
end

"""
    $SIGNATURES

The log pdf for a transformed distribution at `y` in image.

See also [`logpdf_in_domain`](@ref).
"""
function logpdf_in_image(t::TransformDistribution, y)
    x = inverse(t.transformation, y)
    logpdf_in_domain(t, x)
end


# misc utilities
# NOTE: may be experimental, removed/redesigned
# NOTE: code written for simplicity, not speed, benchmark & optimize later

"""
    $SIGNATURES

Make a container for ungrouped results.

The first argument, `Vector` or `Array`, determines the result of ungrouping.

`representative_elt` is a representative element, used to determine element
type and stucture, and `len` is the length of the result.
"""
function make_ungrouped_container(::Type{Vector},
                                  representative_elt::AbstractArray{T},
                                  len) where T
    [Vector{T}(undef, len) for _ in CartesianIndices(axes(representative_elt))]
end

function make_ungrouped_container(::Type{Array},
                                  representative_elt::AbstractArray{T},
                                  len) where T
    Array{T}(undef, len, size(representative_elt)...)
end

make_ungrouped_container(T, representative_elt::Tuple, len) =
    map(x -> make_ungrouped_container(T, x, len), representative_elt)

make_ungrouped_container(_, representative_elt::T, len) where {T <: Real} =
    Vector{T}(undef, len)

"""
    $SIGNATURES

Save `elt` in `container` at index `i`. `T` determines the layout of
`container`, conformably with [`make_ungrouped_container`](@ref).
"""
ungroup_elt!(T, container::NTuple{N, Any}, elt::NTuple{N, Any}, i) where N =
    foreach((c, e) -> ungroup_elt!(T, c, e, i), container, elt)

function ungroup_elt!(::Type{Vector}, container::AbstractArray,
                      elt::AbstractArray, i)
    @argcheck axes(container) == axes(elt)
    for j in 1:length(elt)
        ungroup_elt!(Vector, container[j], elt[j], i)
    end
end

ungroup_elt!(::Type{Array}, container::AbstractArray,
             elt::AbstractArray{T, N}, i) where {T, N} =
    container[i, fill(Colon(), N)...] .= elt

ungroup_elt!(_, container::AbstractVector, elt::Number, i) = container[i] = elt

"""
    $SIGNATURES

Map the vector `A` elementwise using `f`, and collect the results in arrays or
vectors organized as a single result from `f`. For example, if `f` returns
tuples, this function returns tuples of vectors.

The shape of the result is determined by the first argument, which can be
`Vector` or `Array`. The difference is that `Array` stacks arrays into arrays,
while `Vector` does not.

`Array` may be more useful for summaries along dimensions, eg `mean` and
similar.
"""
function ungrouping_map(T, f, A::AbstractVector)
    n = length(A)
    representative_elt = f(A[1])
    container = make_ungrouped_container(T, representative_elt, n)
    ungroup_elt!(T, container, representative_elt, 1)
    for i in 2:n
        ungroup_elt!(T, container, f(A[i]), i)
    end
    container
end

"""
    $SIGNATURES

Test if a numerical argument is made up of finite real numbers.
"""
all_finite(x::Real) = isfinite(x)

all_finite(A::AbstractArray{T}) where {T <: Real} = all(isfinite, A)

all_finite(t::Tuple) = all(all_finite, t)

all_finite(itr) = all(all_finite, itr)

end
