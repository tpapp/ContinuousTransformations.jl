module ContinuousTransformations

using StatsFuns
using ValidatedNumerics
using Parameters
import Base: inv
using ArgCheck

export
    # general
    LOGJAC,
    JAC,
    UnivariateTransformation,
    domain,
    integral_substitution,
    # univariate transformations
    Logistic,
    Logit,
    Exp,
    Log,
    OddsRatio,
    InvOddsRatio,
    Affine,
    Power

######################################################################
# general interface
######################################################################

immutable LogJac end

const LOGJAC = LogJac()

immutable Jac end

const JAC = Jac()

𝕀(T) = zero(T)..one(T)

in𝕀(x) = zero(x) ≤ x ≤ one(x)

ℝ(T) = T(-Inf)..T(Inf)

ℝ⁺(T) = zero(T)..T(Inf)

inℝ⁺(x) = zero(x) ≤ x

abstract UnivariateTransformation

"""
Return the domain of the transformation as an interval, with the given
type (defaults to Float64).
"""
domain(t::UnivariateTransformation) = domain(Float64, t)

"""
Transform an integrand and a domain for integration using `t` as the
substitution. Return the transformed function and the domain.

Example:

```julia
f, D = integral_substitution(InvOddsRatio(), x->exp(-x^2), 0..Inf)
```

will return values such that
``
\int_D f(x) dx = \int_0^\infty exp(-x^2) dx = √π/2
``
"""
function integral_substitution(t, f, domain)
    t⁻ = inv(t)
    function(y)
        x, jac = t⁻(y, JAC)
        f(x)*jac
    end, t(domain)
end

"""
Apply a monotone transformation to an interval by endpoints, using
correct rounding (depending, of course, on `f` respecting rounding
mode).
"""
@inline function map_interval_monotone{T}(f, x::Interval{T}, increasing::Bool = true)
    f_rounded(x, mode) = setrounding(()->f(x), T, mode)
    if increasing
        Interval(f_rounded(x.lo, RoundDown), f_rounded(x.hi, RoundUp))
    else
        Interval(f_rounded(x.hi, RoundDown), f_rounded(x.lo, RoundUp))
    end
end

"""
Given a monotone increasing function `f` that operates on scalars,
define a mode for intervals.
"""
macro lift_monotone_increasing(T)
    quote
        (f::$T)(x::Interval) = map_interval_monotone(f, x)
    end
end

######################################################################
# logistic
######################################################################

"Transform ℝ to (0,1) using the logistic function."
immutable Logistic <: UnivariateTransformation end

domain{T <: Real}(::Type{T}, ::Logistic) = ℝ(T)

(::Logistic)(x) = logistic(x)

@lift_monotone_increasing Logistic

(f::Logistic)(x, ::Jac) = (ℓ = f(x); (ℓ, exp(-x)*ℓ^2))

(f::Logistic)(x, ::LogJac) = f(x), -x-2*log1pexp(-x)

inv(::Logistic) = Logit()

######################################################################
# logit
######################################################################

"""
Transfrom (0,1) to ℝ using the logit function.
"""
immutable Logit <: UnivariateTransformation end

domain{T <: Real}(::Type{T}, ::Logit) = 𝕀(T)

function (::Logit)(x)
    @argcheck in𝕀(x) DomainError()
    logit(x)
end

@lift_monotone_increasing Logit

(f::Logit)(x, ::Jac) = f(x), 1/(x*(1-x))

(f::Logit)(x, ::LogJac) = f(x), -(log(x)+(log(1-x)))

inv(::Logit) = Logistic()

######################################################################
# odds ratio and its inverse
######################################################################

"""
Maps ``(0,1)`` to ``(0, ∞)`` using ``y = x/(1-x)``.
"""
immutable OddsRatio <: UnivariateTransformation end

domain{T <: Real}(::Type{T}, ::OddsRatio) = 𝕀(T)

function (::OddsRatio)(x)
    @argcheck in𝕀(x) DomainError()
    x/(one(x)-x)
end

@lift_monotone_increasing OddsRatio

(f::OddsRatio)(x, ::Jac) = f(x), one(x)/((one(x)-x)^2)

(f::OddsRatio)(x, ::LogJac) = f(x), -2*log(1-x)

inv(::OddsRatio) = InvOddsRatio()

"""
Maps ``(0,∞)`` to ``(0, 1)`` using ``y = x/(1+x)``.
"""
immutable InvOddsRatio <: UnivariateTransformation end

domain{T <: Real}(::Type{T}, ::InvOddsRatio) = ℝ⁺(T)

function (::InvOddsRatio)(x)
    @argcheck inℝ⁺(x) DomainError()
    x == Inf ? one(x) : x/(1+x)
end

@lift_monotone_increasing InvOddsRatio

(f::InvOddsRatio)(x, ::Jac) = f(x), (1+x)^(-2)

(f::InvOddsRatio)(x, ::LogJac) = f(x), -2*log1p(x)

inv(::InvOddsRatio) = OddsRatio()

######################################################################
# exponential and log
######################################################################

"Transform ℝ to the interval (0,∞), using the exponential function."
immutable Exp <: UnivariateTransformation end

domain{T <: Real}(::Type{T}, ::Exp) = ℝ(T)

(::Exp)(x) = exp(x)

(::Exp)(x, ::Jac) = (ϵ = exp(x); (ϵ,ϵ))

(::Exp)(x, ::LogJac) = exp(x), x

inv(::Exp) = Log()

"""
Transform (0,∞) to ℝ  using the logarithm function.
"""
immutable Log <: UnivariateTransformation end

domain{T <: Real}(::Type{T}, ::Log) = ℝ⁺(T)

(::Log)(x) = log(x)

(::Log)(x, ::Jac) = log(x), 1/x

(::Log)(x, ::LogJac) = (ℓ=log(x); (ℓ, -ℓ))

inv(::Log) = Exp()

######################################################################
# affine transformation
######################################################################

"""
Transform ℝ to itself using ``y = α⋅x + β``.
"""
immutable Affine{T <: Real} <: UnivariateTransformation
    α::T
    β::T
    function Affine(α, β)
        @argcheck α ≠ zero(T) DomainError()
        new(α, β)
    end
end

Affine{T}(α::T, β::T) = Affine{T}(α, β)

Affine(α, β) = Affine(promote(α, β)...)

domain{T <: Real}(::Type{T}, ::Affine) = ℝ(T)

(a::Affine)(x) = fma(x, a.α, a.β)

(a::Affine)(x::Interval) = map_interval_monotone(a, x, a.α > 0)

(a::Affine)(x, ::Jac) = a(x), abs(a.α)

(a::Affine)(x, ::LogJac) = a(x), log(abs(a.α))

function inv{T}(a::Affine{T})
    @unpack α, β = a
    Affine(one(T)/α, -β/α)
end

"""
Return an Affine map that maps the first interval to the second.
"""
function Affine(i1::Interval, i2::Interval)
    @argcheck isfinite(i1) && isfinite(i2) "infinite interval(s)"
    d1 = diam(i1)
    d2 = diam(i2)
    @argcheck d1 > 0 && d2 > 0 "empty interval(s)"
    α = d2 / d1
    β = i2.lo - i1.li * α
    Affine(α, β)
end

######################################################################
# power transformation
######################################################################

immutable Power{T <: Real} <: UnivariateTransformation
    γ::T
    function Power(γ)
        @assert γ > zero(γ)
        new(γ)
    end
end

Power{T}(γ::T) = Power{T}(γ)

domain{T}(::Type{T}, ::Power) = ℝ⁺(T)

function (p::Power)(x)
    @argcheck inℝ⁺(x) DomainError()
    x^p.γ
end

@lift_monotone_increasing Power

(p::Power)(x, ::Jac) = p(x), p.γ*x^(p.γ-1)

(p::Power)(x, ::LogJac) = p(x), log(p.γ)+(p.γ-1)*log(x)

inv(p::Power) = Power(1/p.γ)

end # module
