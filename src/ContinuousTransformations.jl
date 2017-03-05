module ContinuousTransformations

using StatsFuns
using ValidatedNumerics
using Parameters
import Base: inv

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

ℝ(T) = T(-Inf)..T(Inf)

ℝ⁺(T) = zero(T)..T(Inf)

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

######################################################################
# logistic
######################################################################

"""
Transform ℝ to (0,1) using the logistic function.
"""
immutable Logistic <: UnivariateTransformation end

domain{T <: Real}(::Type{T}, ::Logistic) = ℝ(T)

(::Logistic)(x) = logistic(x)

(::Logistic)(x, ::Jac) = (ℓ = logistic(x); (ℓ, exp(-x)*ℓ^2))

(::Logistic)(x, ::LogJac) = logistic(x), -x-2*log1pexp(-x)

inv(::Logistic) = Logit()

######################################################################
# logit
######################################################################

"""
Transfrom (0,1) to ℝ using the logit function.
"""
immutable Logit <: UnivariateTransformation end

domain{T <: Real}(::Type{T}, ::Logit) = 𝕀(T)

(::Logit)(x) = logit(x)

(::Logit)(x, ::Jac) = logit(x), 1/(x*(1-x))

(::Logit)(x, ::LogJac) = logit(x), -(log(x)+(log(1-x)))
 
inv(::Logit) = Logistic()

######################################################################
# odds ratio and its inverse
######################################################################

"""
Maps ``(0,1)`` to ``(0, ∞)`` using ``y = x/(1-x)``.
"""
immutable OddsRatio <: UnivariateTransformation end

domain{T <: Real}(::Type{T}, ::OddsRatio) = 𝕀(T)

(::OddsRatio)(x) = x/(one(x)-x)

(::OddsRatio)(x, ::Jac) = x/(one(x)-x), one(x)/((one(x)-x)^2)

(::OddsRatio)(x, ::LogJac) = x/(1-x), -2*log(1-x)

inv(::OddsRatio) = InvOddsRatio()

"""
Maps ``(0,∞)`` to ``(0, 1)`` using ``y = x/(1+x)``.
"""
immutable InvOddsRatio <: UnivariateTransformation end

domain{T <: Real}(::Type{T}, ::InvOddsRatio) = ℝ⁺(T)

(::InvOddsRatio)(x) = x == Inf ? one(x) : x/(1+x)

(::InvOddsRatio)(x::Interval) = Interval(InvOddsRatio()(x.lo),
                                         InvOddsRatio()(x.hi))

(::InvOddsRatio)(x, ::Jac) = x/(1+x), (1+x)^(-2)

(::InvOddsRatio)(x, ::LogJac) = x/(1+x), -2*log1p(x)

inv(::InvOddsRatio) = OddsRatio()

######################################################################
# exponential and log
######################################################################

"""
Transform ℝ to the interval (0,∞), using the exponential function.
"""
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
end

domain{T <: Real}(::Type{T}, ::Affine) = ℝ(T)

(a::Affine)(x) = muladd(x, a.α, a.β)

(a::Affine)(x, ::Jac) = a(x), abs(a.α)

(a::Affine)(x, ::LogJac) = a(x), log(abs(a.α))

function inv{T}(a::Affine{T})
    @unpack α, β = a
    @assert α ≠ 0
    Affine(one(T)/α, -β/α)
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

(p::Power)(x) = (@assert x ≥ 0; x^p.γ)

(p::Power)(x, ::Jac) = p(x), p.γ*x^(p.γ-1)

(p::Power)(x, ::LogJac) = p(x), log(p.γ)+(p.γ-1)*log(x)

inv(p::Power) = Power(1/p.γ)

end # module
