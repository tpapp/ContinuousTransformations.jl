export                     # only export constants for singleton types
    # general
    LOGJAC,
    JAC,
    UnivariateTransformation,
    domain, image, domain_in_image, isincreasing,
    integral_substitution,
    # univariate transformations
    LOGISTIC,
    LOGIT,
    EXP,
    LOG,
    ODDSRATIO,
    INVODDSRATIO,
    REALCIRCLE,
    INVREALCIRCLE,
    Affine,
    Power,
    # composition
    bridge

import Base: inv, show

import Compat: ∘                # replace with Base in 0.6

immutable LogJac end

const LOGJAC = LogJac()

immutable Jac end

const JAC = Jac()

"""
Univariate monotone transformation, either increasing or decreasing on
the whole domain.
"""
abstract UnivariateTransformation

"""
Return the domain of the transformation.
"""
function domain end

"""
Return the image of the transformation.
"""
function image end

"""
Return the largest interval in the domain of `f` such that its image is in `img`.
"""
function domain_in_image(f::UnivariateTransformation, img::AbstractInterval)
    inv(f)(img ∩ image(f))
end

"""
Test whether the transformation is monotone increasing.
"""
function isincreasing end

"""
Test whether the transformation is the identity. May not idenfity all
identity transformations.
"""
isidentity(::UnivariateTransformation) = false

# TODO use when https://github.com/JuliaLang/julia/issues/14919 is fixed, use
#
# function (f::T){T <: UnivariateTransformation}(x::AbstractInterval)
#     monotone_map_interval(f, x, isincreasing(f))
# end
#
# and remove the interval map from below.

"""
Convenience macro that defines methods for univariate transformation
where parameter are not needed for the calculation.

The first argument is of the form `T(x)`, provides the type and the
variable which contains the argument for mappings.

See examples for usage.
"""
macro univariate_transformation_definitions(Tx, keys_and_values)
    @capture Tx T_(x_)
    dict = block_key_value_dict(keys_and_values)
    forms = []
    f = esc(:f)
    maybe_form!(key, f) = haskey(dict, key) && push!(forms, f(dict[key]))
    maybe_form!(:domain, e -> :($(esc(:domain))(::$T) = $(e)))
    maybe_form!(:image, e -> :($(esc(:image))(::$T) = $(e)))
    maybe_form!(:inv, e -> :($(esc(:inv))(::$T) = $(e)))
    maybe_form!(:mapping, e -> :((::$T)($x) = $(e)))
    push!(forms, quote
          function ($f::$T)(x::AbstractInterval)
          monotone_map_interval($f, x, isincreasing($f))
          end
          end)
    maybe_form!(:isincreasing, e -> :($(esc(:isincreasing))($x::$T) = $(e)))
    maybe_form!(:jac, e -> :(($f::$T)($x, ::Jac) = ($f($x), $(e))))
    maybe_form!(:mapping_and_jac, e -> :(($f::$T)($x, ::Jac) = $(e)))
    maybe_form!(:logjac, e -> :(($f::$T)($x, ::LogJac) = ($f($x), $(e))))
    maybe_form!(:mapping_and_logjac, e -> :(($f::$T)($x, ::LogJac) = $(e)))
    maybe_form!(:show, e -> :($(esc(:show))(io::IO, ::$T) = print(io, $(e))))
    quote $(forms...) end
end

"""
Transform an integrand and a domain for integration using `t` as the
substitution. Return the transformed function and the domain.

Example:

```julia
f, D = integral_substitution(INVODDSRATIO, x->exp(-x^2), 0..Inf)
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
Calculate the odds ratio x/(1-x), check for domain.

For internal use. Contains workaround for
https://github.com/dpsanders/ValidatedNumerics.jl/issues/249
"""
@inline function _oddsratio{T <: Real}(x::T)
    @argcheck zero(T) ≤ x ≤ one(T) DomainError()
    x / (one(x)-x)
end

"Transform ℝ to (0,1) using the logistic function."
@define_singleton Logistic <: UnivariateTransformation

@univariate_transformation_definitions Logistic(x) begin
    domain = ℝ
    image = 𝕀
    mapping = one(x) / (one(x) + exp(-x))
    isincreasing = true
    mapping_and_jac = (ℓ = LOGISTIC(x); (ℓ, exp(-x) * ℓ^2))
    logjac = -x-2*log1pexp(-x)
    inv = LOGIT
    show = "x ↦ 1/(1+e⁻ˣ)"
end

"Transfrom (0,1) to ℝ using the logit function."
@define_singleton Logit <: UnivariateTransformation

@univariate_transformation_definitions Logit(x) begin
    domain = 𝕀
    image = ℝ
    mapping = log(_oddsratio(x))
    isincreasing = true
    jac = 1/(x*(1-x))
    logjac = -(log(x)+(log(1-x)))
    inv = LOGISTIC
    show = "x ↦ log(x/(1-x))"
end

"Maps ``(0,1)`` to ``(0, ∞)`` using ``y = x/(1-x)``."
@define_singleton OddsRatio <: UnivariateTransformation

@univariate_transformation_definitions OddsRatio(x) begin
    domain = 𝕀
    image = ℝ⁺
    mapping = _oddsratio(x)
    isincreasing = true
    jac = one(x)/((one(x)-x)^2)
    logjac = -2*log(1-x)
    inv = INVODDSRATIO
    show = "x ↦ x/(1-x)"
end

"Maps ``(0,∞)`` to ``(0, 1)`` using ``y = x/(1+x)``."
@define_singleton InvOddsRatio <: UnivariateTransformation

@univariate_transformation_definitions InvOddsRatio(x) begin
    domain = ℝ⁺
    image = 𝕀
    mapping = begin
        @argcheck x ≥ zero(x) DomainError()
        x == Inf ? one(x) : x/(1+x)
    end
    isincreasing = true
    jac = (1+x)^(-2)
    logjac = -2*log1p(x)
    inv = ODDSRATIO
    show = "x ↦ x/(1+x)"
end

"Map ``ℝ`` to ``(-1,1)`` using ``y = x/√(1+x^2)``."
@define_singleton RealCircle <: UnivariateTransformation

@univariate_transformation_definitions RealCircle(x) begin
    domain = ℝ
    image = -1..1
    mapping = x/√(1+x^2)
    isincreasing = true
    jac = (1+x^2)^(-1.5)
    logjac = -1.5*log1psq(x)
    inv = INVREALCIRCLE
    show = "x ↦ x/√(1+x²)"
end


"Map ``(-1,1)`` to ``ℝ`` using ``y = x/√(1-x²)``."
@define_singleton InvRealCircle <: UnivariateTransformation

@univariate_transformation_definitions InvRealCircle(x) begin
    domain = -1..1
    image = ℝ
    mapping = x/√(1-x^2)
    isincreasing = true
    jac = (1-x^2)^(-1.5)
    logjac = -1.5*(log1p(x)+log1p(-x))
    inv = REALCIRCLE
    show = "x ↦ x/√(1-x²)"
end

"Transform ℝ to the interval (0,∞), using the exponential function."
@define_singleton Exp <: UnivariateTransformation

@univariate_transformation_definitions Exp(x) begin
    domain = ℝ
    image = ℝ⁺
    mapping = exp(x)
    isincreasing = true
    mapping_and_jac = (ϵ = exp(x); (ϵ,ϵ))
    logjac = x
    inv = LOG
    show = "x ↦ eˣ"
end

"Transform (0,∞) to ℝ  using the logarithm function."
@define_singleton Log <: UnivariateTransformation

@univariate_transformation_definitions Log(x) begin
    domain = ℝ⁺
    image = ℝ
    mapping = log(x)
    isincreasing = true
    jac = 1/x
    mapping_and_logjac = (ℓ=log(x); (ℓ, -ℓ))
    inv = EXP
    show = "x ↦ log(x)"
end

"""
Transform ℝ to itself using ``y = α⋅x + β``.
"""
@auto_hash_equals immutable Affine{T <: Real} <: UnivariateTransformation
    α::T
    β::T
    function Affine(α, β)
        @argcheck α ≠ zero(T) DomainError()
        new(α, β)
    end
end

Affine{T}(α::T, β::T) = Affine{T}(α, β)

Affine(α, β) = Affine(promote(α, β)...)

show(io::IO, a::Affine) = print(io, isidentity(a) ? "x ↦ x" : "x ↦ $(a.α)⋅x + $(a.β)")

@univariate_transformation_definitions Affine(x) begin
    domain = ℝ
    image = ℝ
end

(a::Affine)(x) = fma(x, a.α, a.β)

isincreasing{T}(a::Affine{T}) = a.α > zero(T)

(a::Affine)(x, ::Jac) = a(x), abs(a.α)

(a::Affine)(x, ::LogJac) = a(x), log(abs(a.α))

function inv{T}(a::Affine{T})
    @unpack α, β = a
    Affine(one(T)/α, -β/α)
end

isidentity{T}(a::Affine{T}) = a.α == one(T) && a.β == zero(T)

"Transform ℝ⁺ to itself using `y = x^γ`."
@auto_hash_equals immutable Power{T <: Real} <: UnivariateTransformation
    γ::T
    function Power(γ)
        @assert γ > zero(γ)
        new(γ)
    end
end

Power{T}(γ::T) = Power{T}(γ)

show(io::IO, x::Power) = print(io, "x ↦ x^$(x.γ)")

@univariate_transformation_definitions Power(x) begin
    domain = ℝ⁺
    image = ℝ⁺
    isincreasing = true
end

function (p::Power)(x)
    @argcheck x ≥ zero(x) DomainError()
    x^p.γ
end

(p::Power)(x, ::Jac) = p(x), p.γ*x^(p.γ-1)

(p::Power)(x, ::LogJac) = p(x), log(p.γ)+(p.γ-1)*log(x)

inv(p::Power) = Power(1/p.γ)

"""
Compose two univariate transformations. Use the `∘` operator for
construction.
"""
immutable ComposedTransformation{Tf <: UnivariateTransformation,
                                 Tg <: UnivariateTransformation} <:
                                     UnivariateTransformation
    f::Tf
    g::Tg
end

(c::ComposedTransformation)(x) = c.f(c.g(x))

function show(io::IO, c::ComposedTransformation)
    show(io, c.f)
    println(io, " ∘ ")
    show(io, c.g)
end

function (c::ComposedTransformation)(x, ::Jac)
    y, g′x = c.g(x, JAC)
    fy, f′y = c.f(y, JAC)
    fy, f′y * g′x
end

function (c::ComposedTransformation)(x, ::LogJac)
    y, log_g′x = c.g(x, LOGJAC)
    fy, log_f′y = c.f(y, LOGJAC)
    fy, log_f′y + log_g′x
end

function domain(c::ComposedTransformation)
    domain_in_image(c.g, domain(c.f))
end

function image(c::ComposedTransformation)
    # y ∈ image(c) iff ∃ x: f(g(x)) = y and x ∈ domain(g), g(x) ∈ domain(f)
    @unpack f, g = c
    f(image(g) ∩ domain(f))
end

inv(c::ComposedTransformation) = ComposedTransformation(inv(c.g), inv(c.f))

function ∘(f::UnivariateTransformation, g::UnivariateTransformation)
    if isidentity(f)
        g
    elseif isidentity(g)
        f
    elseif f == inv(g)
        Affine(1.0,0.0)
    else
        ComposedTransformation(f, g)
    end
end

∘(f::Affine, g::Affine) = Affine(f.α*g.α, fma(f.α, g.β, f.β))

"""
Return a transformation that would map `dom` to `img`, via `mapping`.

For some cases, `mapping` can be omitted, and a default will be used.
"""
function bridge(dom::AbstractInterval, mapping::UnivariateTransformation,
                img::AbstractInterval)
    m_dom, m_img = domain(mapping), image(mapping)
    bridge(m_img, img) ∘ mapping ∘ bridge(dom, m_dom)
end

function bridge(mapping::UnivariateTransformation, img::AbstractInterval)
    bridge(image(mapping), img) ∘ mapping
end

function bridge(dom::AbstractInterval, mapping::UnivariateTransformation)
    mapping ∘ bridge(dom, domain(mapping))
end


function bridge{Tdom,Timg}(dom::Tdom, img::Timg)
    throw(ArgumentError("Can't bridge a $(Tdom) to a $(Timg) without a transformation."))
end

function bridge(dom::Segment, img::Segment)
    α = width(img) / width(dom)
    β = middle(img) - middle(dom) * α
    Affine(α, β)
end
    
bridge(dom::PositiveRay, img::PositiveRay) = Affine(1.0, img.left - dom.left)
bridge(dom::NegativeRay, img::NegativeRay) = Affine(1.0, img.right - dom.right)
bridge(dom::PositiveRay, img::NegativeRay) = Affine(-1.0, img.right + dom.left)
bridge(dom::NegativeRay, img::PositiveRay) = Affine(-1.0, img.left - dom.right)

bridge(::RealLine, ::RealLine) = Affine(1.0, 0.0)

function bridge{T <: Union{PositiveRay, NegativeRay}}(dom::T, img::Segment)
    bridge(dom, INVODDSRATIO, img)
end

function bridge{T <: Union{PositiveRay, NegativeRay}}(dom::Segment, img::T)
    bridge(dom, ODDSRATIO, img)
end

bridge(::RealLine, img::Segment) = bridge(LOGISTIC, img)
bridge(dom::Segment, ::RealLine) = bridge(dom, LOGIT)

bridge(::RealLine, img::PositiveRay) = bridge(EXP, img)
bridge(dom::PositiveRay, ::RealLine) = bridge(dom, LOG)
