# [Intervals](@id intervals)

The interval types are different from some other interval implementations in Julia. They do not specify if the interval is open or closed at an endpoint, and also encode infiniteness and semi-infiniteness in the type, for type stable code.

```@docs
AbstractInterval
RealLine
ℝ
PositiveRay
ℝ⁺
NegativeRay
Segment
```

Intervals also support the following methods in `Base`: `minimum`, `maximum`, `in`, `isfinite`, `isinf`, `extrema`.

[`Segment`](@ref)s also support `middle`, `linspace`, and

```@docs
width
```

# [Univariate transformations](@id univariate-transformations)
