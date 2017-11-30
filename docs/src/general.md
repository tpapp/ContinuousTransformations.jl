# General interface for transformations

```@meta
CurrentModule = ContinuousTransformations
```

Transformations are [function-like objects](https://docs.julialang.org/en/latest/manual/methods/#Function-like-objects-1), in the sense that they are *callable*. They also support the following general interface.

```@docs
ContinuousTransformation
domain
image
logjac
inverse
```

You can create a transformation using the appropriate constructors, combine [univariate-transformations](@ref), and create a transformation between two [intervals](@ref intervals).

```@docs
bridge
```
