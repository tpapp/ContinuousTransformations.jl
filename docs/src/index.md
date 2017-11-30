# Overview

```@meta
CurrentModule = ContinuousTransformations
```

This package implements some canonically used continuous bijections (also known as a [homeomorphism](https://en.wikipedia.org/wiki/Homeomorphism)) between subsets of ``\mathbb{R}^n``. These are useful if you have a function

```math
f: \mathcal{X} \subset \mathbb{R}^n \to \mathcal{Y}
```

and would like to use it as a building block to define some

```math
g: \mathcal{Z} \subset \mathbb{R}^n \to \mathcal{Y}
```

This package helps you find a function ``h`` such that ``g = f \circ h`` or ``f = g \circ h``.

To make things concrete, consider the following examples.

### Example: Chebyshev polynomials

[Chebyshev polynomials](https://en.wikipedia.org/wiki/Chebyshev_polynomials) are defined on ``(-1, 1)``. If you want to approximate a function on some generic ``(a, b)`` interval, you will need to transform. Usually one uses something like
```math
y = \left(x - \frac{a+b}2\right)\cdot\frac{b-a}2
```
but calculating these things manually is tedious and error prone.

### Example: transformed multivariate normal

You want to characterize the joint distribution of some quantities

```math
x \ge 0,\quad a \le y \le b
```

for a statistical problem. A frequently used approach is to generate a multivariate normal

```math
z \sim N(\mu, \Sigma)
```

and then transform ``z_1`` to ``x``, and ``z_2`` to ``y`` such that the constraints above hold.

### Example: domain transformation for MCMC

You are using Bayesian statistics to estimate a model with a posterior that has constraints, eg for a variance ``\sigma > 0`` is required. You have an [algorithm](https://github.com/tpapp/DynamicHMC.jl) that can perform efficient MCMC for a log posterior

```math
\ell: \mathbb{R}^n \to \mathbb{R}
```

but to apply it, you need to transform from ``\mathbb{R}`` to ``(0, \infty)``. The log posterior should be adjusted by the [log determinant of the transformation's Jacobian](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant#Jacobian_determinant).

This package can help you with all of these.
