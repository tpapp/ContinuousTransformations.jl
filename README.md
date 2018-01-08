# ContinuousTransformations

[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.org/tpapp/ContinuousTransformations.jl.svg?branch=master)](https://travis-ci.org/tpapp/ContinuousTransformations.jl)
[![Coverage Status](https://coveralls.io/repos/tpapp/ContinuousTransformations.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/tpapp/ContinuousTransformations.jl?branch=master)
[![codecov.io](http://codecov.io/github/tpapp/ContinuousTransformations.jl/coverage.svg?branch=master)](http://codecov.io/github/tpapp/ContinuousTransformations.jl?branch=master)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://tpapp.github.io/ContinuousTransformations.jl/latest)

Continuous transformations (or more precisely, homeomorphisms) from ℝ (and two-point compactified version) and ℝⁿ to various open (or closed) sets used in statistics and numerical methods, such as intervals, simplexes, ordered vectors.

**Work in progress, API may change without notice.**

## Overview

This package was born because I was tired of coding the same transformations over and over, with occasional bugs, and wanted something well-tested.

Transformations defined by the package can be

1. called as functions,
2. provide the `logjac(transformation, x)` method for the log Jacobian determinant,
3. prodide the `inverse(transformation, x)` method for the inverse.

Log jacobian determinants and their log are useful for domain transformations in MCMC, among other things.

In addition, the package includes types to represent intervals, and some basic methods of working with them. The concept of intervals is slightly different from [IntervalSet.jl](https://github.com/JuliaMath/IntervalSets.jl) and [ValidatedNumerics.jl](https://github.com/dpsanders/ValidatedNumerics.jl), and as a result not compatible with either.

The convenience function `bridge(dom, img)` figures out the right transformation from `dom` to `img`. Currently implemented for intervals.

Examples:
```julia
using ContinuousTransformations
t = bridge(ℝ, Segment(0.0, 3.0)) # will use a real-circle transformation, stretched
t(0.0)             # 1.5
inverse(t, 1.5)    # ≈ 0.0
logjac(t, 0)       # ≈ 0.405
image(t)           # Segment(0.0, 3.0)
```

`ArrayTransformation(transformation, dimensions...)` transforms a vector of numbers to an array elementwise using `transformation`.

`TransformationTuple(transformations)` can be used for heterogeneous collections of transformations.

`TransformLogLikelihood` wraps a log likelihood function, and `TransformDistribution` transforms a distribution. Both of them take care of the log Jacobian determinant adjustment.

Special transformations, useful for Bayesian methods, are also available (WIP). Feature and pull requests are appreciated.

## Bibliography

- Stan Development Team (2017). "Modeling Language User's Guide and Reference Manual, Version 2.17.0" [(pdf)](http://mc-stan.org/users/documentation/)

- Lewandowski, Daniel, Dorota Kurowicka, and Harry Joe. "Generating random correlation matrices based on vines and extended onion method." Journal of multivariate analysis 100.9 (2009): 1989–2001.
