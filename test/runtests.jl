using ContinuousTransformations
import ForwardDiff: derivative
using Cubature
import Base: rand
using Base.Test
import Compat: ∘

include("test-utilities.jl")

include("test-intervals.jl")
include("test-transformations-basic.jl")
include("test-transformations-composed.jl")
include("test-integral.jl")

