Pkg.clone(pwd())
Pkg.clone("https://github.com/afniedermayer/InferenceUtilities.jl/blob/master/src/InferenceUtilities.jl")
Pkg.build("ContinuousTransformations")
Pkg.test("ContinuousTransformations"; coverage=true)
