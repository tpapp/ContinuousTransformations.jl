Pkg.clone(pwd())
Pkg.clone("https://github.com/afniedermayer/InferenceUtilities.jl.git")
Pkg.build("ContinuousTransformations")
Pkg.test("ContinuousTransformations"; coverage=true)
