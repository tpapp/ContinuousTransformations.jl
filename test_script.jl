Pkg.clone(pwd())
Pkg.build("ContinuousTransformations")
Pkg.test("ContinuousTransformations"; coverage=true)
