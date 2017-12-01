var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "#Overview-1",
    "page": "Overview",
    "title": "Overview",
    "category": "section",
    "text": "CurrentModule = ContinuousTransformationsThis package implements some canonically used continuous bijections (also known as a homeomorphism) between subsets of mathbbR^n. These are useful if you have a functionf mathcalX subset mathbbR^n to mathcalYand would like to use it as a building block to define someg mathcalZ subset mathbbR^n to mathcalYThis package helps you find a function h such that g = f circ h or f = g circ h.To make things concrete, consider the following examples."
},

{
    "location": "#Example:-Chebyshev-polynomials-1",
    "page": "Overview",
    "title": "Example: Chebyshev polynomials",
    "category": "section",
    "text": "Chebyshev polynomials are defined on (-1 1). If you want to approximate a function on some generic (a b) interval, you will need to transform. Usually one uses something likey = left(x - fraca+b2right)cdotfracb-a2but calculating these things manually is tedious and error prone."
},

{
    "location": "#Example:-transformed-multivariate-normal-1",
    "page": "Overview",
    "title": "Example: transformed multivariate normal",
    "category": "section",
    "text": "You want to characterize the joint distribution of some quantitiesx ge 0quad a le y le bfor a statistical problem. A frequently used approach is to generate a multivariate normalz sim N(mu Sigma)and then transform z_1 to x, and z_2 to y such that the constraints above hold."
},

{
    "location": "#Example:-domain-transformation-for-MCMC-1",
    "page": "Overview",
    "title": "Example: domain transformation for MCMC",
    "category": "section",
    "text": "You are using Bayesian statistics to estimate a model with a posterior that has constraints, eg for a variance sigma  0 is required. You have an algorithm that can perform efficient MCMC for a log posteriorell mathbbR^n to mathbbRbut to apply it, you need to transform from mathbbR to (0 infty). The log posterior should be adjusted by the log determinant of the transformation's Jacobian.This package can help you with all of these."
},

{
    "location": "general/#",
    "page": "General API",
    "title": "General API",
    "category": "page",
    "text": ""
},

{
    "location": "general/#ContinuousTransformations.ContinuousTransformation",
    "page": "General API",
    "title": "ContinuousTransformations.ContinuousTransformation",
    "category": "Type",
    "text": "Continuous bijection D  ^n I  ^n or D    I  .\n\n\n\n"
},

{
    "location": "general/#ContinuousTransformations.domain",
    "page": "General API",
    "title": "ContinuousTransformations.domain",
    "category": "Function",
    "text": "domain(transformation)\n\nReturn the domain of the transformation.\n\n\n\n"
},

{
    "location": "general/#ContinuousTransformations.image",
    "page": "General API",
    "title": "ContinuousTransformations.image",
    "category": "Function",
    "text": "image(transformation)\n\nReturn the image of the transformation.\n\n\n\n"
},

{
    "location": "general/#ContinuousTransformations.logjac",
    "page": "General API",
    "title": "ContinuousTransformations.logjac",
    "category": "Function",
    "text": "logjac(t, x)\n\nThe log of the determinant of the Jacobian of t at x. ```\n\n\n\n"
},

{
    "location": "general/#ContinuousTransformations.inverse",
    "page": "General API",
    "title": "ContinuousTransformations.inverse",
    "category": "Function",
    "text": "inverse(t, x)\n\nReturn t(x).\n\ninverse(t)\n\nReturn the transformation t.\n\n\n\n"
},

{
    "location": "general/#ContinuousTransformations.bridge",
    "page": "General API",
    "title": "ContinuousTransformations.bridge",
    "category": "Function",
    "text": "bridge(dom, img, [transformation])\n\nReturn a transformation that maps dom to img.\n\nThe transformation argument may be used to specify a particular transformation family, otherwise default_transformation is used.\n\n\n\n"
},

{
    "location": "general/#General-interface-for-transformations-1",
    "page": "General API",
    "title": "General interface for transformations",
    "category": "section",
    "text": "CurrentModule = ContinuousTransformationsTransformations are function-like objects, in the sense that they are callable. They also support the following general interface.ContinuousTransformation\ndomain\nimage\nlogjac\ninverseYou can create a transformation using the appropriate constructors, combine univariate-transformations, and create a transformation between two intervals.bridge"
},

{
    "location": "univariate/#",
    "page": "Intervals and univariate transformations",
    "title": "Intervals and univariate transformations",
    "category": "page",
    "text": ""
},

{
    "location": "univariate/#ContinuousTransformations.AbstractInterval",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.AbstractInterval",
    "category": "Type",
    "text": "Abstract supertype for all univariate intervals. It is not specified whether they are open or closed.\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.RealLine",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.RealLine",
    "category": "Type",
    "text": "RealLine()\n\nThe real line. Use the constant ℝ.\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.ℝ",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.ℝ",
    "category": "Constant",
    "text": "A constant for the real line.\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.PositiveRay",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.PositiveRay",
    "category": "Type",
    "text": "PositiveRay(left)\n\nThe real numbers above left. See ℝ⁺.\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.ℝ⁺",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.ℝ⁺",
    "category": "Constant",
    "text": "The positive real numbers.\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.NegativeRay",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.NegativeRay",
    "category": "Type",
    "text": "NegativeRay(right)\n\nThe real numbers below right. See ℝ⁻.\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.ℝ⁻",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.ℝ⁻",
    "category": "Constant",
    "text": "The negative real numbers.\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.Segment",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.Segment",
    "category": "Type",
    "text": "Segment(left, right)\n\nThe real numbers between left and right, with -  textleft  textright   enforced.\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.width",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.width",
    "category": "Function",
    "text": "width(s)\n\n\nWidth of a finite interval.\n\n\n\n"
},

{
    "location": "univariate/#intervals-1",
    "page": "Intervals and univariate transformations",
    "title": "Intervals",
    "category": "section",
    "text": "CurrentModule = ContinuousTransformationsThe interval types are different from some other interval implementations in Julia. They do not specify if the interval is open or closed at an endpoint, and also encode infiniteness and semi-infiniteness in the type, for type stable code.AbstractInterval\nRealLine\nℝ\nPositiveRay\nℝ⁺\nNegativeRay\nℝ⁻\nSegmentIntervals also support the following methods in Base: minimum, maximum, in, isfinite, isinf, extrema.Segments also support middle, linspace, andwidth"
},

{
    "location": "univariate/#univariate-transformations-1",
    "page": "Intervals and univariate transformations",
    "title": "Univariate transformations",
    "category": "section",
    "text": ""
},

{
    "location": "univariate/#ContinuousTransformations.UnivariateTransformation",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.UnivariateTransformation",
    "category": "Type",
    "text": "Univariate monotone transformation, either increasing or decreasing on the whole domain (thus, a bijection).\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.isincreasing",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.isincreasing",
    "category": "Function",
    "text": "isincreasing(transformation)\n\nReturn true (false), when the transformation is monotonically increasing (decreasing).\n\n\n\n"
},

{
    "location": "univariate/#General-interface-1",
    "page": "Intervals and univariate transformations",
    "title": "General interface",
    "category": "section",
    "text": "UnivariateTransformation\nisincreasing"
},

{
    "location": "univariate/#ContinuousTransformations.Affine",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.Affine",
    "category": "Type",
    "text": "Affine(α, β)\n\nMapping    using x  x + .\n\n  0 is enforced, see Negation.\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.IDENTITY",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.IDENTITY",
    "category": "Function",
    "text": "Identity (as an affine transformation).\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.Negation",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.Negation",
    "category": "Type",
    "text": "Negation()\n\nMapping    using x  -x.\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.NEGATION",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.NEGATION",
    "category": "Function",
    "text": "Negation()\n\nMapping    using x  -x.\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.Logistic",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.Logistic",
    "category": "Type",
    "text": "Logistic()\n\nMapping   (01) using x  1(1+exp(-x)).\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.LOGISTIC",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.LOGISTIC",
    "category": "Function",
    "text": "Logistic()\n\nMapping   (01) using x  1(1+exp(-x)).\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.Exp",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.Exp",
    "category": "Type",
    "text": "Exp()\n\nMapping    using x  exp(x).\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.EXP",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.EXP",
    "category": "Function",
    "text": "Exp()\n\nMapping    using x  exp(x).\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.Logit",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.Logit",
    "category": "Type",
    "text": "Logit()\n\nMapping (01)   using x  log(x(1-x)).\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.LOGIT",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.LOGIT",
    "category": "Function",
    "text": "Logit()\n\nMapping (01)   using x  log(x(1-x)).\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.Log",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.Log",
    "category": "Type",
    "text": "Log()\n\nMapping    using x  exp(x).\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.LOG",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.LOG",
    "category": "Function",
    "text": "Log()\n\nMapping    using x  exp(x).\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.InvRealCircle",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.InvRealCircle",
    "category": "Type",
    "text": "InvRealCircle()\n\nMapping (-11)   using x  x(1-x^2).\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.INVREALCIRCLE",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.INVREALCIRCLE",
    "category": "Function",
    "text": "InvRealCircle()\n\nMapping (-11)   using x  x(1-x^2).\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.RealCircle",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.RealCircle",
    "category": "Type",
    "text": "RealCircle()\n\nMapping   (-11) using x  x(1+x^2).\n\n\n\n"
},

{
    "location": "univariate/#ContinuousTransformations.REALCIRCLE",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.REALCIRCLE",
    "category": "Function",
    "text": "RealCircle()\n\nMapping   (-11) using x  x(1+x^2).\n\n\n\n"
},

{
    "location": "univariate/#Specific-transformations-1",
    "page": "Intervals and univariate transformations",
    "title": "Specific transformations",
    "category": "section",
    "text": "Affine\nIDENTITY\nNegation\nNEGATION\nLogistic\nLOGISTIC\nExp\nEXP\nLogit\nLOGIT\nLog\nLOG\nInvRealCircle\nINVREALCIRCLE\nRealCircle\nREALCIRCLE"
},

{
    "location": "univariate/#ContinuousTransformations.ComposedTransformation",
    "page": "Intervals and univariate transformations",
    "title": "ContinuousTransformations.ComposedTransformation",
    "category": "Type",
    "text": "ComposedTransformation(f, g)\n\nCompose two univariate transformations, resulting in the mapping fg, or `x ↦ f(g(x)).\n\nUse the ∘ operator for construction.\n\n\n\n"
},

{
    "location": "univariate/#Composing-transformations-1",
    "page": "Intervals and univariate transformations",
    "title": "Composing transformations",
    "category": "section",
    "text": "ComposedTransformation"
},

{
    "location": "grouped/#",
    "page": "Grouped transformations",
    "title": "Grouped transformations",
    "category": "page",
    "text": ""
},

{
    "location": "grouped/#ContinuousTransformations.get_transformation",
    "page": "Grouped transformations",
    "title": "ContinuousTransformations.get_transformation",
    "category": "Function",
    "text": "get_transformation(d)\n\n\nReturn the transformation from a wrapper object, eg TransformLogLikelihood.\n\n\n\n"
},

{
    "location": "grouped/#ContinuousTransformations.TransformationTuple",
    "page": "Grouped transformations",
    "title": "ContinuousTransformations.TransformationTuple",
    "category": "Type",
    "text": "TransformationTuple(transformations::Tuple)\nTransformationTuple(transformations...)\n\nA tuple of ContinuousTransformations. Given a vector of matching length, each takes as many reals as needed, and returns the result as a tuple.\n\n\n\n"
},

{
    "location": "grouped/#ContinuousTransformations.TransformLogLikelihood",
    "page": "Grouped transformations",
    "title": "ContinuousTransformations.TransformLogLikelihood",
    "category": "Type",
    "text": "TransformLogLikelihood(ℓ, transformations::Union{Tuple, TransformationTuple})\nTransformLogLikelihood(ℓ, transformations...)\n\nReturn a callable that\n\ntransforms its vector argument using a transformation tuple (or a tuple of\n\ntransformations, converted as required) to a tuple of values,\n\ncalls ℓ with these, which should return a scalar,\nreturns the result above corrected by the log Jacobians.\n\nUseful when ℓ is a log-likelihood function with a restricted domain, and transformations is used to trasform to this domain from ^n.\n\nSee also get_transformation.\n\n\n\n"
},

{
    "location": "grouped/#Grouped-transformations-1",
    "page": "Grouped transformations",
    "title": "Grouped transformations",
    "category": "section",
    "text": "CurrentModule = ContinuousTransformationsget_transformation\nTransformationTuple\nTransformLogLikelihood"
},

{
    "location": "internals/#",
    "page": "Internals",
    "title": "Internals",
    "category": "page",
    "text": ""
},

{
    "location": "internals/#Various-internals-1",
    "page": "Internals",
    "title": "Various internals",
    "category": "section",
    "text": "CurrentModule = ContinuousTransformations"
},

{
    "location": "internals/#ContinuousTransformations._maybe_segment",
    "page": "Internals",
    "title": "ContinuousTransformations._maybe_segment",
    "category": "Function",
    "text": "_maybe_segment(a, b)\n\n\nHelper function for forming a segment when possible. Internal, not exported.\n\n\n\n"
},

{
    "location": "internals/#Intervals-1",
    "page": "Internals",
    "title": "Intervals",
    "category": "section",
    "text": "_maybe_segment"
},

{
    "location": "internals/#ContinuousTransformations.RRStability",
    "page": "Internals",
    "title": "ContinuousTransformations.RRStability",
    "category": "Type",
    "text": "Trait that is useful for domain and image calculations. See RRStable.\n\n\n\n"
},

{
    "location": "internals/#ContinuousTransformations.RRStable",
    "page": "Internals",
    "title": "ContinuousTransformations.RRStable",
    "category": "Type",
    "text": "Trait that indicates that a univariate transformation\n\nmaps  to ,\nsupports mapping intervals, and\nmaps subtypes of AbstractInterval to the same type.\n\n\n\n"
},

{
    "location": "internals/#ContinuousTransformations.NotRRStable",
    "page": "Internals",
    "title": "ContinuousTransformations.NotRRStable",
    "category": "Type",
    "text": "Trait that indicates that a univariate transformation is not RRStable.\n\n\n\n"
},

{
    "location": "internals/#ContinuousTransformations.RR_stability",
    "page": "Internals",
    "title": "ContinuousTransformations.RR_stability",
    "category": "Function",
    "text": "RR_stability(?)\n\n\nReturn either the trait RRStable and NotRRStable.\n\n\n\n"
},

{
    "location": "internals/#ContinuousTransformations.composed_domain",
    "page": "Internals",
    "title": "ContinuousTransformations.composed_domain",
    "category": "Function",
    "text": "composed_domain(f_RR_stability, g_RR_stability, f, g)\n\n\n\n"
},

{
    "location": "internals/#ContinuousTransformations.composed_image",
    "page": "Internals",
    "title": "ContinuousTransformations.composed_image",
    "category": "Function",
    "text": "composed_image(f_RR_stability, g_RR_stability, f, g)\n\n\n\n"
},

{
    "location": "internals/#ContinuousTransformations.default_transformation",
    "page": "Internals",
    "title": "ContinuousTransformations.default_transformation",
    "category": "Function",
    "text": "default_transformation(dom, img)\n\nReturn a transformation from dom that can be mapped to img using affine_bridge.\n\n\n\n"
},

{
    "location": "internals/#ContinuousTransformations.affine_bridge",
    "page": "Internals",
    "title": "ContinuousTransformations.affine_bridge",
    "category": "Function",
    "text": "affine_bridge(interval1, interval1)\n\nReturn an affine transformation between two intervals of the same type.\n\n\n\n"
},

{
    "location": "internals/#Univariate-transformations-1",
    "page": "Internals",
    "title": "Univariate transformations",
    "category": "section",
    "text": "RRStability\nRRStable\nNotRRStable\nRR_stability\ncomposed_domain\ncomposed_imagedefault_transformation\naffine_bridge"
},

{
    "location": "internals/#ContinuousTransformations.rhs_string",
    "page": "Internals",
    "title": "ContinuousTransformations.rhs_string",
    "category": "Function",
    "text": "rhs_string(transformation, term)\n\nReturn the formula representing the hand side of the transformation, with term as the argument.\n\n\n\n"
},

{
    "location": "internals/#Printing-1",
    "page": "Internals",
    "title": "Printing",
    "category": "section",
    "text": "rhs_string"
},

{
    "location": "internals/#ContinuousTransformations.@define_isapprox",
    "page": "Internals",
    "title": "ContinuousTransformations.@define_isapprox",
    "category": "Macro",
    "text": "@define_isapprox(T, fields)\n\n\nDefine an isapprox method, comparing the given fields in type T.\n\n\n\n"
},

{
    "location": "internals/#ContinuousTransformations.@define_singleton",
    "page": "Internals",
    "title": "ContinuousTransformations.@define_singleton",
    "category": "Macro",
    "text": "Define a singleton type with the given name and supertype (specified as name <: supertype), and a constant which defaults to the name in uppercase.\n\n\n\n"
},

{
    "location": "internals/#ContinuousTransformations._fma",
    "page": "Internals",
    "title": "ContinuousTransformations._fma",
    "category": "Function",
    "text": "_fma(x, y, z)\n\nPlaceholder for Base.fma until https://github.com/JuliaDiff/ReverseDiff.jl/issues/86 is fixed.\n\n\n\n"
},

{
    "location": "internals/#Utilities-1",
    "page": "Internals",
    "title": "Utilities",
    "category": "section",
    "text": "@define_isapprox\n@define_singleton\n_fma"
},

]}
