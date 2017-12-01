using Documenter, ContinuousTransformations

makedocs(modules = [ContinuousTransformations],
         format = :html,
         clean = true,
         sitename = "ContinuousTransformations.jl",
         authors = "Tamás K. Papp",
         checkdocs = :all,
         linkcheck = true,
         assets = ["assets/custom.css"],
         html_prettyurls = haskey(ENV, "TRAVIS"), # clean URLs building on Travis
         pages = [
             "Overview" => "index.md",
             "General API" => "general.md",
             "Intervals and univariate transformations" => "univariate.md",
             "Grouped transformations" => "grouped.md",
             "Wrapped transformations" => "wrapped.md",
             "Internals" => "internals.md"
         ])

deploydocs(repo = "github.com/tpapp/ContinuousTransformations.jl.git",
           target = "build",
           deps = nothing,
           make = nothing,
           julia = "0.6")
