using Documenter
using FixedPointArithmetic

DocMeta.setdocmeta!(FixedPointArithmetic, :DocTestSetup, :(using FixedPointArithmetic); recursive=true)

makedocs(;
    modules=[FixedPointArithmetic],
    authors="Jacob Gummesson Atroshi <jacob.gummesson@uni.lu>",
    sitename="FixedPointArithmetic.jl",
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true", edit_link="main"),
    pages=["Home" => "index.md", "Arithmetic model" => "arithmetic.md", "API reference" => "api.md"],
)

# The repository is taken from GitHub Actions, so the docs deploy from any fork without edits.
deploydocs(; repo="github.com/$(get(ENV, "GITHUB_REPOSITORY", "OWNER/FixedPointArithmetic.jl")).git", devbranch="main")
