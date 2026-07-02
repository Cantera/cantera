# Build the Cantera.jl documentation with Documenter.
#
#   julia --project=docs docs/make.jl
#
# Requires a working libcantera (see docs/src/installation.md); the doctests and
# any `@example` blocks execute against it.

using Documenter
using Cantera

makedocs(
    sitename = "Cantera.jl",
    modules = [Cantera],
    authors = "Cantera Developers",
    pages = [
        "Home" => "index.md",
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    warnonly = true,   # the CLib backend is experimental; don't fail on missing refs
)
