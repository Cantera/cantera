# Build the Cantera.jl documentation with Documenter.
#
#   julia --project=interfaces/julia/docs interfaces/julia/docs/make.jl
#
# The generated HTML is written to `docs/build`, which the Cantera doc build
# (doc/SConscript) copies into `build/doc/html/julia` so it is deployed
# alongside the rest of the Cantera documentation. There is intentionally no
# `deploydocs` call: deployment is handled by the top-level docs workflow.
#
# Requires a working libcantera and the generated CLib bindings (see the
# interface README); `using Cantera` loads the library to introspect docstrings.

using Documenter
using Cantera

makedocs(
    sitename = "Cantera.jl",
    modules = [Cantera],
    authors = "Cantera Developers",
    pages = [
        "Home" => "index.md",
        "API Reference" => "reference.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://cantera.org/dev/julia",
        # The API reference is a single @autodocs page; allow it to exceed the
        # default 200 KiB cap rather than splitting the generated listing.
        size_threshold = 512 * 1024,
    ),
    warnonly = true,   # the CLib backend is experimental; don't fail on missing refs
)
