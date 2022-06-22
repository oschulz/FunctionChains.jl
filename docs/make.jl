# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using FunctionChains

# Doctest setup
DocMeta.setdocmeta!(
    FunctionChains,
    :DocTestSetup,
    :(using FunctionChains);
    recursive=true,
)

makedocs(
    sitename = "FunctionChains",
    modules = [FunctionChains],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://oschulz.github.io/FunctionChains.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    strict = !("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/oschulz/FunctionChains.jl.git",
    forcepush = true,
    push_preview = true,
)
