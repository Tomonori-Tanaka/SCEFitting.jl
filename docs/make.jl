using SCEFitting
import Spglib   # activates the SpglibBackend extension for the executed `@example` blocks
using Documenter

DocMeta.setdocmeta!(SCEFitting, :DocTestSetup, :(using SCEFitting);
                    recursive = true)

makedocs(;
    sitename = "SCEFitting.jl",
    modules = [SCEFitting],
    # Local-only build: there is no published remote yet, so do not try to resolve
    # "edit on GitHub" / source links. Add a `repolink`/`deploydocs` when a remote exists.
    remotes = nothing,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        mathengine = Documenter.MathJax3(),
        edit_link = nothing,
        repolink = "",
        footer = "Built with [Documenter.jl](https://documenter.juliadocs.org).",
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Guide" => [
            "guide/basis.md",
            "guide/fitting.md",
            "guide/io.md",
            "guide/sunny.md",
        ],
        "Tutorials" => [
            "tutorials/index.md",
            "tutorials/heisenberg_chain.md",
            "tutorials/kagome_threebody.md",
            "tutorials/case1_bcc_fe.md",
        ],
        "Theory" => [
            "theory/index.md",
            "theory/sce.md",
            "theory/resolvability.md",
            "theory/architecture.md",
        ],
        "API reference" => "api.md",
    ],
    warnonly = true,
    checkdocs = :exports,
    doctest = false,
)
