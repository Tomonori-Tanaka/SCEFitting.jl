using SCEFitting
import Spglib   # activates the SpglibBackend extension for the executed `@example` blocks
using Documenter
using Documenter: Remotes

DocMeta.setdocmeta!(SCEFitting, :DocTestSetup, :(using SCEFitting);
                    recursive = true)

makedocs(;
    sitename = "SCEFitting.jl",
    modules = [SCEFitting],
    repo = Remotes.GitHub("Tomonori-Tanaka", "SCEFitting.jl"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        mathengine = Documenter.MathJax3(),
        canonical = "https://tomonori-tanaka.github.io/SCEFitting.jl/dev",
        edit_link = "main",
        footer = "Built with [Documenter.jl](https://documenter.juliadocs.org).",
        # api.md is one page listing the whole public surface, so it grows with the
        # API and periodically crosses Documenter's 200 KiB HTML threshold. Raised
        # rather than split: the index at the top is what makes a single page
        # navigable, and splitting it would scatter the cross-references every
        # docstring uses. [SLCE 2ba9140]
        size_threshold = 512 * 2^10,
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
        "Verification" => [
            "verification/angular_momentum.md",
        ],
        "API reference" => "api.md",
    ],
    warnonly = false,   # strict: any @example error / missing docstring fails the build
    checkdocs = :public,   # the unexported `public` surface is API too [SLCE 575a4e3]
    doctest = false,
)

# Publishes to https://tomonori-tanaka.github.io/SCEFitting.jl/ from the
# `documentation build` CI job (which needs `permissions: contents: write`). Outside
# CI this is a no-op, so a local `julia --project=docs docs/make.jl` still just
# builds into `docs/build/`.
deploydocs(;
    repo = "github.com/Tomonori-Tanaka/SCEFitting.jl",
    devbranch = "main",
    push_preview = false,
)
