using Pkg

if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end

if isdir("./docs/build")
    rm("./docs/build",recursive=true)
end
mkdir("./docs/build")

ENV["JULIA_DEBUG"] = "Documenter"

using Documenter, BondGraph

pages = [
    "Home" => "index.md",
    "Installation" => "installation.md",
    "Library" => [
        "Public" => "lib/public.md"
        "Internals" => "lib/internals.md"
    ],
    "Contributing" => "contributing.md"
]

makedocs(modules=[BondGraph],
    format=Documenter.HTML(prettyurls=false),
    sitename="BondGraph.jl",
    # doctest=true,
    # clean=true,
    src="docs/src",
    # force=true,
    pages=pages,
    pagesonly = true,
)