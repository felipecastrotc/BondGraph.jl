using Pkg

if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end

if isdir("./docs/build")
    rm("./docs/build",recursive=true)
end
mkdir("./docs/build")

ENV["JULIA_DEBUG"] = "Documenter"

using Documenter, DocStringExtensions, Literate, BondGraph

function get_nav(filename; suffix="./examples/")
    file = replace(filename, ".jl" => "")
    name = join([uppercase(s[1]) * s[2:end] for s in split(file, "_")], " ")
    path = suffix*replace(filename, ".jl" => ".md")
    return name => path
end

# utility function from https://github.com/JuliaOpt/Convex.jl/blob/master/docs/make.jl
fix_math_md(content) = replace(content, r"\$\$(.*?)\$\$"s => s"```math\1```")
function postprocess(cont)
    """
    The source files for all examples can be found in [/examples](https://github.com/oxfordcontrol/COSMO.jl/tree/master/examples/).
    """ * cont
end

# code adapted from https://github.com/oxfordcontrol/COSMO.jl/blob/master/docs/make.jl
# find all example source files
example_path = joinpath(@__DIR__, "../examples/")
build_path = joinpath(@__DIR__, "src", "examples/")
files = readdir(example_path)
filter!(x -> endswith(x, ".jl"), files)

for file in files
    Literate.markdown(example_path * file, build_path; preprocess=fix_math_md, postprocess=postprocess, documenter=true, credit=true)
end

examples_nav = get_nav.(files; suffix="./examples/")

makedocs(modules=[BondGraph],
    format=Documenter.HTML(prettyurls=false),
    sitename="BondGraph.jl",
    # doctest=true,
    clean=true,
    src="docs/src",
    force=true,
    pages=[
            "Home" => "index.md",
            "Installation" => "man/installation.md",
            "Getting Started" => "man/getting_started.md",
            "Examples" => examples_nav,
            "Library" => [
                "Public" => "lib/public.md"
                "Internals" => "lib/internals.md"
            ],
            "Contributing" => "contributing.md",
            "Changelog" => "changelog.md",
        ],
    pagesonly = true,
)