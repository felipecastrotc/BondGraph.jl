using Pkg

if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end

if isdir("./docs/build")
    rm("./docs/build", recursive = true)
end
mkdir("./docs/build")

ENV["JULIA_DEBUG"] = "Documenter"

using Documenter, DocStringExtensions, Literate, BondGraph

# utility function from https://github.com/JuliaOpt/Convex.jl/blob/master/docs/make.jl
fix_math_md(content) = replace(content, r"\$\$(.*?)\$\$"s => s"```math\1```")

# Add a text after the first heading
function postprocess(cont)
    # Text to insert
    text = """
     *The source files for all examples can be found in [/examples](https://github.com/oxfordcontrol/COSMO.jl/tree/master/examples/).*\n
    """

    # Regular expression pattern to match the first heading
    pattern = r"(?m)^#.*?\n\n"

    # Find the first heading
    heading_match = match(pattern, cont)

    if heading_match !== nothing
        # Get the position after the first heading
        pos = heading_match.offset + length(heading_match.match)
        # Insert the text at the determined position
        cont = string(cont[1:pos-1], text, cont[pos:end])
    end
    return cont
end

function insert_quick_example(quick_example)
    # Load the quick example file
    quickex = read(joinpath(@__DIR__, "src", quick_example[1][2]), String)

    quickex = replace(quickex, "This page was generated using" => "This section was generated using")

    # Regular expression pattern to match the first heading
    pattern = r"(?m)^#.*?\n\n"

    # Find the first heading to clean the text before it
    heading_match = match(pattern, f)
    quickex = quickex[heading_match.offset:end]

    # Load the gettting_started.md
    path_getstarted = joinpath(@__DIR__, "src", "getting_started.md")
    getstarted = read(path_getstarted, String)

    # Find the quick-example position
    pos = findfirst("@quick-example", getstarted)

    cont = string(getstarted[1:pos[1]-1], quickex, getstarted[pos[end] + 1:end])

    write(joinpath(@__DIR__, "src","getting_started_auto.md"), cont)
end

# Generate the example navigation bar 
function get_nav(filename; suffix = "./examples/")
    f = open(suffix * filename, "r")
    name = split(split(readline(f), "[")[2], "]")[1]
    close(f)
    path = suffix * replace(filename, ".jl" => ".md")
    return String(name) => path
end

# code adapted from https://github.com/oxfordcontrol/COSMO.jl/blob/master/docs/make.jl
# find all example source files
example_path = joinpath(@__DIR__, "../examples/")
build_path = joinpath(@__DIR__, "src", "examples/")
files = readdir(example_path)
filter!(x -> endswith(x, ".jl"), files)

for file in files
    Literate.markdown(
        example_path * file,
        build_path;
        preprocess = fix_math_md,
        postprocess = postprocess,
        documenter = true,
        credit = true,
    )
end

# get example paths
examples_nav = get_nav.(files; suffix = "./examples/")
# Get the quick example
quick_example = filter(x -> x[1] == "Quick Example", examples_nav)
insert_quick_example(quick_example)
# Remove the quick example from the examples navigation
filter!(x -> x[1] != "Quick Example", examples_nav)


makedocs(
    modules = [BondGraph],
    format=Documenter.HTML(prettyurls=false, sidebar_sitename=false),
    sitename = "BondGraph.jl",
    # doctest=true,
    clean = true,
    src = "docs/src",
    force = true,
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started_auto.md",
        "User Guide" => [
        ],
        "Examples" => examples_nav,
        "Library" => [
            "Public" => "lib/public.md"
            "Internals" => "lib/internals.md"
        ],
        "Citing BondGraphToolkit" => "citing.md",
        "Contributing" => "contributing.md",
        "License" => "license.md",
        # "Changelog" => "changelog.md",
    ],
    pagesonly = true,
)

# Cleaning
rm(joinpath(@__DIR__, "src", "getting_started_auto.md"))