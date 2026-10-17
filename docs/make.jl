
import Pkg

# Use this checkout of the package
Pkg.develop(path = joinpath(@__DIR__, ".."))
Pkg.instantiate()

# Import required packages
using Documenter, Literate, BondGraphToolkit

# Function to check if display is available
function has_display()
    return haskey(ENV, "DISPLAY") && !isempty(ENV["DISPLAY"])
end

# Set environment variables
if !has_display()
    ENV["GKSwstype"] = "100"  # Set GR to non-interactive mode
end
ENV["JULIA_DEBUG"] = "Documenter"

# Clean docs/build directory
build_dir = joinpath(@__DIR__, "build")
isdir(build_dir) && rm(build_dir, recursive = true)
mkdir(build_dir)


# Utility function to fix math in Markdown content
# https://github.com/JuliaOpt/Convex.jl/blob/master/docs/make.jl
fix_math_md(content) = replace(content, r"\$\$(.*?)\$\$"s => s"```math\1```")

# Function to add text after the first heading in Markdown content
function postprocess(cont)
    # Text to insert
    text = """
     *The source files for all examples can be found in [/examples](https://github.com/felipecastrotc/BondGraph.jl/tree/main/examples/).*\n
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

# Function to generate the example navigation bar
function insert_quick_example(quick_example)
    # Load the quick example file
    quickex = read(joinpath(@__DIR__, "src", quick_example[1][2]), String)

    quickex = replace(quickex, "This page was generated using" => "This section was generated using")

    # Regular expression pattern to match the first heading
    pattern = r"(?m)^#.*?\n\n"

    # Find the first heading to clean the text before it
    heading_match = match(pattern, quickex)
    quickex = quickex[heading_match.offset:end]

    # Load the gettting_started.md
    path_getstarted = joinpath(@__DIR__, "src", "getting_started.md")
    getstarted = read(path_getstarted, String)

    # Find the quick-example position
    pos = findfirst("@quick-example", getstarted)

    cont = string(getstarted[1:pos[1]-1], quickex, getstarted[pos[end]+1:end])

    write(joinpath(@__DIR__, "src", "getting_started_auto.md"), cont)
end


# Function to generate the example navigation bar
function get_nav(filename; suffix="./examples/")
    name = split(split(readline(joinpath(example_path, filename)), "[")[2], "]")[1]
    path = suffix * replace(filename, ".jl" => ".md")
    return String(name) => path
end

# Find all example source files and convert them to Markdown
# Adapted: https://github.com/oxfordcontrol/COSMO.jl/blob/master/docs/make.jl
example_path = joinpath(@__DIR__, "../examples/")
build_path = joinpath(@__DIR__, "src", "examples/")
files = readdir(example_path)
filter!(x -> endswith(x, ".jl"), files)
filter!((x) -> !occursin("tmp_", x),  files)

for file in files
    Literate.markdown(
        example_path * file,
        build_path;
        preprocess=fix_math_md,
        postprocess=postprocess,
        documenter=true,
        credit=true
    )
end

# Get example paths and remove the quick example from navigation
examples_nav = get_nav.(files; suffix="./examples/")
# examples_nav = Dict("Duffing Equation" => "./examples/mechanical_duffing.md")
quick_example = filter(x -> x[1] == "Quick Example", examples_nav)
insert_quick_example(quick_example)
filter!(x -> x[1] != "Quick Example", examples_nav)
quick_md = joinpath(build_path, "electrical_rc.md")
isfile(quick_md) && rm(quick_md)

# Group the examples in the sidebar, in this order; a file not listed goes to "Other"
example_sections = [
    "Systems" => ["mechanical_msd_1dof", "mechanical_msd_2dof", "mechanical_msd_ndof",
                  "mechanical_duffing", "electrical_dc_motor"],
    "Features" => ["two_ports", "signals", "adapters", "elements_oneport",
                   "equations_latex", "drawing"],
    "Applications" => ["control_dc_motor_pi", "robotics_motor_pendulum", "electrical_rlc",
                       "biomedical_windkessel", "hydraulic_coupled_tanks",
                       "chemical_reaction_network"],
]
stem(entry) = splitext(basename(entry[2]))[1]
listed = reduce(vcat, last.(example_sections))
examples_nav = [
    [title => [e for name in names for e in examples_nav if stem(e) == name]
     for (title, names) in example_sections];
    filter(!isempty ∘ last, ["Other" => filter(e -> stem(e) ∉ listed, examples_nav)])
]

# Generate documentation using Documenter
makedocs(
    modules=[BondGraphToolkit],
    # format=Documenter.LaTeX(platform="none"),
    format=Documenter.HTML(prettyurls=get(ENV, "CI", nothing) == "true", sidebar_sitename=false, assets=["assets/favicon.ico"], edit_link="main"),
    sitename="BondGraphToolkit.jl",
    # doctest=true,
    clean=true,
    pages=[
        "Home" => "index.md",
        "Getting started" => "getting_started_auto.md",
        "User Guide" => [
            "Concepts and conventions" => "man/user_guide.md",
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
    pagesonly=true,
    # draft=false,
    # remotes=nothing,
    # workdir = joinpath(@__DIR__, "build/assets")
    # checkbounds=false,  # Ignore broken links
)

# Cleaning
getting_started_auto = joinpath(@__DIR__, "src", "getting_started_auto.md")
isfile(getting_started_auto) && rm(getting_started_auto)

# Deploy documentation to GitHub Pages
deploydocs(repo="github.com/felipecastrotc/BondGraph.jl.git", devbranch="main")