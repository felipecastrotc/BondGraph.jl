"""
Main module for [`BongGraph.jl`](@ref) -- a bond graph modeling toolkit based on [`ModelingToolkit.jl`](https://docs.sciml.ai/ModelingToolkit/stable/) Julia package.

# Exports

- Power
- Se
- Sf
- Dq
- Junction1
- Junction0
- mGY
- mTF
- Mass
- Spring
- Spring3
- Damper
- reducedobs
- isindependent
- generate_graph
- simplifysys

# Imports

$(IMPORTS)

# License

The [`LICENSE`](@ref) abbreviation can be used in the same way for the `LICENSE.md` file.

"""
module BondGraph

using Reexport
using DocStringExtensions

@reexport using ModelingToolkit
# using Symbolics, DifferentialEquations, LinearAlgebra
using DomainSets

export Power,
    Se,
    Sf,
    Dq,
    Junction1,
    Junction0,
    mGY,
    mTF,
    Mass,
    Spring,
    Spring3,
    Damper,
    reducedobs,
    isindependent,
    generate_graph,
    simplifysys

# export ODESystem2D, SgnODESystem2D, Junction12D, Junction02D, mTF2Dtrans, M2Dtrans, get2Dsystem, Mass2D, Spring2D, Damper2D, equations, states, parameters

export renamevars, GenericDamper, latexify, structural_simplify, addsubsys


include("mtk.jl")
include("types.jl")
include("utils/bg_utils.jl")
include("utils/connection_utils.jl")
include("utils/graph.jl")
include("utils/system_utils.jl")
include("bg.jl")
# include("mbg.jl")


end # module