module BondGraph

using Reexport

@reexport using ModelingToolkit
# using Symbolics, DifferentialEquations, LinearAlgebra
using DomainSets

export Power, Se, Sf, Dq, Junction1, Junction0, mGY, mTF, Mass, Spring, Spring3, Damper, reducedobs, isindependent

# export ODESystem2D, SgnODESystem2D, Junction12D, Junction02D, mTF2Dtrans, M2Dtrans, get2Dsystem, Mass2D, Spring2D, Damper2D, equations, states, parameters

export renamevars, GenericDamper, latexify, structural_simplify, addsubsys


include("mtk.jl")
include("types.jl")
include("utils.jl")
include("bg.jl")
include("mbg.jl")


end # module
