import Base: +, -, show
import ModelingToolkit: equations, parameters, states
import ModelingToolkit: structural_simplify
import Symbolics.Latexify: latexify

# =============================================================================
# BgODESystem type

struct BgODESystem <: ModelingToolkit.AbstractODESystem
# struct BgODESystem
    ode::ODESystem
    sgn::Num
    subsys::Array{ODESystem, 1}
    couple::Num
end

function BgODESystem(bg::BgODESystem, couple)
    BgODESystem(bg.ode, bg.sgn, bg.subsys, couple)
end

# function BgODESystem(ode::ODESystem, sgn, subsys::Array{ODESystem, 1})
#     BgODESystem(ode, sgn, subsys, :missing, ode.name)
# end

# function BgODESystem(ode::ODESystem, sgn, subsys::Array{ODESystem, 1}, con::Symbol)
#     BgODESystem(ode, sgn, subsys, con, ode.name)
# end

# function BgODESystem(ode::ODESystem, con::Symbol)
#     BgODESystem(ode, 1, ODESystem[], con, ode.name)
# end

# Set basic functions for the BgODESystem

function Base.show(io::IO, mime::MIME"text/plain", sys::BgODESystem)
    show(io, mime, sys.ode)
end

function Base.getproperty(object::BgODESystem, item::Symbol)
    if item == :sgn || item == :subsys || item == :con || item == :ode || item == :name
        getfield(object, item)
    else
        getproperty(getproperty(object, :ode), item)
    end
end

function equations(sys::BgODESystem)
    equations(sys.ode)
end

function parameters(sys::BgODESystem)
    parameters(sys.ode)
end

function states(sys::BgODESystem)
    states(sys.ode)
end

function latexify(sys::BgODESystem)
    latexify(sys.ode)
end

function structural_simplify(sys::BgODESystem)
    structural_simplify(sys.ode)
end

function reducedobs(sys::BgODESystem)
    reducedobs(sys.ode)
end

function reducedobs(sys::BgODESystem)
    reducedobs(sys.ode)
end

function Base.nameof(item::BgODESystem)
    getproperty(item, :name)
end


-(sys::ODESystem) = BgODESystem(sys, sgn = -1)
+(sys::ODESystem) = BgODESystem(sys, sgn = 1)

-(sys::BgODESystem) = BgODESystem(sys.ode, -1*sys.sgn, sys.subsys)
+(sys::BgODESystem) = BgODESystem(sys.ode, 1*sys.sgn, sys.subsys)

