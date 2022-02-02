using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Symbolics, Symbolics.Latexify
import Base: +, -

include("lib_bg.jl")

# =============================================================================
# Function

struct SgnODESystem
    ode::ODESystem
    sgn::Num
end

-(sys::ODESystem) = SgnODESystem(sys, -1)
+(sys::ODESystem) = SgnODESystem(sys, 1)

function addsgnODE(sys)
    if typeof(sys) == ODESystem
        +sys
    elseif typeof(sys) == SgnODESystem
        sys
    end
end

function rmsgnODE(sys)
    if typeof(sys) == ODESystem
        sys
    elseif typeof(sys) == SgnODESystem
        sys.ode
    end
end

# Deprecated with the usage of hasproperty
function isdefnothing(c::ODESystem, sym::String)
    # sym = Symbol(var)
    st = getproperty(c, sym, namespace = false)
    isnothing(ModelingToolkit.get_defaults(c)[st])
end

function sumvar(con::Vector{SgnODESystem}, var::String)
    sym = Symbol(var)
    C = filter(c -> hasproperty(c.ode, sym), con)
    sum(c -> c.sgn * getproperty(c.ode, sym), C)
end

function equalityeqs(con::Vector{SgnODESystem}, var::String; couple = false, sgn = 1)
    sym = Symbol(var)
    C = filter(c -> hasproperty(c.ode, sym), con)

    eqs = Vector{Equation}(undef, length(C) - 1)
    for i in 1:(length(C)-1)
        # f₁ = C[i].sgn * getproperty(C[i].ode, sym)
        # f₂ = C[i+1].sgn * getproperty(C[i+1].ode, sym)
        f₁ = getproperty(C[i].ode, sym)
        f₂ = getproperty(C[i+1].ode, sym)
        eqs[i] = f₁ ~ f₂
    end

    if couple
        # f₁ = C[end].sgn * getproperty(C[end].ode, sym)
        f₁ = getproperty(C[end].ode, sym)
        f₂ = sgn * getproperty(C[end].ode, sym, namespace = false)
        push!(eqs, f₁ ~ f₂)
    end

    eqs
end

function Junction1(ps...; name, subsys = [], couple = true, sgn = 1)
    # Check subsys type
    sys = subsys isa ODESystem ? [subsys] : subsys

    # Get connections
    ps = Base.Flatten([ps])
    con = addsgnODE.(collect(Base.Flatten([ps, sys])))

    if couple
        @named power = Power()
        @unpack e, f = power
    else
        e, f = 0.0, nothing
    end

    # Σ efforts
    eqs = [0 ~ sumvar(con, "e") + sgn * e]
    # f₁ = f₂ = f₃
    eqs = vcat(eqs, equalityeqs(con, "f", couple = couple))
    # TODO: CHECK WHY DISABILITATING THE SIGN IN THE EQUALITY 
    # THE DC MOTOR MODEL MATCHES WITH THE LITERATURE
    # TODO: CHECK IF THE SIGN CONVENTION IS ONLY FOR THE SUM
    # Example: https://www.20sim.com/webhelp/modeling_tutorial_bond_graphs_frombondgraphtoequations.php
    # The example considers the convention only for the sum in the
    # equality it does not consider the arrow direction

    # Build subsystem
    if couple
        sys = extend(ODESystem(eqs, t, [], []; name = name), power)
    else
        sys = ODESystem(eqs, t, [], []; name = name)
    end


    compose(sys, rmsgnODE.(ps)...)
end

function Junction0(ps...; name, subsys = [], couple = true, sgn = 1)
    # Check subsys type
    sys = subsys isa ODESystem ? [subsys] : subsys

    # Get connections
    ps = Base.Flatten([ps])
    con = addsgnODE.(collect(Base.Flatten([ps, sys])))

    if couple
        @named power = Power()
        @unpack e, f = power
    else
        e, f = nothing, 0.0
    end

    # Σ flows
    eqs = [0 ~ sumvar(con, "f") + sgn * f]
    # e₁ = e₂ = e₃
    eqs = vcat(eqs, equalityeqs(con, "e", couple = couple, sgn = sgn))
    # TODO: CHECK WHY DISABILITATING THE SIGN IN THE EQUALITY 
    # THE DC MOTOR MODEL MATCHES WITH THE LITERATURE
    # TODO: CHECK IF THE SIGN CONVENTION IS ONLY FOR THE SUM
    # Example: https://www.20sim.com/webhelp/modeling_tutorial_bond_graphs_frombondgraphtoequations.php
    # The example considers the convention only for the sum in the
    # equality it does not consider the arrow direction

    # Build subsystem
    if couple
        sys = extend(ODESystem(eqs, t, [], []; name = name), power)
    else
        sys = ODESystem(eqs, t, [], []; name = name)
    end

    compose(sys, rmsgnODE.(ps)...)
end