using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Symbolics, Symbolics.Latexify

include("lib_bg.jl")

# =============================================================================
# Function Junction1


function isdefnothing(c::ODESystem, var::String)
    sym = Symbol(var)
    st = getproperty(c, sym, namespace = false)
    isnothing(ModelingToolkit.get_defaults(c)[st])
end


function sumvar(con::Vector{ODESystem}, var::String)
    sym = Symbol(var)
    C = filter(c -> !isdefnothing(c, sym), con)
    sum(c -> getproperty(c, sym), C)
end


function equalityeqs(con::Vector{ODESystem}, var::String; couple=false)
    sym = Symbol(var)
    C = filter(c -> !isdefnothing(c, sym), con)

    eqs = Vector{Equation}(undef, length(C) - 1)
    for i in 1:(length(C)-1)
        eqs[i] = getproperty(C[i], sym) ~ getproperty(C[i + 1], sym)
    end

    if couple
        V = getproperty(C[end], sym)
        v = getproperty(C[end], sym, namespace = false)
        push!(eqs, V ~ v)
    end

    eqs
end


function Junction1(ps...; name, subsys = [], couple = true)
    # Check subsys type
    sys = subsys isa ODESystem ? [subsys] : subsys

    # Get connections
    ps = Base.Flatten([ps])
    con = collect(Base.Flatten([ps, sys]))

    if couple
        @named power = Power()
        @unpack e, f = power
    else
        e, f = 0.0, nothing
    end

    # Σ efforts
    eqs = [e ~ sumvar(con, "e")]
    # f₁ = f₂ = f₃
    eqs = vcat(eqs, equalityeqs(con, "f", couple=couple))

    # Build subsystem
    if couple
        sys = extend(ODESystem(eqs, t, [], []; name = name), power)
    else
        sys = ODESystem(eqs, t, [], []; name = name)
    end

    compose(sys, ps...)
end

function Junction0(ps...; name, subsys = [], couple = true)
    # Check subsys type
    sys = subsys isa ODESystem ? [subsys] : subsys

    # Get connections
    ps = Base.Flatten([ps])
    con = collect(Base.Flatten([ps, sys]))

    if couple
        @named power = Power()
        @unpack e, f = power
    else
        e, f = nothing, 0.0
    end

    # Σ flows
    eqs = [f ~ sumvar(con, "f")]
    # e₁ = e₂ = e₃
    eqs = vcat(eqs, equalityeqs(con, "e", couple = couple))

    # Build subsystem
    if couple
        sys = extend(ODESystem(eqs, t, [], []; name = name), power)
    else
        sys = ODESystem(eqs, t, [], []; name = name)
    end

    compose(sys, ps...)
end