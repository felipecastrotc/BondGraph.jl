using ModelingToolkit, DifferentialEquations, LinearAlgebra
using Symbolics

include("lib_mtk.jl")
include("lib_utils.jl")

@variables t
D = Differential(t)

@connector function Power(;name, effort=0.0, flow=0.0)
    sts = @variables e(t)=effort f(t)=flow
    ODESystem(Equation[], t, sts, []; name=name)
end

# =============================================================================
# Ports

function Junction1(ps...; name, subsys=[], se=[0.0], couple=true)
    # Check subsys type
    sys = subsys isa ODESystem ? [subsys] : subsys

    # Get connections
    con = collect(Base.Flatten([ps, sys]))
    if couple
        @named power = Power()
        @unpack e, f = power
    
        # Σ efforts
        eqs = [e ~ sum(c -> c.e, con) + sum(c -> c, se)]
        # f₁ = f₂ = f₃
        if length(con) > 1
            for i in 1:(length(con)-1)
                push!(eqs, con[i].f ~ con[i+1].f)
            end
        end
        push!(eqs, con[1].f ~ f)
        # Build subsystem
        ex = extend(ODESystem(eqs, t, [], []; name = name), power)
        if length(ps) > 0
            compose(ex, ps...)
        else
            ex
        end
    else
        # Σ efforts
        eqs = [e ~ sum(c -> c.e, con) + sum(c -> c, force)]
        # f₁ = f₂ = f₃
        if length(con) > 1
            for i in 1:(length(con)-1)
                push!(eqs, con[i].f ~ con[i+1].f)
            end
        end
        # Build system
        ex = ODESystem(eqs, t, [], []; name = name)
        if length(ps) > 0
            compose(ex, ps...)
        else
            ex
        end
    end
end

function Junction0(ps...; name, subsys=[], se=[],couple=true)
    # Check subsys type
    sys = subsys isa ODESystem ? [subsys] : subsys

    # Get connections
    con = collect(Base.Flatten([ps, sys]))
    if couple
        @named power = Power()
        @unpack e, f = power
    
        # Σ flows
        eqs = [f ~ sum(c -> c.f, con)]
        # e₁ = e₂ = e₃
        if length(con) > 1
            for i in 1:(length(con) - 1)
                push!(eqs, con[i].e ~ con[i + 1].e)
            end
        end
        push!(eqs, con[1].e ~ e)
        if length(se) > 0
            push!(eqs, con[1].e ~ se[1])
            for i in 2:length(se)
                push!(eqs, se[i-1] ~ se[i])
            end
        end
        # Build subsystem
        ex = extend(ODESystem(eqs, t, [], []; name=name), power)
        if length(ps) > 0
            compose(ex, ps...)
        else
            ex
        end
    else
        # Σ flows
        eqs = [0 ~ sum(p -> p.f, con)]
        # e₁ = e₂ = e₃
        if length(con) > 1
            for i in 1:(length(con) - 1)
                push!(eqs, con[i].e ~ con[i + 1].e)
            end
        end
        if length(se) > 0
            push!(eqs, con[1].e ~ se[1])
            for i in 2:length(se)
                push!(eqs, se[i-1] ~ se[i])
            end
        end
        # Build system
        ex = ODESystem(eqs, t, [], []; name=name)
        if length(ps) > 0
            compose(ex, ps...)
        else
            ex
        end
    end
end

# =============================================================================
# Elements

function Mass(; name, m = 1.0, u = 0.)
    @named power = Power(flow=u)
    @unpack e, f = power
    ps = @parameters I=m
    eqs = [
            D(f) ~ e/I
        ]
    extend(ODESystem(eqs, t, [], ps; name=name), power)
end

function Spring(; name, k = 10, x = 0.)
    @named power = Power()
    @unpack e, f = power

    @variables q(t)=x
    ps = @parameters C=1/k
    
    eqs = [
            e ~ q/C
            D(q) ~ f
            # D(e) ~ f/C
        ]
    extend(ODESystem(eqs, t, [q], ps; name=name), power)
end

function Spring3(; name, k = 10, x = 0.)
    @named power = Power()
    @unpack e, f = power

    @variables q(t)=x
    ps = @parameters C=1/k
    
    eqs = [
            e ~ q^3/C
            D(q) ~ f
            # D(e) ~ f/C
        ]
    extend(ODESystem(eqs, t, [q], ps; name=name), power)
end

function Damper(; name, c = 10)
    @named power = Power()
    @unpack e, f = power

    ps = @parameters R=c
    eqs = [
            e ~ f*R
        ]
    extend(ODESystem(eqs, t, [], ps; name=name), power)
end
