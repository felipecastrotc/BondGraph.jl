using ModelingToolkit, DifferentialEquations, LinearAlgebra
using Symbolics

include("lib_mtk.jl")
include("lib_utils.jl")

@variables t
D = Differential(t)

@connector function Power(; name, effort = 0.0, flow = 0.0)
    sts = @variables e(t) = effort f(t) = flow
    ODESystem(Equation[], t, sts, []; name = name)
end

# =============================================================================
# Sources

function Se(expr; name)
    sts = @variables e(t) = 0.0
    eqs = [sts[1] ~ expr]
    ODESystem(eqs, t, sts, []; name = name)
end

function Sf(expr; name)
    sts = @variables f(t) = 0.0
    eqs = [sts[1] ~ expr]
    ODESystem(eqs, t, sts, []; name = name)
end

function Dq(; name, x = 0.0)
    @variables f(t) = 0.0
    @variables q(t) = 0.0

    eqs = [
        D(q) ~ f
    ]
    ODESystem(eqs, t, [q, f], []; name = name)
end

# =============================================================================
# Ports

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
    # Remove empty equations
    if couple
        filter!(x -> filterexpr(x, ignore = [e, f]), eqs)
    end
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
    # Remove empty equations
    if couple
        filter!(x -> filterexpr(x, ignore = [e, f]), eqs)
    end
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
# =============================================================================
# Elements

function Mass(; name, m = 1.0, u = 0.0)
    @named power = Power(flow = u)
    @unpack e, f = power
    ps = @parameters I = m
    eqs = [
        D(f) ~ e / I
    ]
    extend(ODESystem(eqs, t, [], ps; name = name), power)
end

function Spring(; name, k = 10, x = 0.0)
    @named power = Power()
    @unpack e, f = power

    @variables q(t) = x
    ps = @parameters C = 1 / k

    eqs = [
        e ~ q / C
        D(q) ~ f
        # D(e) ~ f/C
    ]
    extend(ODESystem(eqs, t, [q], ps; name = name), power)
end

function Spring3(; name, k = 10, x = 0.0)
    @named power = Power()
    @unpack e, f = power

    @variables q(t) = x
    ps = @parameters C = 1 / k

    eqs = [
        e ~ q^3 / C
        D(q) ~ f
        # D(e) ~ f/C
    ]
    extend(ODESystem(eqs, t, [q], ps; name = name), power)
end

function Damper(; name, c = 10)
    @named power = Power()
    @unpack e, f = power

    ps = @parameters R = c
    eqs = [
        e ~ f * R
    ]
    extend(ODESystem(eqs, t, [], ps; name = name), power)
end
