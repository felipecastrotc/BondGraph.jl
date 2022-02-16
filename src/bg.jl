@variables t
D = Differential(t)
It = Integral(t in DomainSets.ClosedInterval(0, t))

# After updating the packages I am unable to create an empty equation ODESystem
# I am replacing the  Power() connector with a simple function to return the
# states e and f
# @connector function Power(; name, effort = 0.0, flow = 0.0)
#     sts = @variables e(t) = effort f(t) = flow
#     ODESystem(Equation[], t, sts, []; name = name)
# end

function Power(; effort = 0.0, flow = 0.0)
    sts = @variables e(t) = effort f(t) = flow
    sts
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
        e, f = Power()
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
        sys = ODESystem(eqs, t, [e, f], []; name = name)
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
        e, f = Power()
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
        sys = ODESystem(eqs, t, [e, f], []; name = name)
    else
        sys = ODESystem(eqs, t, [], []; name = name)
    end

    compose(sys, rmsgnODE.(ps)...)
end

function mGY(subsys...; name, g = 1.0)

    # Get connections
    c = subsys isa ODESystem ? [subsys, nothing] : collect(subsys)

    # Remove nothing from c array
    pos = .!isnothing.(c)
    c = c[pos]
    @assert !isempty(c)

    # If only one subsys is passed it automatically generates an open  
    # connection
    if sum(pos) == 1
        e, f = Power()
        # Set variables according to the position
        if pos[1]
            e₁, f₁ = e, f
            e₂, f₂ = c[1].e, c[1].f
        else
            e₂, f₂ = e, f
            e₁, f₁ = c[1].e, c[1].f
        end
    else
        @assert length(c) == 2
        e₁, f₁ = c[1].e, c[1].f
        e₂, f₂ = c[2].e, c[2].f
    end

    # Gyrator equation
    eqs = [
        e₂ ~ g * f₁,
        e₁ ~ g * f₂,
    ]

    # Check if it is a modulated gyrator
    if isvariable(g)
        sts, ps = [], [g]
    elseif istree(unwrap(g))
        vars = collect(Set(ModelingToolkit.get_variables(r)))
        sts = filter(x -> ~isindependent(Num(x)), vars)
        ps = filter(x -> isindependent(Num(x)), vars)
    else
        sts, ps = [], []
    end

    if (@isdefined e) | (@isdefined f)
        push!(sts, e, f)
    end

    sys = ODESystem(eqs, t, sts, ps; name = name)

    compose(sys, c...)
end


function mTF(subsys...; name, r = 1.0)

    # Get connections
    c = subsys isa ODESystem ? [subsys, nothing] : collect(subsys)

    # Remove nothing from c array
    pos = .!isnothing.(c)
    c = c[pos]
    @assert !isempty(c)

    # If only one subsys is passed it automatically generates an open  
    # connection
    if sum(pos) == 1
        e, f = Power()
        # Set variables according to the position
        if pos[1]
            e₁, f₁ = e, f
            e₂, f₂ = c[1].e, c[1].f
        else
            e₂, f₂ = e, f
            e₁, f₁ = c[1].e, c[1].f
        end
    else
        @assert length(c) == 2
        e₁, f₁ = c[1].e, c[1].f
        e₂, f₂ = c[2].e, c[2].f
    end

    # Transformer equation
    eqs = [
        f₁ * r ~ f₂,
        e₂ * r ~ e₁,
    ]

    # Check if it is a modulated transformer
    if isvariable(r)
        sts, ps = [], [r]
    elseif istree(unwrap(r))
        vars = collect(Set(ModelingToolkit.get_variables(r)))
        sts = filter(x -> ~isindependent(Num(x)), vars)
        ps = filter(x -> isindependent(Num(x)), vars)
    else
        sts, ps = [], []
    end

    if (@isdefined e) | (@isdefined f)
        push!(sts, e, f)
    end

    sys = ODESystem(eqs, t, sts, ps; name = name)

    compose(sys, c...)
end

# =============================================================================
# Elements

function Mass(; name, m = 1.0, u = 0.0)
    e, f = Power(flow = u)
    ps = @parameters I = m
    @variables p(t)

    eqs = [
        # It(e) ~ f*I
        # D(p) ~ e,
        # p ~ I * f,
        # f ~ p/I,
        D(f) ~ e / I
    ]
    # ODESystem(eqs, t, [e, f, p], ps; name = name)
    ODESystem(eqs, t, [e, f], ps; name = name)
    # ODESystem(eqs, t, [], ps; name = name)
    # extend(ODESystem(eqs, t, [p], ps; name = name), power)
end

function Spring(; name, k = 10, x = 0.0)
    e, f = Power()
    @variables q(t) = x
    ps = @parameters C = 1 / k

    eqs = [
        e ~ q / C
        D(q) ~ f
        # D(e) ~ f/C
    ]
    ODESystem(eqs, t, [e, f, q], ps; name = name)
end

function Spring3(; name, k = 10, x = 0.0)
    e, f = Power()
    @variables q(t) = x
    ps = @parameters C = 1 / k

    eqs = [
        e ~ q^3 / C
        D(q) ~ f
        # D(e) ~ f/C
    ]
    ODESystem(eqs, t, [e, f, q], ps; name = name)
end

function Damper(; name, c = 10)
    e, f = Power()

    ps = @parameters R = c
    eqs = [
        e ~ f * R
    ]
    ODESystem(eqs, t, [e, f], ps; name = name)
end
