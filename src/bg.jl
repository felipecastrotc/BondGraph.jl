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

@connector function Power(; name, effort=0.0, flow=0.0, type=op)
    sts = @variables e(t) = effort [connect = bg] f(t) = flow [connect = bg]
    sts = set_bg_metadata.(sts, [[type, bgeffort], [type, bgflow]])
    ODESystem(Equation[], t, sts, []; name=name)
end

# =============================================================================
# Sources

function Se(expr; name)

    @named power = Power(type=op)

    eqs = [power.e ~ expr]
    
    sts = []
    if isvariable(expr)
        ps = [expr]
    elseif istree(unwrap(expr))
        vars = collect(Set(ModelingToolkit.get_variables(g)))
        sts = vcat(sts, filter(x -> ~isindependent(Num(x)), vars))
        ps = filter(x -> isindependent(Num(x)), vars)
    else
        ps = []
    end

    compose(ODESystem(eqs, t, sts, ps; name = name), power)
end

function Sf(expr; name)

    @named power = Power(type=op)

    eqs = [power.f ~ expr]
    
    sts = []
    if isvariable(expr)
        ps = [expr]
    elseif istree(unwrap(expr))
        vars = collect(Set(ModelingToolkit.get_variables(g)))
        sts = vcat(sts, filter(x -> ~isindependent(Num(x)), vars))
        ps = filter(x -> isindependent(Num(x)), vars)
    else
        ps = []
    end

    compose(ODESystem(eqs, t, sts, ps; name = name), power)
end

function Dq(; name, x = 0.0)
    @variables f(t) = 0.0
    @variables q(t) = 0.0

    eqs = [D(q) ~ f]
    ODESystem(eqs, t, [q, f], []; name = name)
end

# =============================================================================
# Ports

function Junction0(ps...; name, couple=true)
    
    @named power = Power(type=j0)
    # Get connections
    ps = collect(Base.Flatten([ps]))

    # Get signs and subsystems 
    subsys, signs = flatinput(ps)

    eqs = Equation[]
    for (s, sign) in zip(subsys, signs)
        if sign == 1
            push!(eqs, connect(s.power, power))
        elseif sign == -1
            push!(eqs, connect(power, s.power))
        end
    end

    # Build subsystem
    sys = ODESystem(eqs, t, [], [], name = name)
    compose(sys, power, subsys...)
end

function Junction1(ps...; name, couple=true)

    @named power = Power(type=j1)
    # Get connections
    ps = collect(Base.Flatten([ps]))

    # Get signs and subsystems 
    subsys, signs = flatinput(ps)

    eqs = Equation[]
    for (s, sign) in zip(subsys, signs)
        if sign == 1
            push!(eqs, connect(s.power, power))
        elseif sign == -1
            push!(eqs, connect(power, s.power))
        end
    end

    # Build subsystem
    sys = ODESystem(eqs, t, [], [], name = name)
    compose(sys, power, subsys...)
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
    eqs = [e₂ ~ g * f₁, e₁ ~ g * f₂]

    # Check if it is a modulated gyrator
    if isvariable(g)
        sts, ps = [], [g]
    elseif istree(unwrap(g))
        vars = collect(Set(ModelingToolkit.get_variables(g)))
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
    eqs = [f₁ * r ~ f₂, e₂ * r ~ e₁]

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
    @named power = Power(flow=u)    
    ps = @parameters I = m

    eqs = [
        D(power.f) ~ power.e / I,
    ]
    compose(ODESystem(eqs, t, [], ps; name = name), power)
end

function Spring(; name, k = 10, x = 0.0)
    @named power = Power()

    @variables q(t) = x
    ps = @parameters C = 1 / k

    eqs = [
        power.e ~ q / C
        D(q) ~ power.f
    ]
    compose(ODESystem(eqs, t, [q], ps; name = name), power)
end

function Spring3(; name, k = 10, x = 0.0)
    # TODO: check if this function is working
    @named power = Power()

    @variables q(t) = x
    ps = @parameters C = 1 / k

    eqs = [
        power.e ~ q^3 / C
        D(q) ~ power.f
    ]
    
    compose(ODESystem(eqs, t, [q], ps; name = name), power)
end

function Damper(; name, c = 10, u=1.0)
    @named power = Power(flow=u)
    
    ps = @parameters R = c
    eqs = [power.e ~ power.f * R]

    compose(ODESystem(eqs, t, [], ps; name = name), power)
end

function GenericDamper(expr; name)
    # TODO: check if this function is working

    @named power = Power(flow=u)

    eqs = [power.e ~ expr]

    # Check the expression
    if isvariable(expr)
        sts, ps = [], [expr]
    elseif istree(unwrap(expr))
        vars = collect(Set(ModelingToolkit.get_variables(expr)))
        sts = filter(x -> ~isindependent(Num(x)), vars)
        ps = filter(x -> isindependent(Num(x)), vars)
    else
        sts, ps = [], []
    end

    compose(ODESystem(eqs, t, sts, ps; name = name), power)
end
