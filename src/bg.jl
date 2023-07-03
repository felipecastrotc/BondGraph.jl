@variables t
D = Differential(t)
It = Integral(t in DomainSets.ClosedInterval(0, t))

# Base

"""
    Power(; name, effort = 0.0, flow = 0.0, type = op)

Creates a power connector for the bond graph modeling.

Arguments:
- name: Name of the power connector.
- effort: Initial effort value. Default is 0.0.
- flow: Initial flow value. Default is 0.0.
- type: Type of the power connector. Default is op (one-port).

Returns:
An ODESystem representing the power connector.

Example:
Power(name = "P1", effort = 1.0, flow = 2.0, type = tpgy)
"""
@connector function Power(; name, effort = 0.0, flow = 0.0, type = op)
    # Define the symbolic variables for effort and flow
    sts = @variables e(t) = effort [connect = bg] f(t) = flow [connect = bg]

    # Set the bond graph metadata for the symbolic variables
    sts = set_bg_metadata.(sts, [[type, bgeffort], [type, bgflow]])

    # Create an ODESystem with empty equations, time variable t, the symbolic variables sts, and optional name
    ODESystem(Equation[], t, sts, []; name = name)
end

# Sources

"""
    Se(expr; name)

Creates a source of effort element in the bond graph model.
TODO: Improve

Arguments:
- expr: Expression representing the source value.
- name: Name of the source element.

Returns:
An ODESystem representing the source element.

Example:
@named s1 = Se(2.0)
"""
function Se(expr; name)
    # Create a Power connector with type op (one-port)
    @named power = Power(type = op)

    # Define the equation between the power.e variable and the source expression
    eqs = [power.e ~ expr]

    sts = []
    if isvariable(expr)
        # If the expression is a variable, add it to the list of states
        ps = [expr]
    elseif istree(unwrap(expr))
        # If the expression is a tree, extract the variables
        vars = collect(Set(ModelingToolkit.get_variables(expr)))
        sts = vcat(sts, filter(x -> ~isindependent(Num(x)), vars))
        ps = filter(x -> isindependent(Num(x)), vars)
    else
        ps = []
    end

    # Compose the ODESystem with the equations, time variable t, states sts, parameter states ps, and name
    compose(ODESystem(eqs, t, sts, ps; name = name), power)
end

"""
    Sf(expr; name)

Creates a source of flow in the bond graph model. CARAIOOOOO
TODO: Improvess

# Arguments:
- expr: Expression representing the flow source value.
- name: Name of the flow source element.

# Returns:
An ODESystem representing the flow source element.

# Example:
@named s2 = Sf(3.5)
"""
function Sf(expr; name)
    # Create a Power connector with type op (one-port)
    @named power = Power(type=op)

    # Define the equation between the power.f variable and the flow source expression
    eqs = [power.f ~ expr]

    sts = []
    if isvariable(expr)
        # If the expression is a variable, add it to the list of states
        ps = [expr]
    elseif istree(unwrap(expr))
        # If the expression is a tree, extract the variables
        vars = collect(Set(ModelingToolkit.get_variables(expr)))
        sts = vcat(sts, filter(x -> ~isindependent(Num(x)), vars))
        ps = filter(x -> isindependent(Num(x)), vars)
    else
        ps = []
    end

    # Compose the ODESystem with the equations, time variable t, states sts, parameter states ps, and name
    compose(ODESystem(eqs, t, sts, ps; name=name), power)
end

function Dq(; name, x = 0.0)
    @variables f(t) = 0.0
    @variables q(t) = 0.0

    eqs = [D(q) ~ f]
    ODESystem(eqs, t, [q, f], []; name = name)
end

# Ports

function Junction0(ps...; name, couple = true)

    @named power = Power(type = j0)
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

function Junction1(ps...; name, couple = true)

    @named power = Power(type = j1)
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

# TODO:ISSUE -> the compose only works when the gyrator systems is the first
function mGY(subsys...; name, g = 1.0, coneqs = [])

    # Generate the in and out connection
    @named pin = Power(type = tpgy)
    @named pout = Power(type = tpgy)

    # Alias for simpler view about the mGY port
    e₁, f₁ = pin.e, pin.f
    e₂, f₂ = pout.e, pout.f

    # Gyrator equation
    eqs = [e₂ ~ -g * f₁, e₁ ~ g * f₂]
    # TODO: check why adding the - on f₁ solves the signal issue on DC motor
    # TODO: I already tried to change the signal on the port type tpgy at 
    #       adjmtx2eqs

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

    sys = compose(ODESystem(eqs, t, sts, ps; name = name), pin, pout)

    gen_tp_con!(coneqs, sys, subsys)

    return sys
end

# TODO:ISSUE -> the compose only works when the transform systems is the first
function mTF(subsys...; name, r = 1.0, coneqs = [])

    # Generate the in and out connection
    @named pin = Power(type = tptf)
    @named pout = Power(type = tptf)

    # Alias for simpler view about the mTF port
    e₁, f₁ = pin.e, pin.f
    e₂, f₂ = pout.e, pout.f

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

    sys = compose(ODESystem(eqs, t, sts, ps; name = name), pin, pout)

    gen_tp_con!(coneqs, sys, subsys)

    return sys
end


# Elements

function Mass(; name, m = 1.0, u = 0.0)
    @named power = Power(flow = u)
    ps = @parameters I = m

    eqs = [D(power.f) ~ power.e / I]
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

function Damper(; name, c = 10, u = 1.0)
    @named power = Power(flow = u)

    ps = @parameters R = c
    eqs = [power.e ~ power.f * R]

    compose(ODESystem(eqs, t, [], ps; name = name), power)
end

function GenericDamper(expr; name)
    # TODO: check if this function is working

    @named power = Power(flow = u)

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
