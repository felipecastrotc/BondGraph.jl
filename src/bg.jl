@variables t
D = Differential(t)
It = Integral(t in DomainSets.ClosedInterval(0, t))

# Base

"""
    Power(; name, effort = 0.0, flow = 0.0, type = op)

Creates a power connector for the bond graph modeling.

# Arguments:
- name: Name of the power connector.
- effort: Initial effort value. Default is 0.0.
- flow: Initial flow value. Default is 0.0.
- type: Type of the power connector. Default is op (one-port).

# Returns:
An ODESystem representing the power connector.

# Example:
```julia-repl
julia> Power(name = "P1", effort = 1.0, flow = 2.0, type = tpgy)
```

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

# Arguments:
- expr: Expression representing the source value.
- name: Name of the source element.

# Returns:
An ODESystem representing the source element.

# Example:
```julia-repl
julia> @named s1 = Se(2.0)
```
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
        # Remove t from parameters
        ps = filter(x -> Symbolics.wrap(x) !== t, ps)
    else
        ps = []
    end

    # Compose the ODESystem with the equations, time variable t, states sts, parameter states ps, and name
    compose(ODESystem(eqs, t, sts, ps; name = name), power)
end

"""
    Sf(expr; name)

Creates a source of flow in the bond graph model
TODO: Improvess

# Arguments:
- `expr`: Expression representing the flow source value.
- name: Name of the flow source element.

# Returns:
An ODESystem representing the flow source element.

# Example:
```julia-repl
julia>  @named s2 = Sf(3.5)
```
"""
function Sf(expr; name)
    # Create a Power connector with type op (one-port)
    @named power = Power(type = op)

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
        # Remove t from parameters
        ps = filter(x -> Symbolics.wrap(x) !== t, ps)
    else
        ps = []
    end

    # Compose the ODESystem with the equations, time variable t, states sts, parameter states ps, and name
    compose(ODESystem(eqs, t, sts, ps; name = name), power)
end

function Dq(; name, x = 0.0)
    @variables f(t) = 0.0
    @variables q(t) = 0.0

    eqs = [D(q) ~ f]
    ODESystem(eqs, t, [q, f], []; name = name)
end

# Ports
"""
    Junction0(ps...; name="", couple=true)

Create a zero-junction element in the bond graph model.

# Arguments
- `ps...`: One or more elements to be connected to the zero-junction.
- `name::String`: The name of the zero-junction element (default: "").

# Returns
An `ODESystem` representing the zero-junction element in the bond graph model.

# Remarks
- The `ps` argument accepts one or more elements to be connected to the zero-junction.
- The function automatically determines the flow directions and subsystems based on the signs of the connections.

# Examples
```julia-repl
# Create a zero-junction connecting multiple elements
julia> junction = Junction0(damper1, damper2, spring1, spring2, name="junction")
```
"""
function Junction0(ps...; name = "")
    # Create a new Power element for the zero-junction
    @named power = Power(type = j0)

    # Flatten the connections to handle variable number of arguments
    ps = collect(Base.Flatten([ps]))

    # Get signs and subsystems
    subsys, signs = flatinput(ps)

    eqs = Equation[]

    # Connect the elements based on the signs
    for (s, sign) in zip(subsys, signs)
        if sign == 1
            push!(eqs, connect(s.power, power))
        elseif sign == -1
            push!(eqs, connect(power, s.power))
        end
    end

    # Build the subsystem with the zero-junction equations
    sys = ODESystem(eqs, t, [], [], name = name)
    compose(sys, power, subsys...)
end

"""
    Junction1(ps...; name="", couple=true)

Create a one-junction element in the bond graph model.

# Arguments
- `ps...`: One or more elements to be connected to the one-junction.
- `name::String`: The name of the one-junction element (default: "").

# Returns
An `ODESystem` representing the one-junction element in the bond graph model.

# Remarks
- The `ps` argument accepts one or more elements to be connected to the one-junction.
- The function automatically determines the flow directions and subsystems based on the signs of the connections.

# Examples
```julia-repl
# Create a one-junction connecting multiple elements
julia> junction = Junction1(damper1, damper2, spring1, spring2, name="junction")
```
"""
function Junction1(ps...; name = "")
    # Create a new Power element for the one-junction
    @named power = Power(type = j1)

    # Flatten the connections to handle variable number of arguments
    ps = collect(Base.Flatten([ps]))

    # Get signs and subsystems
    subsys, signs = flatinput(ps)

    eqs = Equation[]

    # Connect the elements based on the signs
    for (s, sign) in zip(subsys, signs)
        if sign == 1
            push!(eqs, connect(s.power, power))
        elseif sign == -1
            push!(eqs, connect(power, s.power))
        end
    end

    # Build the subsystem with the one-junction equations
    sys = ODESystem(eqs, t, [], [], name = name)
    compose(sys, power, subsys...)
end

"""
    mGY(subsys...; name="", g=1.0, coneqs=[])

Create a gyrator element in the bond graph model.

# Arguments
- `subsys...`: One or more subsystems to be connected to the gyrator.
- `name::String`: The name of the gyrator element (default: "").
- `g::Number|Symbol`: A value or a expression representing the gyrator relationship (default: 1.0).
- `coneqs::Vector{Equation}`: Additional output connection equations of the modulated gyrator (default: []).

# Returns
An `ODESystem` representing the gyrator element in the bond graph model.

# Examples
```julia-repl
# Create a gyrator with a variable value
@variables c
julia> gyrator = mGY(subsys1, subsys2, g=c, name="gyrator")
```
"""
function mGY(subsys...; name = "", g = 1.0, coneqs = nothing)
    # TODO:ISSUE -> the compose only works when the gyrator systems is the first

    # Generate the input and output connections
    @named pin = Power(type = tpgy)
    @named pout = Power(type = tpgy)

    # Alias for simpler view of the mGY port
    e₁, f₁ = pin.e, pin.f
    e₂, f₂ = pout.e, pout.f

    # Gyrator equation
    eqs = [e₂ ~ g * f₁, e₁ ~ g * f₂]
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

    # Create the subsystem with the modulated gyrator equations
    sys = compose(ODESystem(eqs, t, sts, ps; name = name), pin, pout)

    # Generate the additional connection equations
    con = []
    gen_tp_con!(con, sys, subsys)

    if isnothing(coneqs)
        return sys, con
    else
        push!(coneqs, con...)
        return sys
    end
end

"""
    mTF(subsys...; name="", r=1.0, coneqs=[])

Create a transformer element in the bond graph model.

# Arguments
- `subsys...`: One or more subsystems to be connected to the transformer.
- `name::String`: The name of the transformer element (default: "").
- `r::Number|Symbol`: A value or a expression representing the transformer relationship (default: 1.0).
- `coneqs::Vector{Equation}`: Additional output connection equations of the modulated gyrator (default: []).

# Returns
An `ODESystem` representing the transformer element in the bond graph model.

# Examples
```julia-repl
# Create a transformer with a variable value
@variables c
julia> transformer = mTF(subsys1, subsys2, g=c, name="transformer")
```
"""
function mTF(subsys...; name="", r=1.0, coneqs=nothing)
    # TODO:ISSUE -> the compose only works when the transform systems is the first
    # Generate the input and output connections
    @named pin = Power(type=tptf)
    @named pout = Power(type=tptf)

    # Alias for simpler view of the mTF port
    e₁, f₁ = pin.e, pin.f
    e₂, f₂ = pout.e, pout.f

    # Transformer equations
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

    # Create the subsystem with the modulated transformer equations
    sys = compose(ODESystem(eqs, t, sts, ps; name=name), pin, pout)

    # Generate the additional connection equations
    con = []
    gen_tp_con!(con, sys, subsys)

    if isnothing(coneqs)
        return sys, con
    else
        push!(coneqs, con...)
        return sys
    end
end

# Elements

"""
    Mass(; name="", m=1.0, u=0.0)

Create a mass element in the bond graph model.

# Arguments
- `name::String`: The name of the mass element (default: "").
- `m::Float64`: The mass value (default: 1.0).
- `u::Float64`: The initial value for the flow variable (default: 0.0).

# Returns
An `ODESystem` representing the mass element in the bond graph model.

# Examples
```julia-repl
julia> mass = Mass(name="my_mass", m=2.0, u=1.0)
```
"""
function Mass(; name = "", m = 1.0, u = 0.0)
    @named power = Power(flow = u)
    ps = @parameters I = m

    eqs = [D(power.f) ~ power.e / I]
    compose(ODESystem(eqs, t, [], ps; name = name), power)
end


"""
    Spring(; name="", k=10, x=0.0)

Create a spring element in the bond graph model.

# Arguments
- `name::String`: The name of the spring element (default: "").
- `k::Float64`: The spring stiffness (default: 10).
- `x::Float64`: The initial displacement (default: 0.0).

# Returns
An `ODESystem` representing the spring element in the bond graph model.

# Examples
```julia-repl
julia> spring = Spring(name="my_spring", k=5.0, x=0.1)
```
"""
function Spring(; name = "", k = 10, x = 0.0)
    @named power = Power()

    @variables q(t) = x
    ps = @parameters C = 1 / k

    eqs = [power.e ~ q / C, D(q) ~ power.f]
    compose(ODESystem(eqs, t, [q], ps; name = name), power)
end

"""
    Damper(; name="", c=10, u=1.0)

Create a damper element in the bond graph model.

# Arguments
- `name::String`: The name of the damper element (default: "").
- `c::Float64`: The damping coefficient (default: 10).
- `u::Float64`: The initial flow variable value (default: 1.0).

# Returns
An `ODESystem` representing the damper element in the bond graph model.

# Examples
```julia-repl
julia> damper = Damper(name="my_damper", c=5.0, u=0.1)
```
"""
function Damper(; name = "", c = 10, u = 1.0)
    @named power = Power(flow = u)

    ps = @parameters R = c
    eqs = [power.e ~ power.f * R]

    compose(ODESystem(eqs, t, [], ps; name = name), power)
end


"""
    GenericDamper(expr; name="")

Create a generic damper element in the bond graph model.

# Arguments
- `expr`: The expression defining the behavior of the damper element.
- `name::String`: The name of the damper element (default: "").

# Returns
An `ODESystem` representing the generic damper element in the bond graph model.

# Remarks
- The `expr` can be any valid expression defining the behavior of the damper element.
- The function automatically determines the state variables and parameters based on the `expr`.
"""
function GenericDamper(expr; name = "")
    # TODO: check if this function is working

    @named power = Power()

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
