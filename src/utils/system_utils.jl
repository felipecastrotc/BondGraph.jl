import ModelingToolkit: structural_simplify, expand_connections


"""
    simplifysys(sys::ODESystem; name)

Simplify an ODE system by expanding connections, performing structural simplification and reducing the observed variables.

# Arguments:
- sys::ODESystem: The ODE system to simplify.
- name: Optional name for the simplified system.

# Returns:
- The simplified ODE system.
"""
function simplifysys(sys::ODESystem; name)
    exp_sys = expand_connections(sys)
    sys = structural_simplify(exp_sys)
    reducedobs(sys; name = name)
end


"""
    reducedobs(sys::ODESystem; name)

Create a reduced ODESystem by substituting observed variables in the equations.

# Arguments:
- sys::ODESystem: The original ODESystem object.
- name: Optional name for the reduced ODESystem.

# Returns:
A reduced ODESystem object with observed variables substituted.

# Example:
```julia-repl
julia> sys = ODESystem([eq1, eq2], t)
julia> reduced_sys = reducedobs(sys; name = "Reduced System")
```
"""
function reducedobs(sys::ODESystem; name)
    eqs = map(eq -> substitute_observed(eq, sys), equations(sys))
    # Create and return a new ODESystem with substituted equations
    ODESystem(eqs; name = name)
end


"""
    substitute_dict(O, expr, sub)

Substitute the variables in `O` with their corresponding values from the dictionary `sub` in the expression `expr`.

# Arguments
- `O`: The variable or expression to be substituted.
- `expr`: The expression in which the substitution should be performed.
- `sub`: A dictionary mapping variables to their corresponding values.

# Returns
The expression `expr` with the variables in `O` substituted using the values from `sub`.
"""
function substitute_dict(O, expr, sub)
    O = unwrap(O)  # Unwrap the variable if it is a ModelingToolkit variable

    if O in keys(sub)
        # If the variable is in the substitution dictionary, substitute it with the corresponding value
        return substitute(expr, Dict(O => sub[O]))
    elseif istree(O)
        # If the variable is a tree, recursively substitute its arguments
        A = arguments(O)
        for a in A
            expr = substitute_dict(a, expr, sub)
        end
    end

    return expr
end

"""
    substitute_observed(eq::Equation, sys::ODESystem)

Substitute observed variables in the equation with their assigned values.

# Arguments:
- eq::Equation: The equation to be processed.
- sys::ODESystem: The ODESystem object containing the observed variables.

# Returns:
An equation with observed variables substituted by their assigned values.

# Example:
```julia-repl
julia> eq = @eq x'(t) = a * x(t) + b * y(t)
julia> sys = ODESystem([eq], t)
julia> substituted_eq = substitute_observed(eq, sys)
```
"""
function substitute_observed(eq::Equation, sys::ODESystem)

    O = eq.rhs

    # Create a dictionary mapping observed variables to their assigned values
    sub = Dict(o.lhs => o.rhs for o in observed(sys))

    # Substitute observed variables in the equation with their assigned values
    expr = substitute_dict(O, O, sub)

    # Repeat the substitution until no further changes are made
    while !isequal(O, expr)
        O = expr
        expr = substitute_dict(O, O, sub)
    end

    # Return the equation with substituted observed variables
    return eq.lhs ~ expr
end
