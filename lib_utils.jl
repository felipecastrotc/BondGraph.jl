using ModelingToolkit, Symbolics


function extract_vars(O)
    O = unwrap(O)
    if istree(O)
        reduce(vcat, [extract_vars(o) for o in arguments(O) if !(o isa Number)])
    elseif isvariable(O)
        O
    end
end

# =============================================================================
# Functions to replace variables on equations by the observed variables

function substitute_dict(O, expr, sub)
    O = unwrap(O)

    if O in keys(sub)

        return substitute(expr, Dict(O => sub[O]))

    elseif istree(O)
        A = arguments(O)
        for a in A
            expr = substitute_observed(a, expr, sub)
        end
    end

    return expr
end

function substitute_observed(eq::Equation, sys::ODESystem)

    O = eq.rhs
    sub = Dict(o.lhs => o.rhs for o in observed(sys))

    expr = substitute_dict(O, O, sub)
    while !isequal(O, expr)
        O = expr
        expr = substitute_dict(O, O, sub)
    end

    return eq.lhs ~ expr
end

function reducedobs(sys::ODESystem; name)
    eqs = map(eq -> substitute_observed(eq, sys), equations(sys))
    ODESystem(eqs; name = name)
end
