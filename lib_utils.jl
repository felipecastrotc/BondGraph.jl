using ModelingToolkit, Symbolics
import Base: +, -

# Deprecated -> use ModelingToolkit.get_variables
function extract_vars(O)
    O = unwrap(O)
    if istree(O)
        reduce(vcat, [extract_vars(o) for o in arguments(O) if !(o isa Number)])
    elseif isvariable(O)
        O
    end
end

# =============================================================================
# Support junction functions

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
    if length(C) == 0
        0.0
    else
        sum(c -> c.sgn * getproperty(c.ode, sym), C)
    end
end

function equalityeqs(con::Vector{SgnODESystem}, var::String; couple = false, sgn = 1)
    sym = Symbol(var)
    C = filter(c -> hasproperty(c.ode, sym), con)

    if length(C) > 0
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

    else
        Vector{Equation}(undef, 0)
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
            expr = substitute_dict(a, expr, sub)
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


# =============================================================================
# multibond graph 


