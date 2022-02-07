# Extending ModelingToolkit
import Symbolics: Symbolic
# equations functions
import SymbolicUtils: FnType
import ModelingToolkit: namespace_equation, namespace_expr
import ModelingToolkit: independent_variables, unwrap
import ModelingToolkit: isvariable, renamespace
import ModelingToolkit: istree, arguments
import ModelingToolkit: symtype, operation
import ModelingToolkit: similarterm, getmetadata
import ModelingToolkit: getname, rename, setmetadata
import ModelingToolkit: namespace_equations, equations


function namespace_equation(eq::Equation, sys, couple)
    _lhs = namespace_expr(eq.lhs, sys, couple)
    _rhs = namespace_expr(eq.rhs, sys, couple)
    _lhs ~ _rhs
end

function namespace_expr(O, sys, couple) where {T}
    ivs = independent_variables(sys)
    O = unwrap(O)
    if any(isequal(O), ivs)
        return O
    elseif isvariable(O)
        renamespace(sys, O, couple)
    elseif istree(O)
        renamed = map(a -> namespace_expr(a, sys, couple), arguments(O))
        if symtype(operation(O)) <: FnType
            renamespace(sys, O, couple)
        else
            similarterm(O, operation(O), renamed)
        end
    elseif O isa Array
        map(Base.Fix2(namespace_expr, sys), O)
    else
        O
    end
end

function renamespace(sys, x, couple)
    x = unwrap(x)
    if x isa Symbolic
        let scope = getmetadata(x, SymScope, LocalScope())
            if scope isa LocalScope
                sys_name = getname(sys)
                var_name = getname(x)

                chk = filter(x -> x ≠ string(sys_name), couple)
                lvl = split(string(var_name), "₊")

                if lvl[1] in chk
                    x
                    # var_name = Symbol(join(deleteat!(lvl, 1), "₊"))
                else
                    rename(x, renamespace(sys_name, var_name, couple))
                end
            elseif scope isa ParentScope
                setmetadata(x, SymScope, scope.parent)
            else # GlobalScope
                x
            end
        end
    else
        Symbol(getname(sys), :₊, x)
    end
end

function namespace_equations(sys::ODESystem, couple)
    eqs = equations(sys)
    isempty(eqs) && return Equation[]
    map(eq -> namespace_equation(eq, sys, couple), eqs)
end

function equations(sys::Array{ODESystem,1})

    # Systems to be coupled
    couple = [string(s.name) for s in sys]

    eqs = reduce(vcat, namespace_equations.(sys, Ref(couple));
        init = Equation[])

    return eqs
end