# Extending ModelingToolkit
import Symbolics: Symbolic
# equations functions
import SymbolicUtils: FnType
import ModelingToolkit: expand_connections, generate_connection_set
import ModelingToolkit: generate_connection_equations_and_stream_connections
import ModelingToolkit: expand_instream, @set!, AbstractSystem, ConnectionSet
import ModelingToolkit: generate_connection_set!
import ModelingToolkit: is_domain_connector, domain_defaults, get_defaults

import ModelingToolkit: namespace_equation, namespace_expr
import ModelingToolkit: independent_variables, unwrap
import ModelingToolkit: isvariable, renamespace
import ModelingToolkit: istree, arguments
import ModelingToolkit: symtype, operation
import ModelingToolkit: similarterm, getmetadata
import ModelingToolkit: getname, rename, setmetadata
import ModelingToolkit: namespace_equations, equations
import ModelingToolkit: single_named_expr, split_assign, check_name
import ModelingToolkit: gensym, _named


# Connection functions

function expand_connections(
    sys::AbstractSystem,
    find = nothing,
    replace = nothing;
    debug = false,
    tol = 1e-10,
)
    sys, (csets, domain_csets) = generate_connection_set(sys, find, replace)
    bgeqs = generate_bg_eqs!(domain_csets)
    ceqs, instream_csets = generate_connection_equations_and_stream_connections(csets)
    _sys = expand_instream(instream_csets, sys; debug = debug, tol = tol)
    sys = flatten(sys, true)
    @set! sys.eqs = [equations(_sys); ceqs; bgeqs]
    d_defs = domain_defaults(sys, domain_csets)
    @set! sys.defaults = merge(get_defaults(sys), d_defs)
end

function generate_connection_set(sys::AbstractSystem, find = nothing, replace = nothing)
    connectionsets = ConnectionSet[]
    domain_csets = ConnectionSet[]
    sys = generate_connection_set!(connectionsets, domain_csets, sys, find, replace)

    # Add the bgconnectionsets
    bgconnectionsets = get_bg_connection_set!(connectionsets)

    csets = merge(connectionsets)
    domain_csets = merge([csets; domain_csets], true)

    sys, (csets, vcat(domain_csets, bgconnectionsets))
end


# Named replace function

function split_assign(expr)
    # specific case for var1, var2 = function(args...)
    if expr isa Expr && expr.head === :tuple && expr.args[2].args[2].head === :call
        name = :($(expr.args[1]), $(expr.args[2].args[1]))
        call = expr.args[2].args[2]
    elseif !(expr isa Expr && expr.head === :(=) && expr.args[2].head === :call)
        throw(ArgumentError("expression should be of the form `sys = foo(a, b)`"))
    else
        name, call = expr.args
    end
    return name, call
end

function single_named_expr(expr)
    name, call = split_assign(expr)
    if Meta.isexpr(name, :ref)
        name, idxs = name.args
        check_name(name)
        var = gensym(name)
        ex = quote
            $var = $(_named(name, call))
            $name = map(i -> $rename($var, Symbol($(Meta.quot(name)), :_, i)), $idxs)
        end
        ex
    else
        if typeof(name) === Expr
            varname, sysname = name, name.args[1]
        else
            varname, sysname = name, name
        end
        check_name(sysname)
        :($varname = $(_named(sysname, call)))
    end
end


# TODO: Check for legacy functions

function namespace_equation_couple(eq::Equation, sys, couple)
    _lhs = namespace_expr_couple(eq.lhs, sys, couple)
    _rhs = namespace_expr_couple(eq.rhs, sys, couple)
    _lhs ~ _rhs
end

function namespace_expr_couple(O, sys, couple) where {T}
    ivs = independent_variables(sys)
    O = unwrap(O)
    if any(isequal(O), ivs)
        return O
    elseif isvariable(O)
        renamespace(sys, O, couple)
    elseif istree(O)
        renamed = map(a -> namespace_expr_couple(a, sys, couple), arguments(O))
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
                if isa(couple, Symbol)
                    couple = [couple]
                end

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
    map(eq -> namespace_equation_couple(eq, sys, couple), eqs)
end

function equations(sys::Array{ODESystem,1})

    # Systems to be coupled
    couple = [string(s.name) for s in sys]

    eqs = reduce(vcat, namespace_equations.(sys, Ref(couple)); init = Equation[])

    return eqs
end


# Rename variables

function namespace_equation_ren(
    eq::Equation,
    sys,
    old::Union{Num,String},
    new::Union{Num,Symbol},
)
    _lhs = namespace_expr_ren(eq.lhs, sys, old, new)
    _rhs = namespace_expr_ren(eq.rhs, sys, old, new)
    _lhs ~ _rhs
end

function namespace_expr_ren(
    O,
    sys,
    old::Union{Num,String},
    new::Union{Num,Symbol},
) where {T}
    ivs = independent_variables(sys)
    O = unwrap(O)
    if any(isequal(O), ivs)
        return O
    elseif isvariable(O)
        renamespace(sys, O, old, new)
    elseif istree(O)
        renamed = map(a -> namespace_expr_ren(a, sys, old, new), arguments(O))
        if symtype(operation(O)) <: FnType
            renamespace(sys, O, old, new)
        else
            similarterm(O, operation(O), renamed)
        end
    elseif O isa Array
        map(Base.Fix2(namespace_expr_ren, sys), O)
    else
        O
    end
end

function renamespace(sys, x, old::Num, new::Num)
    x = unwrap(x)
    if isequal(x, old)
        unwrap(new)
    else
        x
    end
end

function renamespace(sys, x, old::String, new::Symbol)
    x = unwrap(x)
    if isequal(string(x), old)
        rename(x, new)
    else
        x
    end
end


function renamevars(sys::ODESystem, pairs::Dict)
    eqs = equations(sys)

    if isempty(eqs)
        return sys
    else
        for (old, new) in pairs
            eqs = map(eq -> namespace_equation_ren(eq, sys, old, new), eqs)
        end
    end

    return ODESystem(eqs; name = sys.name)
end
