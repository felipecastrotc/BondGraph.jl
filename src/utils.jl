# =============================================================================
# Deprecate

# Consider deprecate
function filterexpr(expr::Equation; ignore::Vector{Num} = Num[])

    expr_rhs = Set(ModelingToolkit.get_variables(expr.rhs))
    expr_lhs = Set(ModelingToolkit.get_variables(expr.lhs))
    expr_var = union(expr_rhs, expr_lhs)

    vars = collect(setdiff(expr_var, Set(ignore)))

    length(vars) > 0
end

# Deprecated with the usage of hasproperty
function isdefnothing(c::ODESystem, sym::String)
    # sym = Symbol(var)
    st = getproperty(c, sym, namespace = false)
    isnothing(ModelingToolkit.get_defaults(c)[st])
end

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

function addsgn(sys)
    if typeof(sys) == ODESystem
        +sys
    elseif typeof(sys) == BgODESystem
        sys
    end
end

function rmsgn(sys)
    if typeof(sys) == ODESystem
        sys
    elseif typeof(sys) == BgODESystem
        sys.ode
    end
end

function getsign(sys)
    if hasproperty(sys, :sgn)
        return getproperty(sys, :sgn)
    else
        return 1
    end
end

function getsym(sys, sym::Symbol)
    if hasproperty(sys, :power)
        return getproperty(getproperty(sys, :power), sym)
    else
        return getproperty(sys,sym)
    end
end

function getoneport(con)
    filter(c -> !hasproperty(c, :power), con)
end

function getmultiport(con)
    filter(c -> hasproperty(c, :power), con)
end

function sumvar(con, sym::Symbol)
    sum(c -> getsign(c)*getsym(c, sym), con; init=0.0)
end

function sumcoupling(vars, signs::Vector)
    if length(signs) > 1
        return sum([sgn*c for (c, sgn) in zip(vars, signs)])
    elseif length(signs) == 1
        return signs[1]*vars
    else
        return 0.0
    end
end

function equalityeqs(con, sym::Symbol, signs, vars)
    
    vars = isa(vars, AbstractArray) ? vars : [vars]

    start_con = 1
    start_vars = 1
    
    base = 0.0
    if length(con) > 0
        base += getsym(con[1], sym) 
        start_con = 2
    elseif length(sign) > 0
        base += sign[1]*vars[1]
        start_vars = 2
    else
        return []
    end

    eqs = Equation[]

    for i in start_con:length(con)
        push!(eqs, base ~ getsym(con[i], sym))
    end

    for i in start_vars:length(signs)
        push!(eqs, base ~ signs[i]*vars[i])
    end

    return eqs
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

function isindependent(var::Num)
    if isvariable(var)
        !istree(unwrap(var))
    else
        true
    end
end

# =============================================================================
# Function to handle subsystem

function addsubsys(sys, pss)
    
    ps = pss isa BgODESystem ? [pss] : collect(pss)

    if BgODESystem in typeof.(ps)
        
        eqs, names = Equation[], Symbol[]

        rename = Dict()

        filter!(x -> x isa BgODESystem, ps);

        for p in ps
            for s in p.subsys

                if !(s.name in names)
                    push!(names, s.name)
                    eqs = vcat(eqs, equations(s))
                end

                pname = string(p.name)*"₊"
                sname = string(s.name)*"₊"
                pname = pname*sname

                for st in s.states
                    str_st = string(st)
                    sym_st = split(str_st, "(")[1]
                    rename[pname*str_st] = Symbol(sname*sym_st)
                end
            end
        end

        sys = renamevars(sys, rename)
        return sys
        # Remove duplicate equations
        # eqs = collect(Set(vcat(equations(sys), eqs)))
    #    return ODESystem(eqs; name = sys.name)
    else
        return sys
    end
end

# =============================================================================
# Function to handle bond graph connections Deprecated
using LinearAlgebra

function csets2adjlist_bckp(csets)

    ale = Dict{String, Vector{String}}()    # Adjacency list for effort variable
    alf = Dict{String, Vector{String}}()    # Adjacency list for effort variable
    str2con = Dict{String, ModelingToolkit.ConnectionElement}()

    # Generate adjacency list from csets
    for cset in csets
        if get_connection_type(cset.set[1].v) === bg
            n = Vector{String}()
            for e in cset.set
                @unpack sys, v, isouter = e
                k = String(nameof(sys)) * "₊" * String(Symbol(v))
                push!(n, k)
                e = deepcopy(e)
                ModelingToolkit.@set! e.isouter = false
                str2con[k] = e
            end
            N = Set(n)
            for e in cset.set
                @unpack sys, v, isouter = e
                al = String(Symbol(v)) == "f(t)" ? alf : ale
                k = String(nameof(sys)) * "₊" * String(Symbol(v))
                idx = get(al, k, nothing)
                if idx === nothing
                    al[k] = String[]
                end
                push!(al[k], collect(setdiff(N, [k]))...)
            end
        end
    end

    return ale, alf, str2con
end

function mtx2list(am, idx2k::Vector)
    # Generate an adjacency list from the adjacency matrix

    al = Dict{String, Vector{String}}()

    for i in 1:size(am, 1)
        k = idx2k[i][2]
        # Initialize adjacency list key
        idx = get(al, k, nothing)
        if idx === nothing
            al[k] = String[]
        end
        # Add connections to the adjacency list
        for j in 1:size(am, 2)
            if am[i, j] > 0
                push!(al[k], idx2k[j][2])
            end
        end
        if length(al[k]) == 0
            delete!(al, k)
        end
    end

    return al
end

function mtx2list(am, al::Dict)
    idx2k = [i => k for (i, (k, vs)) in enumerate(al)]
    return mtx2list(am, idx2k)
end

function am2dam(am)
    # Convert undirected graph to directed graph by
    # transforming the adjacency matrix in a upper 
    # triangular matrix
    # println(am)
    return tril(am)'
end

# =============================================================================
# Function to handle bond graph connections

function name2sys(csets)

    str2con = Dict{String, ModelingToolkit.ConnectionElement}()

    # Generate name to systems from csets
    for cset in csets
        if get_connection_type(cset.set[1].v) === bg
            n = Vector{String}()
            for e in cset.set
                @unpack sys, v, isouter = e
                k = String(nameof(sys)) * "₊" * String(Symbol(v))
                push!(n, k)
                e = deepcopy(e)
                ModelingToolkit.@set! e.isouter = false
                str2con[k] = e
            end
        end
    end

    return str2con
end

function csets2adjlist(csets, var="")

    al = Dict{String, Vector{String}}()    # Adjacency list for effort variable

    # Generate adjacency list from csets
    for cset in csets
        h = ""
        for (i, e) in enumerate(cset.set)
            @unpack sys, v, isouter = e
            k = String(nameof(sys)) * "₊" * String(Symbol(v))
            if (var == String(Symbol(v))) || length(var) == 0
                if i == 1
                    h = k
                    idx = get(al, k, nothing)
                    if idx === nothing
                        al[k] = String[]
                    end
                else
                    push!(al[h], k)
                end
            end
        end
    end

    return al

end

function csets2adjmtx(csets, str2con)

    # Generate adjacency matrix from csets
    k2idx = Dict(k => i for (i, k) in enumerate(keys(str2con)))

    am = zeros(Int, length(str2con), length(str2con))
    for cset in csets
        h = 0
        for (i, e) in enumerate(cset.set)
            @unpack sys, v, isouter = e
            k = String(nameof(sys)) * "₊" * String(Symbol(v))
            if i == 1
                h = k2idx[k]
            else
                am[h, k2idx[k]] = 1
            end
        end
    end
    return am
end

function csets2adjmtx(csets)
    str2con = name2sys(csets)
    return csets2adjmtx(csets, str2con)
end

function list2mtx(al)
    # Generate adjacency matrix from the adjacency list
    # k2idx = Dict(k => i for (i, (k, vs)) in enumerate(al))
    K = collect(Set(vcat([k for (k, vs) in al], [vs for (k, vs) in al]...)))
    k2idx = Dict(k => i for (i, k) in enumerate(K))

    am = zeros(Int, length(k2idx), length(k2idx))

    for (k, vs) in al
        for v in vs
            if v in al[k]
                am[k2idx[k], k2idx[v]] = 1
            end
        end
    end

    return am, k2idx
end

function adjlist2csets(al, str2con)
    mcsets = ModelingToolkit.ConnectionSet[]
    for (k, v) in al
        vsys = [str2con[j] for j in v]
        push!(mcsets, ModelingToolkit.ConnectionSet(vcat([str2con[k]], vsys)))
    end
    return mcsets
end

function generate_graph(mdl, var=:e)
    connectionsets = ModelingToolkit.ConnectionSet[]

    sys = ModelingToolkit.generate_connection_set!(connectionsets, mdl)

    varin = var == :e ? "e(t)" : "f(t)"
    al = csets2adjlist(connectionsets, varin);
    
    am, k2idx = list2mtx(al)
    
    idx2k = [v => k for (k, v) in k2idx]
    nm = [p[2] for p in idx2k]

    graphplot(am, names=nm, nodeshape=:rect, size=(600, 700), method=:stress)


end