using GraphRecipes

import ModelingToolkit: ConnectionSet, ConnectionElement
import ModelingToolkit: namespaced_var, get_connection_type, getname
import ModelingToolkit: rename, generate_connection_set!

# =============================================================================
# Depecrate function

function directcon(con)
    filter(c -> !hasproperty(c, :power), con)
end

function haspower(con)
    filter(c -> hasproperty(c, :power), con)
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
# Function to handle bond graph connections
# =============================================================================
# Interface functions
function check_bg_con(connectionset)
    ele = namespaced_var(connectionset.set[1])
    return get_connection_type(ele) === bg
end

function get_bg_connection_set!(connectionsets)
    bgconnectionsets = filter(check_bg_con, connectionsets)
    filter!(x -> !check_bg_con(x), connectionsets)
    return bgconnectionsets
end

# =============================================================================
# BG functions
function equalityeqs(con)
    
    base = 0.0
    if length(con) > 1
        base += con[1]
    else
        return []
    end

    eqs = Equation[]
    for i in 2:length(con)
        push!(eqs, base ~ con[i])
    end

    return eqs
end

function sumvar(con)
    if length(con) > 0 
        return 0 ~ sum(con)
    else
        return []
    end
end

function flatinput(ps)

    subsys = ModelingToolkit.AbstractSystem[]
    signs = Int[]

    for (i, p) in enumerate(ps)
        if isa(p, ModelingToolkit.AbstractSystem)
            push!(subsys, p)
            push!(signs, 1)
        elseif isa(p, AbstractVector)
            if (length(p) == 2) && isa(p[1], Int) && isa(p[2], ModelingToolkit.AbstractSystem)
                push!(subsys, p[2])
                push!(signs, p[1])
            else
                throw(DomainError("The input number "*string(i)*" is not valid, a valid element has a two elements AbstractVector, where the first is an integer and the second is an AbstractSystem"))
            end
        else
            throw(DomainError(p, "The input number "*string(i)*" is not valid. The valid are: AbstractSystem or a two element AbstractVector."))
        end
    end

    return subsys, signs
end

# =============================================================================
# Util functions

function get_var(idx, idx2k, str2con)
    return namespaced_var(str2con[idx2k[idx]])
end

function add_idx(var, idx)
    v = deepcopy(var)
    vstr = String(getname(v))
    newname = Symbol(vstr*string(idx))
    return rename(v, newname)
end

# =============================================================================
# Algorithm functions

function csets2dict(csets)
    
    str2con = Dict{String, ConnectionElement}()

    # Generate name to systems from csets
    for cset in csets
        n = Vector{String}()
        for ele in cset.set
            k = String(Symbol(namespaced_var(ele)))
            str2con[k] = ele
        end
    end

    return str2con
end

function csets2adjmtx(csets, str2con; filterstr="f(t)", filterflow=false)

    # Filter str2con
    if filterflow
        filter!(d -> occursin(filterstr, d[1]), str2con)
    end

    # Generate adjacency matrix from csets
    k2idx = Dict(k => i for (i, k) in enumerate(keys(str2con)))

    am = zeros(Int, length(str2con), length(str2con))
    for cset in csets
        h = 0
        for (i, e) in enumerate(cset.set)
            @unpack sys, v, isouter = e
            k = String(nameof(sys)) * "₊" * String(Symbol(v))
            if String(Symbol(v)) != filterstr && filterflow
                continue
            end
            if i == 1
                h = k2idx[k]
            else
                am[h, k2idx[k]] = 1
            end
        end
    end

    return am
end

# TODO: add support adding signs it has issues on the algorithm
function adjmtx2eqs(am, str2con)
    
    idx2k = Dict(i => k for (i, k) in enumerate(keys(str2con)))

    in_vec = sum(am, dims=1)[:]
    out_vec = sum(am, dims=2)[:]

    idx_con = (1:size(am, 1))
    in_con = idx_con[in_vec .> 0]

    eqs = Equation[]
    for i in in_con
        v = get_var(i, idx2k, str2con)
        # Get variables leaving the node
        vout = Num[]
        if out_vec[i] == 1
            push!(vout, -get_var(i, idx2k, str2con))
        elseif out_vec[i] > 1
            for j in 1:out_vec[i]
                push!(vout, -add_idx(v, j))
            end
        end
        # Get variables being added to the node
        vin = Num[]
        # Iterate over the connections
        for j in idx_con[am[:, i] .> 0]
            # Check if the input has multiple outputs
            msk= am[j, :] .> 0
            chk = sum(msk) > 1 ? findfirst(sort(idx_con[msk]) .== i) : nothing
            if !isnothing(chk)
                # Create the variable for the i node
                push!(vin, add_idx(get_var(j, idx2k, str2con), chk))
            else
                push!(vin, get_var(j, idx2k, str2con))
            end
        end
        # Apply the junction type
        jtype = get_bg_junction(v)[1]
        vtype = get_bg_junction(v)[2]
        if (jtype == j0  && vtype === bgeffort) || (jtype === j1 && vtype === bgflow)
            push!(eqs, equalityeqs(vcat(vout, vin))...)
        elseif (jtype == j0  && vtype === bgflow) || (jtype === j1 && vtype === bgeffort)
            push!(eqs, sumvar(vcat(vout, vin)))
        elseif (jtype === tpgy) || (jtype === tptf)
            lvin, lvout = length(vin), length(vout)
            if (lvin <= 1) && (lvout <= 1) && (lvout + lvout == 1)
                throw(DomainError("The "*string(jtype)*" type only allows one connection in and one connection out"))
            end
            if lvin == 1
                push!(eqs, equalityeqs(vcat(vin, v))...)
            elseif lvout == 1
                push!(eqs, equalityeqs(vcat(-vout, v))...)
            end
        end
    end

    return eqs
end

# =============================================================================
# Main functions
function generate_bg_eqs!(connectionsets)

    bgconnectionsets = get_bg_connection_set!(connectionsets)
    str2con = csets2dict(bgconnectionsets)
    am = csets2adjmtx(bgconnectionsets, str2con)
    eqs = adjmtx2eqs(am, str2con)
    return eqs
end

# Support functions
function generate_graph(mdl, var=:e)

    connectionsets = ConnectionSet[]
    sys = generate_connection_set!(connectionsets, mdl)

    bgconnectionsets = get_bg_connection_set!(connectionsets)
    str2con = csets2dict(bgconnectionsets)

    varin = var == :e ? "e(t)" : "f(t)"
    am = csets2adjmtx(bgconnectionsets, str2con, filterstr=varin, filterflow=true)
    
    idx2k = Dict(i => k for (i, k) in enumerate(keys(str2con)))
    nm = [idx2k[i] for i in 1:length(idx2k)]

    graphplot(am, names=nm, nodeshape=:rect, size=(600, 700), method=:stress)

end