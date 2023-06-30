using ModelingToolkit
using GraphRecipes
using Plots
using LinearAlgebra

# ============================================================================
# 1st trial
# Generate connection set
connectionsets = ModelingToolkit.ConnectionSet[]
sys = ModelingToolkit.generate_connection_set!(connectionsets, mdl)

csets = deepcopy(connectionsets)
mcsets = ModelingToolkit.ConnectionSet[]
ale = Dict{String,Vector{String}}()    # Adjacency list for effort variable
alf = Dict{String,Vector{String}}()    # Adjacency list for effort variable
str2con = Dict{String,ModelingToolkit.ConnectionElement}()

vtype = get_connection_type(csets[1].set[1].v)

# Generate adjacency list from csets
for cset in csets
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

ale
alf

alef = filter(x -> length(x[2]) > 1, ale)
alef = filter(x -> length(x[2]) > 1, alf)

al = ale
al = alf

alef
# Generate adjacency matrix from the adjacency list
k2idx = Dict(k => i for (i, (k, vs)) in enumerate(al))
am = zeros(Int, length(al), length(al))

for (k, vs) in alef
    for v in vs
        if v in al[k]
            am[k2idx[k], k2idx[v]] = 1
        end
    end
end
am

idx2k
am
idx2k = [i => k for (i, (k, vs)) in enumerate(al)]
nm = [p[2] for p in idx2k]
# Plot connection graph
graphplot(am, names = nm, nodeshape = :rect, size = (600, 700), method = :stress)

# Get only one direction graph
dam = tril(am)'
graphplot(dam, names = nm, nodeshape = :rect, size = (600, 700), method = :stress)

# Generate an adjacency list from the directed graph adjacency matrix
ald = Dict{String,Vector{String}}()
for i = 1:size(am, 1)
    k = idx2k[i][2]
    # Initialize adjacency list key
    idx = get(ald, k, nothing)
    if idx === nothing
        ald[k] = String[]
    end
    # Add connections to the adjacency list
    for j = 1:size(dam, 2)
        if dam[i, j] > 0
            push!(ald[k], idx2k[j][2])
        end
    end
    if length(ald[k]) == 0
        delete!(ald, k)
    end
end
ald

mcsets = ModelingToolkit.ConnectionSet[]
for (k, v) in ald
    vsys = [str2con[j] for j in v]
    push!(mcsets, ModelingToolkit.ConnectionSet(vcat([str2con[k]], vsys)))
end

mcsets

sys, csets = generate_connection_set(sys)
csets = mcsets

tol = 1e-10
debug = false

ceqs, instream_csets =
    ModelingToolkit.generate_connection_equations_and_stream_connections(csets)
_sys = ModelingToolkit.expand_instream(instream_csets, sys; debug = debug, tol = tol)
sys = flatten(sys, true)
ModelingToolkit.@set! sys.eqs = [equations(_sys); ceqs]

equations(sys)
ceqs

equations(alias_elimination(sys))

# ============================================================================
# 2nd trial

# Generate connection set
connectionsets = ModelingToolkit.ConnectionSet[]
sys = ModelingToolkit.generate_connection_set!(connectionsets, mdl)

csets = deepcopy(connectionsets)
mcsets = ModelingToolkit.ConnectionSet[]
al = Dict{String,Vector{String}}()    # Adjacency list for effort variable
str2con = Dict{String,ModelingToolkit.ConnectionElement}()

vtype = get_connection_type(csets[1].set[1].v)

# Generate name to systems from csets
for cset in csets
    n = Vector{String}()
    for ele in cset.set
        # @unpack sys, v, isouter = e
        k = String(Symbol(ModelingToolkit.namespaced_var(ele)))
        # push!(n, k)
        # e = deepcopy(e)
        # ModelingToolkit.@set! e.isouter = false
        str2con[k] = ele
    end
end

# filter out flow
filterflow = true
filterstr = "f(t)"
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

# Generate adjacency list from csets
for cset in csets
    h = ""
    for (i, e) in enumerate(cset.set)
        @unpack sys, v, isouter = e
        k = String(nameof(sys)) * "₊" * String(Symbol(v))
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

# Visualize
idx2k = Dict(i => k for (i, k) in enumerate(keys(str2con)))
nm = [idx2k[i] for i = 1:length(idx2k)]
# Plot connection graph
graphplot(am, names = nm, nodeshape = :rect, size = (700, 700), method = :stress)

function get_var(idx, idx2k, str2con)
    return ModelingToolkit.namespaced_var(str2con[idx2k[idx]])
end

function add_idx(var, idx)
    v = deepcopy(var)
    vstr = String(ModelingToolkit.getname(v))
    newname = Symbol(vstr * string(idx))
    return ModelingToolkit.rename(v, newname)
end


in_vec = sum(am, dims = 1)[:]
out_vec = sum(am, dims = 2)[:]

idx_con = (1:size(am, 1))
in_con = idx_con[in_vec.>0]
out_con = idx_con[out_vec.>0]

eqs = Equation[]
for i in in_con
    v = get_var(i, idx2k, str2con)
    # Get variables leaving the node
    vout = Num[]
    if out_vec[i] == 1
        push!(vout, -get_var(i, idx2k, str2con))
    elseif out_vec[i] > 1
        for j = 1:out_vec[i]
            push!(vout, -add_idx(v, j))
        end
    end
    # Get variables being added to the node
    vin = Num[]
    for (j, k) in enumerate(am[:, i])
        # k is a value checking if there is a connection or not
        if k > 0
            # Check if the input has multiple outputs
            msk = am[j, :] .> 0
            chk = sum(msk) > 1 ? findfirst(sort(idx_con[msk]) .== i) : nothing
            if !isnothing(chk)
                # Create the variable for the i node
                push!(vin, add_idx(get_var(j, idx2k, str2con), chk))
            else
                push!(vin, get_var(j, idx2k, str2con))
            end
        end
    end
    # Apply the junction type
    jtype = get_bg_junction(v)[1]
    vtype = get_bg_junction(v)[2]
    if jtype === j0
        if vtype === bgeffort
            push!(eqs, equalityeqs(vcat(vout, vin))...)
        elseif vtype === bgflow
            push!(eqs, sumvar(vcat(vout, vin)))
        end
    elseif jtype === j1
        if vtype === bgeffort
            push!(eqs, sumvar(vcat(vout, vin)))
        elseif vtype === bgflow
            push!(eqs, equalityeqs(vcat(vout, vin))...)

        end
    end
end

eqs

ceqs = deepcopy(eqs)

tol = 1e-10
debug = false

sys, csets = ModelingToolkit.generate_connection_set(mdl)
ceqs, instream_csets =
    ModelingToolkit.generate_connection_equations_and_stream_connections(csets)
ceqs
_sys = ModelingToolkit.expand_instream(instream_csets, sys; debug = debug, tol = tol)
sys = flatten(sys, true)
ModelingToolkit.@set! sys.eqs = [equations(_sys); ceqs]

states(sys)
equations(sys)

equations(alias_elimination(sys))

@named syss = reducedobs(structural_simplify(sys))

# ============================================================================
# Tiding up the 2nd trial

# Generate connection set
connectionsets = ModelingToolkit.ConnectionSet[]
sys = ModelingToolkit.generate_connection_set!(connectionsets, mdl)

csets = deepcopy(connectionsets)
mcsets = ModelingToolkit.ConnectionSet[]
str2con = Dict{String,ModelingToolkit.ConnectionElement}()

# Generate name to systems from csets
for cset in csets
    n = Vector{String}()
    for ele in cset.set
        k = String(Symbol(ModelingToolkit.namespaced_var(ele)))
        str2con[k] = ele
    end
end

# filter out flow
filterflow = false
filterstr = "f(t)"
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

# Visualize
idx2k = Dict(i => k for (i, k) in enumerate(keys(str2con)))
nm = [idx2k[i] for i = 1:length(idx2k)]
# Plot connection graph
graphplot(am, names = nm, nodeshape = :rect, size = (700, 700), method = :stress)

function get_var(idx, idx2k, str2con)
    return ModelingToolkit.namespaced_var(str2con[idx2k[idx]])
end

function add_idx(var, idx)
    v = deepcopy(var)
    vstr = String(ModelingToolkit.getname(v))
    newname = Symbol(vstr * string(idx))
    return ModelingToolkit.rename(v, newname)
end

sm = deepcopy(am)

# Identify the one port elements to apply the signal
cord = findall(am .> 0)
filter!(x -> get_bg_junction(get_var(x[2], idx2k, str2con))[1] === op, cord)
# Convert the connections and apply the signal matrix sm
for c in cord
    flipc = CartesianIndex(c[2], c[1])
    am[c], am[flipc] = 0, 1
    sm[c], sm[flipc] = 0, -1
end

in_vec = sum(am, dims = 1)[:]
out_vec = sum(am, dims = 2)[:]

idx_con = (1:size(am, 1))
in_con = idx_con[in_vec.>0]

eqs = Equation[]
for i in in_con
    v = get_var(i, idx2k, str2con)
    # Get variables leaving the node
    vout = Num[]
    if out_vec[i] == 1
        push!(vout, -get_var(i, idx2k, str2con))
    elseif out_vec[i] > 1
        for j = 1:out_vec[i]
            push!(vout, -add_idx(v, j))
        end
    end
    # Get variables being added to the node
    vin = Num[]
    # Iterate over the connections
    for j in idx_con[am[:, i].>0]
        # Check if the input has multiple outputs
        msk = am[j, :] .> 0
        chk = sum(msk) > 1 ? findfirst(sort(idx_con[msk]) .== i) : nothing
        jtype = get_bg_junction(get_var(j, idx2k, str2con))[1]
        if !isnothing(chk)
            # Create the variable for the i node
            push!(vin, add_idx(get_var(j, idx2k, str2con), chk))
        elseif (jtype === tpgy) || (jtype === tptf)
            push!(vin, -get_var(j, idx2k, str2con))
        elseif (jtype === op)
            # Apply the signal matrix
            sgn = sm[j, i]
            push!(vin, sgn * get_var(j, idx2k, str2con))
        else
            push!(vin, get_var(j, idx2k, str2con))
        end
    end

    # Apply the junction type
    jtype = get_bg_junction(v)[1]
    vtype = get_bg_junction(v)[2]
    if (jtype === j0 && vtype === bgeffort) || (jtype === j1 && vtype === bgflow)
        push!(eqs, equalityeqs(vcat(vout, vin))...)
    elseif (jtype === j0 && vtype === bgflow) || (jtype === j1 && vtype === bgeffort)
        push!(eqs, sumvar(vcat(vout, vin)))
    elseif (jtype === tpgy) || (jtype === tptf)
        lvin, lvout = length(vin), length(vout)
        if (lvin <= 1) && (lvout <= 1) && (lvout + lvout == 1)
            throw(
                DomainError(
                    "The " *
                    string(jtype) *
                    " type only allows one connection in and one connection out",
                ),
            )
        end
        if lvin == 1
            push!(eqs, equalityeqs(vcat(vin, v))...)
        elseif lvout == 1
            push!(eqs, equalityeqs(vcat(vout, v))...)
        end
    end
end

eqs

eqs
tol = 1e-10
debug = false

sys, csets = ModelingToolkit.generate_connection_set(mdl)
ceqs, instream_csets =
    ModelingToolkit.generate_connection_equations_and_stream_connections(csets)
equations(sys)
_sys = ModelingToolkit.expand_instream(instream_csets, sys; debug = debug, tol = tol)
equations(_sys)
sys = flatten(sys, true)
ceqs = deepcopy(eqs)
ModelingToolkit.@set! sys.eqs = [equations(_sys); ceqs]

states(sys)
equations(sys)

equations(alias_elimination(sys))

@named syss = reducedobs(structural_simplify(sys))

equations(syss)
