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
ale = Dict{String, Vector{String}}()    # Adjacency list for effort variable
alf = Dict{String, Vector{String}}()    # Adjacency list for effort variable
str2con = Dict{String, ModelingToolkit.ConnectionElement}()

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
graphplot(am, names=nm, nodeshape=:rect, size=(600, 700), method=:stress)

# Get only one direction graph
dam = tril(am)'
graphplot(dam, names=nm, nodeshape=:rect, size=(600, 700), method=:stress)

# Generate an adjacency list from the directed graph adjacency matrix
ald = Dict{String, Vector{String}}()
for i in 1:size(am, 1)
    k = idx2k[i][2]
    # Initialize adjacency list key
    idx = get(ald, k, nothing)
    if idx === nothing
        ald[k] = String[]
    end
    # Add connections to the adjacency list
    for j in 1:size(dam, 2)
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

ceqs, instream_csets = ModelingToolkit.generate_connection_equations_and_stream_connections(csets)
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
al = Dict{String, Vector{String}}()    # Adjacency list for effort variable
str2con = Dict{String, ModelingToolkit.ConnectionElement}()

vtype = get_connection_type(csets[1].set[1].v)

# Generate name to systems from csets
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
end

# filter out flow
filterflow = true
if filterflow
    filter!(d -> occursin("e(t)", d[1]), str2con)
end

# Generate adjacency matrix from csets
k2idx = Dict(k => i for (i, k) in enumerate(keys(str2con)))

am = zeros(Int, length(str2con), length(str2con))
for cset in csets
    h = 0
    for (i, e) in enumerate(cset.set)
        @unpack sys, v, isouter = e
        k = String(nameof(sys)) * "₊" * String(Symbol(v))
        if String(Symbol(v)) != "e(t)" && filterflow
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
idx2k = [i => k for (i, k) in enumerate(keys(str2con))]
nm = [idx2k[i][2] for i in 1:length(idx2k)]
# Plot connection graph
graphplot(am, names=nm, nodeshape=:rect, size=(700, 700), method=:stress)

al = Dict{String, Vector{String}}()    # Adjacency list for effort variable
am_con = sum(am, dims=1)
for (i, v) in idx2k
    i, v = idx2k[1]
    if am_con[i] > 0
        idx = get(al, v, nothing)
        if idx === nothing
            al[v] = String[]
        end
        push!(al[v], [k for (j, k) in idx2k[am[:, i] .> 0]])
    end
end



am

sum(am, dims=1)

idx2k[1]

csets = deepcopy(connectionsets)
# al = csets2adjlist(csets, "f(t)");
al = csets2adjlist(csets, "");
if filterflow
    filter!(x -> occursin("e(t)", x[1]), al)
end
csets = adjlist2csets(al, str2con)
al

csets










al

tol = 1e-10
debug = false

sys = mdl
ceqs, instream_csets = ModelingToolkit.generate_connection_equations_and_stream_connections(csets)
ceqs
_sys = ModelingToolkit.expand_instream(instream_csets, sys; debug = debug, tol = tol)
sys = flatten(sys, true)
ModelingToolkit.@set! sys.eqs = [equations(_sys); ceqs]

equations(sys)

equations(alias_elimination(sys))

generate_graph(mdl, :f)
