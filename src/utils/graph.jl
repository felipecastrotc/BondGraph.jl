using Plots
using GraphRecipes
import ModelingToolkit: ConnectionSet, ConnectionElement
import ModelingToolkit: namespaced_var, generate_connection_set!


# Main functions

"""
    generate_bg_eqs!(connectionsets)

Generates bond graph equations from connection sets.

# Arguments:
- connectionsets: Array of connection sets.

# Returns:
- eqs: Array of equations representing the bond graph.
"""
function generate_bg_eqs!(connectionsets)
    # Get the bond graph connection sets
    bgconnectionsets = get_bg_connection_set!(connectionsets)

    # Create a dictionary mapping connection set names to connection elements
    str2con = csets2dict(bgconnectionsets)

    # Generate the adjacency matrix from the bond graph connection sets
    am = csets2adjmtx(bgconnectionsets, str2con)

    # Generate the bond graph equations from the adjacency matrix and connection elements
    eqs = adjmtx2eqs(am, str2con)

    return eqs
end


"""
    generate_graph(mdl, var=:e; method=:stress)

Generates a graph visualization of the bond graph model.

# Arguments:
- mdl: Bond graph model.
- var: Variable to visualize connections for. Can be :e (effort) or :f (flow). Default is :e.
- method: Graph layout method. Default is :stress.

# Example:
```julia-repl
julia> generate_graph(mdl, var=:e, method=:stress)
```
"""
function generate_graph(mdl, var=:e; method=:stress)

    connectionsets = ConnectionSet[]

    # Generate the bond graph connection sets
    sys = generate_connection_set!(connectionsets, mdl, nothing, nothing)

    # Get the bond graph connection sets
    bgconnectionsets = get_bg_connection_set!(connectionsets)

    # Create a dictionary mapping connection set names to connection elements
    str2con = csets2dict(bgconnectionsets)

    # Filter the connection sets based on the variable
    varin = var == :e ? "e(t)" : "f(t)"
    am = csets2adjmtx(bgconnectionsets, str2con, filterstr=varin, filterflow=true)

    # Create a dictionary mapping index to connection set names
    idx2k = Dict(i => k for (i, k) in enumerate(keys(str2con)))

    # Generate node names for visualization
    # TODO: check -> Indexing with indices obtained from `length`, `size` etc is discouraged. Use `eachindex` or `axes` 
    nm = [idx2k[i] for i = 1:length(idx2k)]

    # Generate the graph visualization using GraphPlot.jl
    graphplot(
        am,
        names=nm,
        nodeshape=:rect,
        size=(800, 600),
        method=method,
        arrow=arrow(:simple, :head, 1, 1),
        curves=false,
    )
end


# Graph-based algorithm functions

"""
    csets2dict(csets)

Converts connection sets to a dictionary.

# Arguments:
- csets: Array of connection sets.

# Returns:
- str2con: Dictionary mapping connection set names to connection elements.
"""
function csets2dict(csets)
    str2con = Dict{String,ConnectionElement}()

    for cset in csets
        for ele in cset.set
            k = String(Symbol(namespaced_var(ele)))
            str2con[k] = ele
        end
    end

    return str2con
end


"""
    csets2adjmtx(csets, str2con; filterstr="f(t)", filterflow=false)

Generates an adjacency matrix from connection sets and a dictionary.

# Arguments:
- csets: Array of connection sets.
- str2con: Dictionary mapping connection set names to connection elements.
- filterstr: String to filter connections by name (default: "f(t)").
- filterflow: Boolean indicating whether to filter connections by flow (default: false).

# Returns:
- am: Adjacency matrix representing the connections between connection elements.

# Example:
```julia-repl
julia> csets = [ConnectionSet([:a], [:b]), ConnectionSet([:c], [:d])]
julia> str2con = csets2dict(csets)
julia> am = csets2adjmtx(csets, str2con)
```
"""
function csets2adjmtx(csets, str2con; filterstr="f(t)", filterflow=false)

    # Filter str2con
    if filterflow
        filter!(d -> occursin(filterstr, d[1]), str2con)
    end

    # Create a dictionary to map connection element keys to indices
    k2idx = Dict(k => i for (i, k) in enumerate(keys(str2con)))

    # Initialize the adjacency matrix with zeros
    am = zeros(Int, length(str2con), length(str2con))

    # Iterate over each connection set
    for cset in csets
        h = 0
        for (i, e) in enumerate(cset.set)
            @unpack sys, v, isouter = e
            k = String(nameof(sys)) * "₊" * String(Symbol(v))

            if String(Symbol(v)) != filterstr && filterflow
                continue
            end

            # Assign the index of the "from" connection element if it's the first element
            if i == 1
                h = k2idx[k]
            else
                # Set the corresponding entry in the adjacency matrix to 1
                am[h, k2idx[k]] = 1
            end
        end
    end

    return am
end

"""
    get_sm_mtx!(am, idx2k, str2con)

Generates the signal matrix by modifying the adjacency matrix.

# Arguments:
- am: Adjacency matrix representing the connections between connection elements.
- idx2k: Dictionary mapping indices to connection element keys.
- str2con: Dictionary mapping connection set names to connection elements.

# Returns:
- sm: Signal matrix representing the connections with modified entries.

# Example:
```julia-repl
julia> am = [0 1 0; 0 0 1; 0 0 0]
julia> idx2k = Dict(1 => "a", 2 => "b", 3 => "c")
julia> str2con = Dict("a" => con1, "b" => con2, "c" => con3)
julia> sm = get_sm_mtx!(am, idx2k, str2con)
```
"""
function get_sm_mtx!(am, idx2k, str2con)
    # Get the signal matrix and remove the one-port elements as receiving node
    sm = deepcopy(am)

    # Find the connections where the receiving node is an one-port element
    cord = findall(am .> 0)
    filter!(x -> get_bg_junction(get_var(x[2], idx2k, str2con))[1] === op, cord)

    # Convert the connections and modify the signal matrix (sm)
    for c in cord
        flipc = CartesianIndex(c[2], c[1])
        am[c], am[flipc] = 0, 1
        sm[c], sm[flipc] = 0, -1
    end

    return sm
end


"""
    adjmtx2eqs(am, str2con)

Generates a set of equations from an adjacency matrix and connection element dictionary.

# Arguments:
- am: Adjacency matrix representing the connections between connection elements.
- str2con: Dictionary mapping connection set names to connection elements.

# Returns:
- eqs: Array of equations generated from the adjacency matrix and connection elements.
"""

function adjmtx2eqs(am, str2con)
    # Create a dictionary mapping indices to connection element keys
    idx2k = Dict(i => k for (i, k) in enumerate(keys(str2con)))

    # Generate the signal matrix
    sm = get_sm_mtx!(am, idx2k, str2con)

    # Calculate the input and output vectors
    in_vec = sum(am, dims=1)[:]
    out_vec = sum(am, dims=2)[:]

    idx_con = 1:size(am, 1)
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
            elseif jtype === tpgy || jtype === tptf
                push!(vin, get_var(j, idx2k, str2con))
            elseif jtype === op
                # Apply the signal matrix
                sgn = sm[j, i]
                var = get_var(j, idx2k, str2con)
                # The negative direction in BG refers to the product e*f
                sgn = get_bg_junction(var)[2] === bgflow ? abs(sgn) : sgn
                push!(vin, sgn * var)
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
                        "The $jtype type only allows one connection in and one connection out",
                    ),
                )
            end

            if lvin == 1
                # The negative direction in BG refers to the product e*f
                sgn = get_bg_junction(vin[1])[2] === bgflow ? -1 : 1
                push!(eqs, equalityeqs(vcat(sgn * vin, v))...)
            elseif lvout == 1
                push!(eqs, equalityeqs(vcat(vout, v))...)
            end
        end
    end

    return eqs
end
