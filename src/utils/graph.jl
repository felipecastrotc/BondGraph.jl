using Plots
using SparseArrays
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

    # Lists to store the rows, columns, and values for the sparse matrix
    rows = Int[]
    cols = Int[]
    vals = Int[]

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
                # Add the row, column, and value data for the sparse matrix
                push!(rows, h)
                push!(cols, k2idx[k])
                push!(vals, 1)
            end
        end
    end
    # Construct and return the sparse matrix
    return sparse(rows, cols, vals, length(str2con), length(str2con))
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

function adjmtx2eqs(adjacency_matrix, str2con)
    # Create a dictionary mapping indices to connection element keys
    idx2k = Dict(i => k for (i, k) in enumerate(keys(str2con)))

    # Initialize list to store generated equations
    generated_equations = []

    # Loop through the adjacency matrix
    for i in 1:size(adjacency_matrix, 1)
        var, node_type, node_property = retrieve_node_details(i, idx2k, str2con)

        # Check applicability of rules
        apply_equalityeqs = (node_type == j0 && node_property == bgeffort) ||
                            (node_type == j1 && node_property == bgflow)

        apply_sumvar = (node_type == j0 && node_property == bgflow) ||
                       (node_type == j1 && node_property == bgeffort)
        check_tp_con = (node_type === tpgy) || (node_type === tptf)

        # Skip iteration if no rules apply
        if !apply_equalityeqs && !apply_sumvar && !check_tp_con
            continue
        end

        # Identify connected nodes and respective signs
        connected_nodes = []
        connection_signs = []
        for j in 1:size(adjacency_matrix, 2)
            if adjacency_matrix[i, j] == 1
                push!(connected_nodes, j)
                push!(connection_signs, -1)  # Negative sign for leaving node
            elseif adjacency_matrix[j, i] == 1
                push!(connected_nodes, j)
                push!(connection_signs, 1)   # Positive sign for arriving node
            end
        end

        # Gather variables associated with connected nodes
        vars = []
        for (n, s) in zip(connected_nodes, connection_signs)
            aux_var, aux_type, _ = retrieve_node_details(n, idx2k, str2con)
            # Generate the connection name
            if aux_type == j0 || aux_type == j1
                aux_var = gen_con_name(var, aux_var, s)
            end
            # The sign assumption is for the product e(t)*f(t). Therefore, 
            # one of the power variables should not follow the sign assumption.
            # In this code I arbitrary defined the equality equations to not 
            # follow. It will be a mix and match of e(t) and f(t). I 
            # particularly find the sum of the variable easier to check for 
            # inconsistencies in the algorithm.
            s = apply_equalityeqs ? 1 : s
            push!(vars, s * aux_var)
        end

        # Apply rules and generate equations
        if apply_equalityeqs
            push!(generated_equations, equalityeqs(vars)...)
            # push!(generated_equations, equalityeqs(vars)...)
        elseif apply_sumvar
            # push!(generated_equations, sumvar(vcat(var, vars...)))
            push!(generated_equations, sumvar(vars))
        elseif check_tp_con
            if length(vars) != 1
                throw(
                    DomainError(
                        "The $node_type type only allows one connection in and one connection out",
                    ),
                )
            end
        end
    end

    return generated_equations
end


"""
`retrieve_node_details(i, idx2k, str2con)` fetches details about a node 
based on its index in the adjacency matrix.

Inputs:
    - i: Index of the node in the adjacency matrix.
    - idx2k: Dictionary mapping indices to connection element keys.
    - str2con: Dictionary mapping strings to connection information.
Outputs:
    - var: Extracted variable based on the given index.
    - node_type: Type of the node (e.g., j0, j1).
    - node_property: Property of the node (e.g., bgeffort, bgflow).
"""
function retrieve_node_details(i, idx2k, str2con)
    var = get_var(i, idx2k, str2con)
    node_type, node_property = get_bg_junction(var)
    return var, node_type, node_property
end


"""
`gen_con_name(base_var, con_var, con_sign)` generates a connection name based 
on the given inputs. It creates a new name that represents a directed connection
between the two input variable names.

Inputs:
    - base_var: The base variable name.
    - con_var: The connecting variable name.
    - con_sign: Integer indicating the direction of the connection.
                Positive indicates con_var -> base_var,
                Negative indicates base_var -> con_var.
Outputs:
    - Symbol representing the new connection name.
"""
function gen_con_name(base_var, con_var, con_sign)

    function extract_var_name(var_name)
        """
        Extracts the main part of the variable name and bg type (effort or flow).
        """
        var_str = String(getname(var_name))

        if occursin("power", var_str)
            name = chop(split(var_str, "power")[1], tail=1)
        else
            name = chop(var_str, tail=2)
        end
        return name, last(var_str)
    end

    # Check the validity of con_sign
    if con_sign == 0
        throw(ArgumentError("con_sign should be either positive or negative, but got 0."))
    end

    # Notes:
    # If con_sign > 0  the bg is: con_var -> base_var
    # If con_sign < 0  the bg is: base_var -> con_var
    # Determine connection direction based on con_sign
    left_var, right_var = con_sign > 0 ? (con_var, base_var) : (base_var, con_var)

    # Extract variable names
    left_name = extract_var_name(left_var)[1]
    right_name, last_name = extract_var_name(right_var)

    # Construct and return the new name
    return rename(left_var, Symbol(left_name * "▶" * right_name * "₊" * last_name))
end
