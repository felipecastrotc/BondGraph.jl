import ModelingToolkit: ConnectionSet, rename
import ModelingToolkit: namespaced_var, get_connection_type, getname


# Function to handle bond graph connections

"""
    check_bg_con(connectionset)

    Check if the connection set contains a bond-graph (bg) connection.

    Arguments:
    - connectionset: ConnectionSet object.

    Returns:
    - Boolean indicating whether the connection set contains a bond-graph connection.
"""
function check_bg_con(connectionset)
    ele = namespaced_var(connectionset.set[1])
    return get_connection_type(ele) === bg
end


"""
    get_bg_connection_set!(connectionsets)

    Extracts the bond-graph (bg) connection sets from a list of connection sets.

    Arguments:
    - connectionsets: Array of ConnectionSet objects.

    Returns:
    - Array of ConnectionSet objects containing only the bond-graph connection sets.
"""
function get_bg_connection_set!(connectionsets)
    bgconnectionsets = filter(check_bg_con, connectionsets)  # Filter connection sets with bond-graph connections
    filter!(x -> !check_bg_con(x), connectionsets)  # Remove the bond-graph connection sets from the original list
    return bgconnectionsets
end

# General connection util functions

function get_var(idx, idx2k, str2con)
    return namespaced_var(str2con[idx2k[idx]])
end

"""
    add_idx(var, idx)

    Adds an index to a variable by renaming it with the index appended to its name.

    Arguments:
    - var: Variable to be indexed.
    - idx: Index to be added to the variable name.

    Returns:
    - Indexed variable with the index appended to its name.

    Example:
    x = Variable(:x)
    indexed_var = add_idx(x, 1)  # Renames "x" to "x1"
"""
function add_idx(var, idx)
    v = deepcopy(var)
    vstr = String(getname(v))
    newname = Symbol(vstr * string(idx))
    return rename(v, newname)  # Rename the variable with the indexed name
end


"""
    gen_tp_con!(eqs, sys, subsys)

    Generates bond graph connections between a system and its sub-systems for the junction 0 and 1 (two-port elements).

    Arguments:
    - eqs: Array of equations to which the connections will be added.
    - sys: The main system to which the connections will be made.
    - subsys: Sub-systems to be connected to the main system.

    Returns:
    - None

    Example:
    eqs = Equation[]
    sys = ODESystem([:x, :y], [dx, dy])
    subsys = [ODESystem([:a], [da]), ODESystem([:b], [db])]
    gen_tp_con!(eqs, sys, subsys)
"""
function gen_tp_con!(eqs, sys, subsys)
    # Get connections
    c = subsys isa ODESystem ? [subsys, nothing] : collect(subsys)

    # Remove nothing from c array
    pos = .!isnothing.(c)
    c = c[pos]

    # Apply connections
    if length(c) > 0
        if sum(pos) == 1
            if pos[1]
                push!(eqs, connect(c[1].power, sys.pin))  # Connect sub-system's power port to the main system's input port
            else
                push!(eqs, connect(sys.pout, c[1].power))  # Connect main system's output port to the sub-system's power port
            end
        elseif length(c) == 2
            push!(eqs, connect(c[1].power, sys.pin), connect(sys.pout, c[2].power))  # Connect both sub-systems' power ports to the main system's input and output ports
        end
    end
end
