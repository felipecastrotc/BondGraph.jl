# The bond graph of a model, read back from its systems: junctions (which record their
# members), elements (a `power` port), two-ports (`p1`, `p2`) and the `connect`
# statements between free ports, two-ports and elements.

const LETTERS = Dict("Resistance" => "R", "Compliance" => "C", "Inertia" => "I")

# An element's or two-port's bond-graph type: the constructor that built it.
function element_type(sys)
    meta = ModelingToolkit.get_gui_metadata(sys)
    type = meta === nothing ? "?" : last(split(string(meta.type), '.'))
    get(LETTERS, type, type)
end

subpath(path, name) = isempty(path) ? string(name) : string(path, "₊", name)

"""
    bondgraph(sys)

The bond graph of the model `sys` as data: `(; nodes, bonds)`.

- `nodes`: one `(; id, kind, label)` per junction (`kind = :junction`, label `"0"` or
  `"1"`), element (`:element`) and two-port (`:twoport`), labelled `type:name` as in
  Karnopp's figures (`"C:s"`, `"TF:lever"`). `id` is the path to the subsystem.
- `bonds`: one `from => to` pair of node ids per bond, following its half-arrow.

A `connect` of more than two power ports, or of two that both expect the bond to point
into them (e.g. two elements: a plain `connect` to an element is a bond pointing into
it), is a 0-junction and gets its own node.
"""
function bondgraph(sys)
    nodes = NamedTuple{(:id, :kind, :label), Tuple{String, Symbol, String}}[]
    bonds = Pair{String, String}[]
    ports = Dict{String, Tuple{String, Symbol}}()   # port path => (node id, :into | :outof)
    connections = Vector{String}[]
    walk!(sys, path) = begin
        id = isempty(path) ? string(nameof(sys)) : path
        subs = ModelingToolkit.get_systems(sys)
        names = nameof.(subs)
        meta = ModelingToolkit.getmetadata(sys, JunctionMembers, nothing)
        if meta !== nothing
            kind, members = meta
            push!(nodes, (; id, kind = :junction, label = string(kind)))
            for (sign, name, free) in members
                p = subpath(path, name)
                if free
                    ports[p] = (id, sign == 1 ? :into : :outof)
                else
                    walk!(subs[findfirst(==(name), names)], p)
                    push!(bonds, sign == 1 ? p => id : id => p)
                end
            end
        elseif :power in names
            push!(nodes, (; id, kind = :element, label = "$(element_type(sys)):$(nameof(sys))"))
            ports[subpath(path, :power)] = (id, :into)
        elseif :p1 in names && :p2 in names
            push!(nodes, (; id, kind = :twoport, label = "$(element_type(sys)):$(nameof(sys))"))
            ports[subpath(path, :p1)] = (id, :into)
            ports[subpath(path, :p2)] = (id, :outof)
        else
            foreach(s -> walk!(s, subpath(path, nameof(s))), subs)
            for eq in ModelingToolkit.get_eqs(sys)
                c = Symbolics.value(eq.rhs)
                c isa ModelingToolkit.Connection &&
                    push!(connections, [subpath(path, nameof(s)) for s in c.systems])
            end
        end
    end
    walk!(sys, "")
    for ends in connections
        ends = [ports[e] for e in ends if haskey(ports, e)]   # power ports only
        length(ends) < 2 && continue
        roles = last.(ends)
        if length(ends) == 2 && Set(roles) == Set([:into, :outof])
            push!(bonds, first(ends[roles .== :outof][1]) => first(ends[roles .== :into][1]))
        else
            j0 = "0_$(length(nodes) + 1)"
            push!(nodes, (; id = j0, kind = :junction, label = "0"))
            append!(bonds, [role == :into ? j0 => node : node => j0 for (node, role) in ends])
        end
    end
    (; nodes, bonds)
end

"""
    bondgraph_dot(sys)

The bond graph of `sys` (see [`bondgraph`](@ref)) as Graphviz DOT text: laid out left
to right, junctions drawn as `0`/`1`, elements and two-ports as `type:name`, bonds as
half-arrows. Render it with any Graphviz tool, or with [`bondgraph_svg`](@ref).
"""
bondgraph_dot(sys) = dot_text(bondgraph(sys), n -> string('"', n.id, '"'))

# DOT text of a bond graph `g`, naming each node `name(node)`.
function dot_text(g, name)
    io = IOBuffer()
    println(io, "digraph bondgraph {\n  rankdir=LR;\n  splines=false;\n  node [shape=plaintext];")
    foreach(n -> println(io, "  ", name(n), " [label=\"", n.label, "\"];"), g.nodes)
    # Bonds between junctions and two-ports are the main power path: keep them straight.
    node = Dict(n.id => n for n in g.nodes)
    for (a, b) in g.bonds
        weight = node[a].kind != :element && node[b].kind != :element ? ", weight=10" : ""
        println(io, "  ", name(node[a]), " -> ", name(node[b]), " [arrowhead=lvee", weight, "];")
    end
    print(io, "}")
    String(take!(io))
end

"""
    bondgraph_svg(sys)

The bond graph of `sys` rendered by Graphviz, displayed as SVG in notebooks and in the
docs. Needs `using Graphviz_jll`, which loads the package extension.
"""
function bondgraph_svg end

"""
    bondgraph_tikz(sys)

The bond graph of `sys` as a standalone LaTeX/TikZ document, placed with Graphviz's
layout and drawn in the style of this package's documentation figures (half-arrow
bonds, `0`/`1` junctions, `type:name` labels). Compile it with `pdflatex`, or copy its
`tikzpicture` (and the `bgarrow` declaration) into a paper. Needs `using Graphviz_jll`.
"""
function bondgraph_tikz end
