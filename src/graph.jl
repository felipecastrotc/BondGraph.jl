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

# A bond-graph node: its subsystem path, `:junction`/`:element`/`:twoport`, and label.
const Node = NamedTuple{(:id, :kind, :label), Tuple{String, Symbol, String}}

subpath(path, name) = isempty(path) ? string(name) : string(path, "₊", name)

"""
    BondGraph

The bond graph of a model, returned by [`bondgraph`](@ref): fields `nodes` and `bonds`.
Displayed as the list of its bonds; draw it with [`bgplot`](@ref).
"""
struct BondGraph
    nodes::Vector{Node}
    bonds::Vector{Pair{String, String}}
end

"""
    bondgraph(sys)

The bond graph of the model `sys` as data, a [`BondGraph`](@ref) with fields `nodes`
and `bonds`, displayed as the list of its bonds (`Se:V ⇀ 1:je`).

- `nodes`: one `(; id, kind, label)` per junction (`kind = :junction`, label `"0"` or
  `"1"`), element (`:element`) and two-port (`:twoport`), labelled `type:name` as in
  Karnopp's figures (`"C:s"`, `"TF:lever"`). `id` is the path to the subsystem.
- `bonds`: one `from => to` pair of node ids per bond, following its half-arrow.

A `connect` of more than two power ports, or of two that both expect the bond to point
into them (e.g. two elements: a plain `connect` to an element is a bond pointing into
it), is a 0-junction and gets its own node.

A subsystem with a `power` port is drawn as one element, even when it is a composite
model inside: that's how hand-written elements are recognised.
"""
function bondgraph(sys)
    nodes = Node[]
    bonds = Pair{String, String}[]
    ports = Dict{String, Tuple{String, Symbol}}()   # port path => (node id, :into | :outof)
    connections = Vector{String}[]
    function walk!(sys, path)
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
    BondGraph(nodes, bonds)
end

# A node as written in the text display: `type:name`, junctions `1:name`; the 0-junction
# of a plain `connect` has no name.
function text_label(n)
    n.kind != :junction && return n.label
    startswith(n.id, "0_") ? n.label : string(n.label, ":", last(split(n.id, "₊")))
end

function Base.show(io::IO, ::MIME"text/plain", g::BondGraph)
    name = Dict(n.id => text_label(n) for n in g.nodes)
    print(io, "BondGraph with ", length(g.nodes), " nodes and ", length(g.bonds), " bonds:")
    foreach(((a, b),) -> print(io, "\n  ", name[a], " ⇀ ", name[b]), g.bonds)
end

# DOT text of a bond graph `g`, naming each node `name(node)`.
function dot_text(g, name)
    io = IOBuffer()
    println(io, "digraph bondgraph {\n  rankdir=LR;\n  splines=false;\n  nodesep=0.6;\n  ranksep=0.6;\n  node [shape=plaintext];")
    foreach(n -> println(io, "  ", name(n), " [label=\"", n.label, "\"];"), g.nodes)
    # Textbook layout: junctions and two-ports form a straight left-to-right power path, and
    # each element sits beside its node: in a free horizontal slot, then above, then below.
    node = Dict(n.id => n for n in g.nodes)
    leaf(id) = node[id].kind == :element
    taken = Dict(n.id => Symbol[] for n in g.nodes)
    for (a, b) in g.bonds
        leaf(a) || leaf(b) || (push!(taken[a], :right); push!(taken[b], :left))
    end
    edge(u, v, bond, attrs = "") = println(io, "  ", name(node[u]), " -> ", name(node[v]),
        bond == (u => v) ? " [arrowhead=lvee" : " [dir=back, arrowtail=lvee", attrs, "];")
    for (a, b) in g.bonds
        if !(leaf(a) || leaf(b))
            edge(a, b, a => b, ", weight=10")
            continue
        end
        el, j = leaf(a) ? (a, b) : (b, a)
        k = findfirst(∉(taken[j]), [:left, :right, :above, :below])
        slot = k === nothing ? :none : [:left, :right, :above, :below][k]
        slot === :none || push!(taken[j], slot)
        # dot puts an edge's head after its tail: to the right, or below within a rank.
        u, v = slot in (:left, :above) ? (el, j) : slot === :none ? (a, b) : (j, el)
        edge(u, v, a => b)
        slot in (:above, :below) && println(io, "  {rank=same; ", name(node[u]), "; ", name(node[v]), "}")
    end
    print(io, "}")
    String(take!(io))
end

# Graphviz's rendering of DOT text in `format` (svg, png, pdf, plain), as bytes.
graphviz(dot, format) = read(pipeline(`$(Graphviz_jll.dot()) -T$format`; stdin = IOBuffer(dot)))

dot_source(g) = dot_text(g, n -> string('"', n.id, '"'))

# The half-arrow of the docs' figures (docs/src/misc/tikz/base_tikzset.tex).
const TIKZ_PREAMBLE = raw"""
\documentclass[tikz]{standalone}
\usetikzlibrary{arrows.meta}
\pgfdeclarearrow{
  name = bgarrow,
  parameters = { \the\pgfarrowlength },
  setup code = {
    \pgfarrowssettipend{0\pgfarrowlength}
    \pgfarrowssetlineend{0\pgfarrowlength}
    \pgfarrowssetvisualbackend{.5\pgfarrowlength}
    \pgfarrowssavethe\pgfarrowlength
  },
  drawing code = {
    \pgfpathmoveto{\pgfpoint{0\pgfarrowlength}{0pt}}
    \pgfpathlineto{\pgfpoint{-10pt}{-0.45\pgfarrowlength}}
    \pgfusepathqstroke
  },
  defaults = { length = 5mm }
}
"""

# `type:name` → `type: $name$`, the name in math mode as in the docs' figures (Greek too).
function tikz_label(label)
    occursin(':', label) || return label
    type, name = split(label, ':'; limit = 2)
    string(type, ": ", latexify(Symbol(name)))
end

# A standalone TikZ document, placed with Graphviz's layout.
function tikz_source(g)
    names = Dict(n.id => "n$k" for (k, n) in enumerate(g.nodes))
    # `dot -Tplain` lines: `node <name> <x> <y> ...`, positions in inches.
    plain = String(graphviz(dot_text(g, n -> names[n.id]), "plain"))
    pos = Dict(f[2] => (f[3], f[4]) for f in split.(split(plain, '\n')) if !isempty(f) && f[1] == "node")
    io = IOBuffer()
    print(io, TIKZ_PREAMBLE, "\\begin{document}\n\\begin{tikzpicture}[x=1in, y=1in, line width=0.25mm]\n")
    for n in g.nodes
        x, y = pos[names[n.id]]
        println(io, "  \\node[inner sep=4pt] (", names[n.id], ") at (", x, ", ", y, ") {", tikz_label(n.label), "};")
    end
    foreach(((a, b),) -> println(io, "  \\draw[-bgarrow] (", names[a], ") -- (", names[b], ");"), g.bonds)
    print(io, "\\end{tikzpicture}\n\\end{document}\n")
    String(take!(io))
end

# A bond graph in `format`, as a string or bytes.
function source(g, format)
    format === :dot ? dot_source(g) : format === :tikz ? tikz_source(g) : graphviz(dot_source(g), format)
end

"""
    BondGraphPlot

What [`bgplot`](@ref) returns: a bond graph and the format to show it in. Save it with
[`savebg`](@ref); `print` writes its source.
"""
struct BondGraphPlot
    graph::BondGraph
    format::Symbol
end

graph(g::BondGraph) = g
graph(p::BondGraphPlot) = p.graph
graph(sys) = bondgraph(sys)

"""
    bgplot(x; format = :svg)

Draw the bond graph of `x`, a model or a [`BondGraph`](@ref), the way `plot` draws a
solution. The drawing follows the textbook layout: the power path through junctions and
two-ports runs left to right, and each element sits beside its node, with half-arrow
bonds, `0`/`1` junctions and `type:name` labels.

- `format = :svg`: the drawing. It appears in the VS Code plot pane, Jupyter, Pluto and the
  docs; in a terminal REPL it opens in the default viewer.
- `format = :tikz`: a standalone LaTeX/TikZ document in the style of the docs' figures,
  shown as text. Compile it with `pdflatex`, or copy its `tikzpicture` (and the `bgarrow`
  declaration) into a paper.
- `format = :dot`: the Graphviz DOT source, shown as text.

Save any of them with [`savebg`](@ref).
"""
function bgplot(x; format = :svg)
    format in (:svg, :tikz, :dot) ||
        throw(ArgumentError("bgplot: format must be :svg, :tikz or :dot, not :$format"))
    BondGraphPlot(graph(x), format)
end

Base.print(io::IO, p::BondGraphPlot) = write(io, source(p.graph, p.format))
Base.showable(::MIME"image/svg+xml", p::BondGraphPlot) = p.format === :svg
Base.show(io::IO, ::MIME"image/svg+xml", p::BondGraphPlot) = print(io, p)
function Base.show(io::IO, ::MIME"text/plain", p::BondGraphPlot)
    if p.format === :svg
        print(io, "BondGraphPlot (SVG): savebg(\"file.svg\", p) to save it")
    else
        print(io, p)
    end
end

# A terminal can't show SVG: open it in the default viewer, as Plots opens a GR window.
function Base.display(d::REPL.REPLDisplay, p::BondGraphPlot)
    p.format === :svg || return display(d, MIME("text/plain"), p)
    path = savebg(tempname() * ".svg", p)
    cmd = Sys.isapple() ? `open $path` : Sys.iswindows() ? `cmd /c start "" $path` : `xdg-open $path`
    run(cmd; wait = false)
    nothing
end

const EXTENSIONS = Dict(".svg" => :svg, ".png" => :png, ".pdf" => :pdf, ".tex" => :tikz, ".dot" => :dot)

"""
    savebg(path, x)

Save the bond graph of `x` (a model, a [`BondGraph`](@ref) or a [`bgplot`](@ref)) to
`path`, in the format of its extension: `.svg`, `.png`, `.pdf` (rendered by Graphviz),
`.tex` (TikZ, see [`bgplot`](@ref)) or `.dot`. Returns `path`.
"""
function savebg(path, x)
    ext = lowercase(splitext(path)[2])
    haskey(EXTENSIONS, ext) ||
        throw(ArgumentError("savebg: unknown extension \"$ext\"; use .svg, .png, .pdf, .tex or .dot"))
    write(path, source(graph(x), EXTENSIONS[ext]))
    path
end
