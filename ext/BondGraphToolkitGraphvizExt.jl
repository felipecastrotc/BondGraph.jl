module BondGraphToolkitGraphvizExt

using BondGraphToolkit, Graphviz_jll
using BondGraphToolkit: bondgraph, dot_text, latexify

"""The SVG text of a rendered bond graph; it displays as an image."""
struct BondGraphSVG
    svg::String
end
Base.show(io::IO, ::MIME"image/svg+xml", g::BondGraphSVG) = print(io, g.svg)

graphviz(dot, format) = read(pipeline(`$(Graphviz_jll.dot()) -T$format`; stdin = IOBuffer(dot)), String)

BondGraphToolkit.bondgraph_svg(sys) = BondGraphSVG(graphviz(bondgraph_dot(sys), "svg"))

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

function BondGraphToolkit.bondgraph_tikz(sys)
    g = bondgraph(sys)
    names = Dict(n.id => "n$k" for (k, n) in enumerate(g.nodes))
    # `dot -Tplain` lines: `node <name> <x> <y> ...`, positions in inches.
    plain = graphviz(dot_text(g, n -> names[n.id]), "plain")
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

end
