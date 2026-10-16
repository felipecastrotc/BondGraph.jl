module BondGraphToolkitGraphvizExt

using BondGraphToolkit, Graphviz_jll

"""The SVG text of a rendered bond graph; it displays as an image."""
struct BondGraphSVG
    svg::String
end
Base.show(io::IO, ::MIME"image/svg+xml", g::BondGraphSVG) = print(io, g.svg)

graphviz(dot, format) = read(pipeline(`$(Graphviz_jll.dot()) -T$format`; stdin = IOBuffer(dot)), String)

BondGraphToolkit.bondgraph_svg(sys) = BondGraphSVG(graphviz(bondgraph_dot(sys), "svg"))

end
