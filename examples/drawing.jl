# # [Drawing Bond Graphs](@id drawing)
#
# A model built with BondGraphToolkit knows its own bond graph: junctions record their bonds and signs, and elements and two-ports are recognised by their ports. `bondgraph` returns it as data, and three functions draw it:
#
# - `bondgraph_dot`: Graphviz DOT text, for any Graphviz tool;
# - `bondgraph_svg`: the drawing itself, shown in notebooks and in these docs;
# - `bondgraph_tikz`: a standalone LaTeX/TikZ document in the style of these docs' figures, for a report or paper.
#
# The last two need Graphviz, which comes with `Graphviz_jll`: `using Graphviz_jll` loads them.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t
using Graphviz_jll

# The [DC motor](@ref dc-motor):
@named L = Inertia(I = 0.5)
@named R = Resistance(R = 1.0)
@named V = Se(12.0)
@named J = Inertia(I = 0.01)
@named c = Resistance(R = 0.1)
@named τ = Se(1.0)
@named je = Junction1(V, [-1, R], [-1, L], [-1, :gy])
@named jm = Junction1(τ, [-1, c], [-1, J], :gy)
@named gy = GY(r = 0.01)
@named motor = System([connect(je.gy, gy.p1), connect(gy.p2, jm.gy)], t; systems = [je, gy, jm])

# ## The bond graph as data
#
# Nodes are labelled `type:name` as in Karnopp, Margolis & Rosenberg's figures, and each bond is a `from => to` pair following its half-arrow:
g = bondgraph(motor)
g.nodes
#-
g.bonds

# ## Drawing it
#
# Bonds are drawn as half-arrows, left to right, with the power path through junctions and two-ports kept straight:
bondgraph_svg(motor)

# A `connect` of more than two free ports is itself a 0-junction, and gets its own `0` node in the drawing.
#
# ## TikZ for a paper
#
# `bondgraph_tikz` places the nodes with Graphviz's layout and writes a document that compiles with `pdflatex`. Copy its `tikzpicture`, with the `bgarrow` declaration, into your own document, and adjust positions or labels as you like:
tex = bondgraph_tikz(motor)
print(tex)
