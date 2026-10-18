# # [Drawing Bond Graphs](@id drawing)
#
# A model built with BondGraphToolkit knows its own bond graph: junctions record their bonds and signs, and elements and two-ports are recognised by their ports. Three functions work with it, the way `solve`, `plot` and `savefig` work with a simulation:
#
# - `bondgraph`: the bond graph as data, shown as the list of its bonds;
# - `bgplot`: the drawing, or its TikZ or Graphviz DOT source;
# - `savebg`: saves it to a file, in the format of the file's extension.

using BondGraphToolkit
using BondGraphToolkit.ModelingToolkit: t_nounits as t

# The [DC motor](@ref dc-motor):
@named L = Inertia(I = 0.5)
@named R = Resistance(R = 1.0)
@named V = Se(12.0)
@named J = Inertia(I = 0.01)
@named c = Resistance(R = 0.1)
@named τ = Se(1.0)
@named je = Junction1(V => :pout, R, L, :gy)
@named jm = Junction1(τ => :pout, c, J, :gy => :in)
@named gy = GY(r = 0.01)
@named motor = System([connect(je.gy, gy.p1), connect(gy.p2, jm.gy)], t; systems = [je, gy, jm])

# ## The bond graph as data
#
# `bondgraph` lists the bonds, each along its half-arrow, with nodes labelled `type:name` as in Karnopp, Margolis & Rosenberg's figures:
g = bondgraph(motor)

# Its `nodes` and `bonds` are plain data, for your own checks or drawings:
g.bonds

# ## Drawing it
#
# `bgplot` draws the bond graph in the textbook layout: the power path through junctions and two-ports runs left to right, and each element sits beside its junction. The drawing appears in the VS Code plot pane, in notebooks and in these docs; in a terminal REPL it opens in your default viewer.
bgplot(motor)

# A `connect` of more than two free ports is itself a 0-junction, and gets its own `0` node in the drawing.
#
# ## TikZ for a paper
#
# `format = :tikz` gives a standalone LaTeX/TikZ document in the style of these docs' figures, placed with the same layout; it compiles with `pdflatex`. Copy its `tikzpicture`, with the `bgarrow` declaration, into your own document, and adjust positions or labels as you like. `format = :dot` gives the Graphviz DOT source.
bgplot(motor; format = :tikz)

# ## Saving
#
# `savebg` takes the format from the file extension: `.svg`, `.png`, `.pdf`, `.tex` (TikZ) or `.dot`.
#
# ```julia
# savebg("motor.pdf", motor)
# savebg("motor.tex", motor)
# ```
