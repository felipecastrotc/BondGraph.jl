# The bond graph of a model as data, as Graphviz DOT text, and rendered to SVG.
using BondGraphToolkit.ModelingToolkit: t_nounits as t
using Graphviz_jll

# Se ⇀ 1 ⇀ R, and the junction's free port `out` ⇀ TF ⇀ I.
function tf_model()
    @named se = Se(1.0)
    @named R = Resistance()
    @named tf = TF(m = 2.0)
    @named I = Inertia()
    @named j = Junction1(se, [-1, R], [-1, :out])
    @named model = System([connect(j.out, tf.p1), connect(tf.p2, I.power)], t; systems = [j, tf, I])
end

@testset "bondgraph: elements, junctions, two-ports and half-arrows" begin
    g = bondgraph(tf_model())
    @test Dict(n.id => n.label for n in g.nodes) ==
          Dict("j" => "1", "j₊se" => "Se:se", "j₊R" => "R:R", "tf" => "TF:tf", "I" => "I:I")
    @test Set(g.bonds) == Set(["j₊se" => "j", "j" => "j₊R", "j" => "tf", "tf" => "I"])
end

@testset "bondgraph: a connect of several free ports is a 0-junction" begin
    @named a = Se(1.0)
    @named b = Resistance()
    @named c = Compliance()
    @named j1 = Junction1(a, [-1, :out])   # bond leaves j1 ...
    @named j2 = Junction1(:in, [-1, b])    # ... and enters j2 and j3
    @named j3 = Junction1(:in, [-1, c])
    @named model = System([connect(j1.out, j2.in, j3.in)], t; systems = [j1, j2, j3])
    g = bondgraph(model)
    zero = only(n.id for n in g.nodes if n.label == "0")
    @test Set(g.bonds) == Set(["j1₊a" => "j1", "j1" => zero, zero => "j2", zero => "j3",
                               "j2" => "j2₊b", "j3" => "j3₊c"])
end

@testset "bondgraph_dot and bondgraph_svg" begin
    dot = bondgraph_dot(tf_model())
    @test occursin("rankdir=LR", dot)
    @test count("->", dot) == 4                    # one edge per bond
    @test count("arrowhead=lvee", dot) == 4        # drawn as half-arrows
    @test occursin("label=\"TF:tf\"", dot)
    svg = sprint(show, MIME("image/svg+xml"), bondgraph_svg(tf_model()))
    @test occursin("<svg", svg) && occursin("TF:tf", svg)
end

@testset "bondgraph_tikz" begin
    tex = bondgraph_tikz(tf_model())
    @test occursin("\\begin{tikzpicture}", tex)
    @test count("\\node", tex) == 5              # one per junction, element, two-port
    @test count("\\draw[-bgarrow]", tex) == 4    # one half-arrow per bond
    @test occursin("{TF: \$tf\$}", tex) && occursin("{1}", tex)
end
