# The bond graph of a model as data and text, drawn with bgplot, and saved with savebg.
using BondGraphToolkit.ModelingToolkit: t_nounits as t

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

@testset "bondgraph shows as text" begin
    text = repr("text/plain", bondgraph(tf_model()))
    @test occursin("Se:se ⇀ 1:j", text) && occursin("1:j ⇀ TF:tf", text) && occursin("TF:tf ⇀ I:I", text)
end

@testset "bgplot: SVG drawing, DOT and TikZ source" begin
    p = bgplot(tf_model())
    @test showable("image/svg+xml", p)
    svg = repr("image/svg+xml", p)
    @test occursin("<svg", svg) && occursin("TF:tf", svg)
    @test bgplot(bondgraph(tf_model())) isa typeof(p)   # a BondGraph works too

    dot = sprint(print, bgplot(tf_model(); format = :dot))
    @test occursin("rankdir=LR", dot)
    @test count("->", dot) == 4                    # one edge per bond
    @test count("lvee", dot) == 4                  # drawn as half-arrows
    @test occursin("label=\"TF:tf\"", dot)

    tp = bgplot(tf_model(); format = :tikz)
    @test !showable("image/svg+xml", tp)
    @test repr("text/plain", tp) == sprint(print, tp)   # shown as its source
    tex = sprint(print, tp)
    @test occursin("\\begin{tikzpicture}", tex)
    @test count("\\node", tex) == 5              # one per junction, element, two-port
    @test count("\\draw[-bgarrow]", tex) == 4    # one half-arrow per bond
    @test occursin("{TF: \$tf\$}", tex) && occursin("{1}", tex)
    # Textbook layout: Se left of the junction on the power path, R above it.
    at = Dict(m[3] => parse.(Float64, (m[1], m[2]))
              for m in eachmatch(r"at \(([-\d.]+), ([-\d.]+)\) \{(.+)\};", tex))
    j, se, R = at["1"], at["Se: \$se\$"], at["R: \$R\$"]
    @test se[2] ≈ j[2] && se[1] < j[1]
    @test R[1] ≈ j[1] && R[2] > j[2]

    @test_throws ArgumentError bgplot(tf_model(); format = :jpg)
end

@testset "savebg: format from the file extension" begin
    mktempdir() do dir
        magic = Dict(".svg" => "<?xml", ".png" => "\x89PNG", ".pdf" => "%PDF",
                     ".tex" => "\\documentclass", ".dot" => "digraph")
        for (ext, head) in magic
            path = savebg(joinpath(dir, "m" * ext), tf_model())
            @test startswith(read(path, String), head)
        end
        @test isfile(savebg(joinpath(dir, "g.svg"), bondgraph(tf_model())))
        @test isfile(savebg(joinpath(dir, "p.tex"), bgplot(tf_model())))
        @test_throws ArgumentError savebg(joinpath(dir, "m.jpg"), tf_model())
    end
end
