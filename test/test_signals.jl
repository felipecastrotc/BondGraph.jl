# Signals (activated bonds): sources and moduli driven by MTK standard-library blocks,
# and detectors measuring a bond without carrying power.

# `model` plus a block whose output drives `input`, as one system.
with_block(model, block, input) =
    System([connect(block.output, input)], t; name = :top, systems = [model, block])

@testset "Se(): effort from a signal (MSe)" begin
    function series_rc(se)
        @named R = Resistance(R = 1.0)
        @named C = Compliance(C = 1.0)
        Junction1(se => :pout, R, C; name = :m)
    end
    @named se = Se(sin(2π * t))
    reference = series_rc(se)
    @named se = Se()
    @named sine = Blocks.Sine(frequency = 1.0)
    m = series_rc(se)
    top = with_block(m, sine, m.se.input)
    @test solve_tp(top)(TT; idxs = top.m.C.q) ≈ solve_tp(reference)(TT; idxs = reference.C.q) atol = 1e-8
end

@testset "Sf(): flow from a signal (MSf)" begin
    @named sf = Sf()
    @named C = Compliance(C = 1.0)
    @named m = Junction0(sf => :pout, C)
    @named k = Blocks.Constant(k = 2.0)
    top = with_block(m, k, m.sf.input)
    @test solve_tp(top)(TT; idxs = top.m.C.q) ≈ 2 .* TT atol = 1e-8   # q = F t
end

@testset "TF(): modulus from a signal (MTF)" begin
    @named tf = TF()
    @named I = Inertia(I = 1.0)
    m = driven(tf, I; E = 3.0)
    @named k = Blocks.Constant(k = 2.0)
    top = with_block(m, k, m.tf.input)
    @test solve_dae(top)(TT; idxs = top.m.I.p) ≈ 1.5 .* TT atol = 1e-8   # as TF(m = 2)
end

@testset "Df measures a flow, De an effort, without changing the model" begin
    @named se = Se(2.0)
    @named R = Resistance(R = 4.0)
    @named df = Df()
    @named m = Junction1(se => :pout, R, df)   # Df sits on a 1-junction
    sol = solve_tp(m)
    @test sol(TT; idxs = m.df.output.u) ≈ fill(0.5, length(TT)) atol = 1e-10   # f = E/R
    @named se = Se(2.0)
    @named de = De()
    @named R = Resistance(R = 4.0)
    @named m = Junction0(se => :pout, R, de)   # De sits on a 0-junction
    @test solve_tp(m)(TT; idxs = m.de.output.u) ≈ fill(2.0, length(TT)) atol = 1e-10
end

@testset "Closed loop: De → blocks → GY modulus" begin
    # Se(2) ⇀ 1 ⇀ GY ⇀ 0 ⇀ C with modulus r = 2 + 0.5 e_C read by a detector, so
    # q̇ = 2/(2 + 0.5 q) and q = -4 + 2√(4 + 2t).
    function loop(gy)
        @named C = Compliance(C = 1.0)
        @named de = De()
        @named out = Junction0(:in => :in, C, de)
        @named se = Se(2.0)
        @named src = Junction1(se => :pout, :out)
        System([connect(src.out, gy.p1), connect(gy.p2, out.in)], t;
               name = :m, systems = [src, gy, out])
    end
    @named gy = GY()
    m = loop(gy)
    @named gain = Blocks.Gain(k = 0.5)
    @named offset = Blocks.Constant(k = 2.0)
    @named add = Blocks.Add()
    @named top = System([connect(m.out.de.output, gain.input), connect(gain.output, add.input1),
                         connect(offset.output, add.input2), connect(add.output, m.gy.input)], t;
                        systems = [m, gain, offset, add])
    @test solve_dae(top)(TT; idxs = top.m.out.C.q) ≈ -4 .+ 2 .* sqrt.(4 .+ 2 .* TT) atol = 1e-8
end
