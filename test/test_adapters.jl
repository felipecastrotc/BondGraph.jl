# Adapters to MTK standard-library networks, each checked against a reference model.
# An adapter's bond is drawn into it (`[-1, adapter]` at a junction): the MTKSL network
# is the element at the end of that bond.
import ModelingToolkitStandardLibrary.Electrical as El
import ModelingToolkitStandardLibrary.Mechanical.Rotational as Rot
import ModelingToolkitStandardLibrary.Mechanical.Translational as Tr
import ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition as TP
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as Hy

# Reference: Se(E) ⇀ 1 ⇀ R, C (series), all bond graph.
function series_rc(; E, R, C)
    @named se = Se(E)
    @named R = Resistance(R = R)
    @named C = Compliance(C = C)
    Junction1(se, [-1, R], [-1, C]; name = :ref)
end

@testset "ElectricalAdapter: MTKSL load" begin
    # Se ⇀ 1 ⇀ adapter, with an MTKSL resistor and capacitor in series behind it.
    @named se = Se(3.0)
    @named ea = ElectricalAdapter()
    @named j = Junction1(se, [-1, ea])
    @named res = El.Resistor(R = 2.0)
    @named cap = El.Capacitor(C = 0.5, v = 0.0)
    @named gnd = El.Ground()
    @named m = System([connect(j.ea.p, res.p), connect(res.n, cap.p), connect(cap.n, j.ea.n, gnd.g)], t;
                      systems = [j, res, cap, gnd])
    ref = series_rc(E = 3.0, R = 2.0, C = 0.5)
    s, r = solve_dae(m), solve_dae(ref)
    @test s(TT; idxs = m.cap.v) ≈ r(TT; idxs = ref.C.power.e) atol = 1e-8
    @test s(TT; idxs = m.res.i) ≈ r(TT; idxs = ref.R.power.f) atol = 1e-8
end

@testset "ElectricalAdapter: MTKSL source gives negative bond power" begin
    # The MTKSL voltage source supplies the bond graph, so with the bond drawn into the
    # adapter the bond's flow and the charge are the reference's with the opposite sign.
    @named k = Blocks.Constant(k = 3.0)
    @named src = El.Voltage()
    @named gnd = El.Ground()
    @named ea = ElectricalAdapter()
    @named R = Resistance(R = 2.0)
    @named C = Compliance(C = 0.5)
    @named j = Junction1([-1, ea], [-1, R], [-1, C])
    @named m = System([connect(k.output, src.V), connect(src.p, j.ea.p), connect(src.n, j.ea.n, gnd.g)], t;
                      systems = [k, src, gnd, j])
    ref = series_rc(E = 3.0, R = 2.0, C = 0.5)
    s, r = solve_dae(m), solve_dae(ref)
    @test s(TT; idxs = m.j.C.q) ≈ -r(TT; idxs = ref.C.q) atol = 1e-8
    @test all(<=(1e-12), s(TT; idxs = m.j.ea.power.e) .* s(TT; idxs = m.j.ea.power.f))
end

@testset "RotationalAdapter: DC motor driving an MTKSL load" begin
    # The DC motor of the examples, with the rotor inertia and friction moved to MTKSL.
    function motor(load)
        @named L = Inertia(I = 0.5)
        @named R = Resistance(R = 1.0)
        @named V = Se(12.0)
        @named τ = Se(1.0)
        @named je = Junction1(V, [-1, R], [-1, L], [-1, :gy])
        @named jm = Junction1(τ, :gy, load...)
        @named gy = GY(r = 0.01)
        System([connect(je.gy, gy.p1), connect(gy.p2, jm.gy)], t; name = :motor, systems = [je, gy, jm])
    end
    @named J = Inertia(I = 0.01)
    @named c = Resistance(R = 0.1)
    ref = motor([[-1, c], [-1, J]])
    @named ra = RotationalAdapter(phi = 0.0)
    inner = motor([[-1, ra]])
    @named inertia = Rot.Inertia(J = 0.01, w = 0.0)
    @named damper = Rot.Damper(d = 0.1)
    @named fixed = Rot.Fixed()
    @named m = System([connect(inner.jm.ra.flange, inertia.flange_a), connect(inertia.flange_b, damper.flange_b),
                       connect(damper.flange_a, fixed.flange)], t; systems = [inner, inertia, damper, fixed])
    s, r = solve_dae(m), solve_dae(ref)
    @test s(TT; idxs = m.inertia.w) ≈ r(TT; idxs = ref.jm.J.power.f) atol = 1e-8
    @test s(TT; idxs = m.motor.je.L.power.f) ≈ r(TT; idxs = ref.je.L.power.f) atol = 1e-8
end

@testset "ElectricalAdapter: MTKSL armature circuit driving a bond-graph motor" begin
    # Voltage, resistor and inductor in MTKSL feed GY ⇀ 1ₘ. GY's port 1 and the adapter
    # both expect the bond to point into them: a 0-junction where the MTKSL side supplies.
    @named k = Blocks.Constant(k = 12.0)
    @named src = El.Voltage()
    @named gnd = El.Ground()
    @named res = El.Resistor(R = 1.0)
    @named ind = El.Inductor(L = 0.5, i = 0.0)
    @named ea = ElectricalAdapter()
    @named gy = GY(r = 0.01)
    @named J = Inertia(I = 0.01)
    @named c = Resistance(R = 0.1)
    @named τ = Se(1.0)
    @named jm = Junction1(τ, [-1, c], [-1, J], :gy)
    @named m = System([connect(k.output, src.V), connect(src.p, res.p), connect(res.n, ind.p),
                       connect(ind.n, ea.p), connect(ea.n, src.n, gnd.g),
                       connect(ea.power, gy.p1), connect(gy.p2, jm.gy)], t;
                      systems = [k, src, gnd, res, ind, ea, gy, jm])
    # Reference: the all-bond-graph DC motor of the examples.
    @named L = Inertia(I = 0.5)
    @named R = Resistance(R = 1.0)
    @named V = Se(12.0)
    @named je = Junction1(V, [-1, R], [-1, L], [-1, :gy])
    @named J = Inertia(I = 0.01)
    @named c = Resistance(R = 0.1)
    @named τ = Se(1.0)
    @named jm = Junction1(τ, [-1, c], [-1, J], :gy)
    @named gy = GY(r = 0.01)
    @named ref = System([connect(je.gy, gy.p1), connect(gy.p2, jm.gy)], t; systems = [je, gy, jm])
    s, r = solve_dae(m), solve_dae(ref)
    @test s(TT; idxs = m.jm.J.power.f) ≈ r(TT; idxs = ref.jm.J.power.f) atol = 1e-8
    @test s(TT; idxs = m.ind.i) ≈ r(TT; idxs = ref.je.L.power.f) atol = 1e-8
end

@testset "Translational adapters: a force accelerates an MTKSL mass" begin
    # Se(1) ⇀ 1 ⇀ adapter ⇀ MTKSL mass of 2 kg: v = t/2, as with an MTKSL force source.
    @named k = Blocks.Constant(k = 1.0)
    @named force = Tr.Force()
    @named mass = Tr.Mass(m = 2.0)
    @named ref = System([connect(k.output, force.f), connect(force.flange, mass.flange)], t;
                        systems = [k, force, mass])
    @named se = Se(1.0)
    @named ta = TranslationalAdapter()
    @named j = Junction1(se, [-1, ta])
    @named m = System([connect(j.ta.flange, mass.flange)], t; systems = [j, mass])
    # Tr.Mass has no initial values; in `ref` mass.v is observed, so start from the force's.
    trsolve(x, v) = (sys = mtkcompile(x);
                     solve(ODEProblem(sys, [v(sys) => 0.0, sys.mass.s => 0.0], (0.0, 1.0)), Rodas5P();
                           abstol = 1e-10, reltol = 1e-10))
    s, r = trsolve(m, x -> x.mass.v), trsolve(ref, x -> x.force.flange.v)
    @test s(TT; idxs = m.mass.v) ≈ TT ./ 2 atol = 1e-8
    @test s(TT; idxs = m.mass.v) ≈ r(TT; idxs = ref.mass.v) atol = 1e-8
    @named force = TP.Force()
    @named mass = TP.Mass(m = 2.0, s = 0.0, v = 0.0)
    @named ref = System([connect(k.output, force.f), connect(force.flange, mass.flange)], t;
                        systems = [k, force, mass])
    @named se = Se(1.0)
    @named tpa = TranslationalPositionAdapter()
    @named j = Junction1(se, [-1, tpa])
    @named m = System([connect(j.tpa.flange, mass.flange)], t; systems = [j, mass])
    s, r = solve_dae(m), solve_dae(ref)
    @test s(TT; idxs = m.mass.v) ≈ TT ./ 2 atol = 1e-8
    @test s(TT; idxs = m.mass.s) ≈ r(TT; idxs = ref.mass.s) atol = 1e-8
end

@testset "HydraulicAdapter: MTKSL pressure source" begin
    # As the electrical source case: the MTKSL side supplies, so the volume in C is the
    # reference's with the opposite sign. Volumes are ~1e-4 m³, hence the tight abstol.
    p0, Rh, Ch = 2.0e5, 1.0e9, 2.0e-9
    @named k = Blocks.Constant(k = p0)
    @named src = Hy.Pressure()
    @named fluid = Hy.HydraulicFluid()
    @named ha = HydraulicAdapter()
    @named R = Resistance(R = Rh)
    @named C = Compliance(C = Ch)
    @named j = Junction1([-1, ha], [-1, R], [-1, C])
    @named m = System([connect(k.output, src.p), connect(fluid, src.port), connect(src.port, j.ha.port)], t;
                      systems = [k, src, fluid, j])
    ref = series_rc(E = p0, R = Rh, C = Ch)
    T = range(0.0, 10.0, length = 11)
    hsolve(x) = solve(ODEProblem(mtkcompile(x), [], (0.0, 10.0)), Rodas5P(); abstol = 1e-18, reltol = 1e-10)
    @test hsolve(m)(T; idxs = m.j.C.q) ≈ -hsolve(ref)(T; idxs = ref.C.q) rtol = 1e-6
end
