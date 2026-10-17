# Two masses and a coupling spring-damper, joined by a native 0-junction.
for (case, force, q) in (("mechanical_msd_2dof_unforced", 0.0, 1.0),
                         ("mechanical_msd_2dof_forced", sin(t), 0.0))
    @named m = Inertia(I = 1.0)
    @named s = Compliance(C = 1.0)
    @named d = Resistance(R = 1.0)
    @named f = Se(force)
    @named m1 = Junction1([-1, m], [-1, d], [-1, s], [-1, :out])
    @named m2 = Junction1([-1, m], f, :in)
    @named s2 = Compliance(C = 1.0, q = q)
    @named sd = Junction1([-1, s2], [-1, d], :in)
    @named model = System([connect(m1.out, m2.in, sd.in)], t; systems = [m1, m2, sd])
    check_case(case, model, [
        "m1_j1₊m₊power₊f(t)" => model.m1.m.power.f,
        "m1_j1₊s₊q(t)" => model.m1.s.q,
        "m2_j1₊m₊power₊f(t)" => model.m2.m.power.f,
        "m2_sd₊s₊q(t)" => model.sd.s2.q,
    ])
end
