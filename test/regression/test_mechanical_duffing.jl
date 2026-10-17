# Linear and cubic compliances share the 1-junction.

@oneport function CubicCompliance(; name, C = 0.2, q = 0.0)
    @parameters C = C
    @variables q(t) = q
    [power.e ~ q^3 / C, D(q) ~ power.f]
end

for (case, force, displacement) in
    (("mechanical_duffing_unforced", 0.0, 1.0),
     ("mechanical_duffing_forced", 8.0cos(0.5t), 0.0))
    @named m = Inertia(I = 1.0)
    @named s = Compliance(C = 1.0, q = displacement)
    @named s3 = CubicCompliance(C = 0.2, q = displacement)
    @named d = Resistance(R = 0.02)
    @named f = Se(force)
    @named model = Junction1(f, [-1, m], [-1, d], [-1, s], [-1, s3])
    check_case(case, model, ["m₊power₊f(t)" => model.m.power.f,
                             "s₊q(t)" => model.s.q,
                             "s3₊q(t)" => model.s3.q])
end
