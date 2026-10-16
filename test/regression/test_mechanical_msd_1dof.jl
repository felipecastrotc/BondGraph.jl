# The source bond enters the 1-junction; its I, C and R bonds leave it.
using BondGraphToolkit.ModelingToolkit: t_nounits as t

for (case, force, displacement) in
    (("mechanical_msd_1dof_unforced", 0.0, 1.0),
     ("mechanical_msd_1dof_forced", sin(t), 0.0))
    @named m = Inertia(I = 1.0)
    @named s = Compliance(C = 1.0, q = displacement)
    @named d = Resistance(R = 1.0)
    @named f = Se(force)
    @named model = Junction1(f, [-1, m], [-1, d], [-1, s])
    check_case(case, model, ["m₊power₊f(t)" => model.m.power.f,
                             "s₊q(t)" => model.s.q])
end
