# Model from examples/mechanical_msd_2dof.jl (unforced and forced cases)
function build_mechanical_msd_2dof()
    @named m = Mass(m = 1.0)
    @named s = Spring(k = 1.0)
    @named d = Damper(c = 1.0)
    @parameters F(t)
    @named f = Se(F)
    @named m1_j1 = Junction1([-1, m], [-1, d], [-1, s])
    @named m2_j1 = Junction1([-1, m], f)
    @named m2_j0 = Junction0()
    @named m2_sd = Junction1([-1, s], [-1, d])
    cons = [
        connect(m2_j0.power, m2_j1.power),
        connect(m2_j0.power, m2_sd.power),
        connect(m1_j1.power, m2_j0.power),
    ]
    @named mdl = ODESystem(cons, t)
    mdl = compose(mdl, m1_j1, m2_j0, m2_j1, m2_sd)
    @named sys = simplifysys(mdl)
    unforced = substitute(sys, Dict(m2_j1.f.F => 0.0))
    forced = substitute(sys, Dict(m2_j1.f.F => sin(t)))
    return sys, [
        "mechanical_msd_2dof_unforced" =>
            ODEProblem(unforced, [m2_j1.m.power.f => 0.0, m2_sd.s.q => 1.0], (0, 10)),
        "mechanical_msd_2dof_forced" =>
            ODEProblem(forced, [m2_j1.m.power.f => 0.0, m2_sd.s.q => 0.0], (0, 10)),
    ]
end

check_example("mechanical_msd_2dof", build_mechanical_msd_2dof)
