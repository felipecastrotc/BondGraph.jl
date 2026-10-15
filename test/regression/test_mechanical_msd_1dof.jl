# Model from examples/mechanical_msd_1dof.jl (unforced and forced cases)
function build_mechanical_msd_1dof()
    @named m = Mass(m = 1.0)
    @named s = Spring(k = 1.0)
    @named d = Damper(c = 1.0)
    @parameters F(t)
    @named f = Se(F)
    @named mdl = Junction1(m, d, s, [-1, f])
    @named sys = simplifysys(mdl)
    unforced = substitute(sys, Dict(f.F => 0.0))
    forced = substitute(sys, Dict(f.F => sin(t)))
    return sys, [
        "mechanical_msd_1dof_unforced" =>
            ODEProblem(unforced, [m.power.f => 0.0, s.q => 1.0], (0, 10)),
        "mechanical_msd_1dof_forced" =>
            ODEProblem(forced, [m.power.f => 0.0, s.q => 0.0], (0, 10)),
    ]
end

check_example("mechanical_msd_1dof", build_mechanical_msd_1dof)
