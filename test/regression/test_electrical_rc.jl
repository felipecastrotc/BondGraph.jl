# Model from examples/electrical_rc.jl
function build_electrical_rc()
    @named C = Spring(k = 1 / 1.0)
    @named R = Damper(c = 1.0)
    @named model = Junction0(R, C)
    @named sys = simplifysys(model)
    prob = ODEProblem(sys, [R.power.f => 0.0, C.q => 1.0], (0, 10))
    return sys, ["electrical_rc" => prob]
end

check_example("electrical_rc", build_electrical_rc)
