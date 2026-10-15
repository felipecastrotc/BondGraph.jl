# Bond-graph model from examples/electrical_dc_motor.jl (the hand-written
# reference ODE in that example is not captured; it does not use BondGraph).
function build_electrical_dc_motor()
    @named L = Mass(m = 0.5)
    @named R = Damper(c = 1.0)
    @named V = Se(12.0)
    @named J = Mass(m = 0.01)
    @named c = Damper(c = 0.1)
    @named τ = Se(1.0)
    @named jm = Junction1(τ, [-1, c], [-1, J])
    @named je = Junction1(V, [-1, R], [-1, L])
    eqs = []
    @named gy = mGY(je, jm, g = 0.01, coneqs = eqs)
    @named mdl = ODESystem(eqs, t)
    mdl = compose(mdl, gy, je, jm)
    @named sys = simplifysys(mdl)
    prob = ODEProblem(sys, [], (0, 10))
    return sys, ["electrical_dc_motor" => prob]
end

check_example("electrical_dc_motor", build_electrical_dc_motor)
