# Model from examples/mechanical_duffing.jl (unforced and forced cases)
@oneport function Spring3(; name, k = 1.0, x = 0.0)
    @variables q(t) = x
    @parameters C = 1 / k
    [
        power.e ~ q^3 / C
        D(q) ~ power.f
    ]
end

function build_mechanical_duffing()
    α, β, δ, γ, ω = 1.0, 5.0, 0.02, 8.0, 0.5
    @named m = Mass(m = 1.0)
    @named s = Spring(k = α)
    @named s3 = Spring3(k = β)
    @named d = Damper(c = δ)
    @parameters F(t)
    @named f = Se(F)
    @named mdl = Junction1([-1, m], [-1, d], [-1, s], [-1, s3], f)
    @named sys = simplifysys(mdl)
    unforced = substitute(sys, Dict(f.F => 0.0))
    forced = substitute(sys, Dict(f.F => γ * cos(ω * t)))
    return sys, [
        "mechanical_duffing_unforced" =>
            ODEProblem(unforced, [m.power.f => 0.0, s.q => 1.0, s3.q => 1.0], (0, 10)),
        "mechanical_duffing_forced" =>
            ODEProblem(forced, [m.power.f => 0.0, s.q => 0.0], (0, 10)),
    ]
end

check_example("mechanical_duffing", build_mechanical_duffing)
