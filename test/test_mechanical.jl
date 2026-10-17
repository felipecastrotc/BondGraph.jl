@testset "Mechanical shorthands match role-named elements" begin
    function oscillator(make_elements)
        I, C, R = make_elements()
        @named se = Se(2.0)
        Junction1(se, [-1, I], [-1, C], [-1, R]; name = :model)
    end

    role = oscillator() do
        (Inertia(I = 2.0, p = 1.0, name = :I),
         Compliance(C = 0.25, q = 0.3, name = :C),
         Resistance(R = 0.5, name = :R))
    end
    shorthand = oscillator() do
        (Mass(m = 2.0, u = 0.5, name = :I),
         Spring(k = 4.0, x = 0.3, name = :C),
         Damper(c = 0.5, name = :R))
    end
    times = range(0.0, 1.0, length = 11)
    role_sol = solve(ODEProblem(mtkcompile(role), [], (0.0, 1.0)), Tsit5(); abstol = 1e-10, reltol = 1e-10)
    short_sol = solve(ODEProblem(mtkcompile(shorthand), [], (0.0, 1.0)), Tsit5(); abstol = 1e-10, reltol = 1e-10)
    @test short_sol(times; idxs = shorthand.I.power.f) ≈ role_sol(times; idxs = role.I.power.f) atol = 1e-8
    @test short_sol(times; idxs = shorthand.C.q) ≈ role_sol(times; idxs = role.C.q) atol = 1e-8

    default_role = oscillator() do
        (Inertia(I = 1.0, p = 0.0, name = :I),
         Compliance(C = 0.1, q = 0.0, name = :C),
         Resistance(R = 10, name = :R))
    end
    default_short = oscillator() do
        (Mass(name = :I), Spring(name = :C), Damper(name = :R))
    end
    default_role_sol = solve(ODEProblem(mtkcompile(default_role), [], (0.0, 1.0)), Tsit5())
    default_short_sol = solve(ODEProblem(mtkcompile(default_short), [], (0.0, 1.0)), Tsit5())
    @test default_short_sol(times; idxs = default_short.I.power.f) ≈
          default_role_sol(times; idxs = default_role.I.power.f) atol = 1e-8
    @test default_short_sol(times; idxs = default_short.C.q) ≈
          default_role_sol(times; idxs = default_role.C.q) atol = 1e-8
end

@testset "Inertia power changes stored energy" begin
    @named source = Se(3.0)
    @named inertia = Inertia(I = 2.0, p = 1.0)
    @named model = Junction1(source, [-1, inertia])
    times = range(0.0, 1.0, length = 101)
    sol = solve(ODEProblem(mtkcompile(model), [], (0.0, 1.0)), Tsit5(); abstol = 1e-10, reltol = 1e-10)
    p = sol(times; idxs = model.inertia.p)
    energy = p .^ 2 ./ 4   # E = p²/(2I), Karnopp §3.1
    power = sol(times; idxs = model.inertia.power.e) .* sol(times; idxs = model.inertia.power.f)
    dt = step(times)
    @test maximum(abs, (energy[3:end] .- energy[1:end-2]) ./ (2dt) .- power[2:end-1]) < 1e-8
end

@testset "Two inertias driven by an effort share one flow (DAE)" begin
    @named source = Se(3.0)
    @named first_inertia = Inertia(I = 1.0)
    @named second_inertia = Inertia(I = 2.0)
    @named model = Junction1(source, [-1, first_inertia], [-1, second_inertia])
    times = range(0.0, 1.0, length = 11)
    # Closed form: (I1 + I2) df/dt = E with f(0) = 0, so f = t, p1 = I1 f, p2 = I2 f.
    # One momentum follows the other (derivative causality), so both default initial
    # values over-determine initialization; they agree (p1 = p2 = 0).
    prob = ODEProblem(mtkcompile(model), [], (0.0, 1.0); warn_initialize_determined = false)
    sol = solve(prob, Rodas5P(); abstol = 1e-10, reltol = 1e-10)
    @test sol(times; idxs = model.first_inertia.p) ≈ times atol = 1e-8
    @test sol(times; idxs = model.second_inertia.p) ≈ 2 .* times atol = 1e-8
    @test sol(times; idxs = model.first_inertia.power.f) ≈ times atol = 1e-8
end

@testset "Two inertias share a 1-junction flow" begin
    @named source = Sf(t)
    @named first_inertia = Inertia(I = 1.0)
    @named second_inertia = Inertia(I = 2.0)
    @named model = Junction1(source, [-1, first_inertia], [-1, second_inertia])
    times = range(0.0, 1.0, length = 11)
    # Both momenta are fixed by the source (derivative causality), so their default
    # initial values over-determine initialization; they agree with F(0) = 0.
    prob = ODEProblem(mtkcompile(model), [], (0.0, 1.0); warn_initialize_determined = false)
    sol = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10)
    @test sol(times; idxs = model.first_inertia.p) ≈ times atol = 1e-8
    @test sol(times; idxs = model.second_inertia.p) ≈ 2 .* times atol = 1e-8
end
