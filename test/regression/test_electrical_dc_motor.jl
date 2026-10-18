# Model from examples/electrical_dc_motor.jl: V ⇀ 1ₑ ⇀ GY ⇀ 1ₘ, with R, L on the
# electrical 1-junction and the load torque τ, damping c and rotor inertia J on the
# mechanical one.

function dc_motor(; modulated = false)
    @named L = Inertia(I = 0.5)
    @named R = Resistance(R = 1.0)
    @named V = Se(12.0)
    @named J = Inertia(I = 0.01)
    @named c = Resistance(R = 0.1)
    @named τ = Se(1.0)
    @named je = Junction1(V => :pout, R, L, :gy)
    @named jm = Junction1(τ => :pout, c, J, :gy => :in)
    if modulated   # MGY: the ratio reads the shaft speed
        @named gy = GY(r = 0.01 + 0.001 * jm.J.f)
    else
        @named gy = GY(r = 0.01)
    end
    @named motor = System([connect(je.gy, gy.p1), connect(gy.p2, jm.gy)], t; systems = [je, gy, jm])
end

let motor = dc_motor()
    check_case("electrical_dc_motor", motor,
               ["je₊L₊power₊f(t)" => motor.je.L.f, "jm₊J₊power₊f(t)" => motor.jm.J.f])
end

# No golden data for a modulated gyrator (the ratio reads the shaft speed): it must
# still conserve power, e₁f₁ = e₂f₂ with f₂ = -p2.f.
@testset "modulated GY conserves power" begin
    motor = dc_motor(modulated = true)
    sol = solve(ODEProblem(mtkcompile(motor), [], (0.0, 10.0)), Tsit5(); abstol = 1e-10, reltol = 1e-10)
    gy = motor.gy
    p_in = sol(TS; idxs = gy.p1.e) .* sol(TS; idxs = gy.p1.f)
    p_out = sol(TS; idxs = gy.p2.e) .* -sol(TS; idxs = gy.p2.f)
    @test maximum(abs, p_in - p_out) < 1e-8
    @test maximum(abs, p_in) > 0.1   # power does flow through the gyrator
end
