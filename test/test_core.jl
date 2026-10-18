# Core elements and junctions through the public API, checked against closed-form
# solutions and the junction laws (Karnopp Eqs. 3.16-3.21). All parameters are 1.

# e·f on an element's bond. (MTK 11 can't index a solution by a product of observed
# variables, so the product is taken after sampling.)
power(sol, x) = sol(T1; idxs = x.e) .* sol(T1; idxs = x.f)

@testset "Se drives R and C in series on a 1-junction" begin
    @named se = Se(1.0)
    @named R = Resistance(R = 1.0)
    @named C = Compliance(C = 1.0)
    @named m = Junction1(se => :pout, R, C)
    sol = simulate(m)
    @test sol(T1; idxs = m.C.q) ≈ 1 .- exp.(-T1) atol = 1e-8   # q = CE(1 - e^{-t/RC})
    # power in from the source = power out to R and C (bᵢ = +1, -1, -1)
    @test maximum(abs, power(sol, m.se) - power(sol, m.R) - power(sol, m.C)) < 1e-8
end

@testset "Sf fills C on a 0-junction" begin
    @named sf = Sf(2.0)
    @named C = Compliance(C = 1.0)
    @named m = Junction0(sf => :pout, C)
    @test simulate(m)(T1; idxs = m.C.q) ≈ 2 .* T1 atol = 1e-8    # q = F t
end

@testset "Inertia stores momentum p" begin
    @named se = Se(3.0)
    @named I = Inertia(I = 2.0)
    @named m = Junction1(se => :pout, I)
    sol = simulate(m)
    @test sol(T1; idxs = m.I.p) ≈ 3 .* T1 atol = 1e-8            # ṗ = e
    @test sol(T1; idxs = m.I.f) ≈ 1.5 .* T1 atol = 1e-8    # f = p / I
end

@testset "RC on a 0-junction: power balance and a reversed resistor" begin
    @named R = Resistance(R = 1.0)
    @named C = Compliance(C = 1.0, q = 1.0)
    @named m = Junction0(R, C)
    sol = simulate(m)
    @test sol(T1; idxs = m.C.q) ≈ exp.(-T1) atol = 1e-8
    @test maximum(abs, -power(sol, m.R) - power(sol, m.C)) < 1e-8
    # R's bond drawn into the junction: its law stays e = R f in the bond's
    # variables, so it delivers power, a negative resistor (ticket 06, Q2)
    @named R = Resistance(R = 1.0)
    @named C = Compliance(C = 1.0, q = 1.0)
    @named m = Junction0(R => :pout, C)
    sol = simulate(m)
    @test sol(T1; idxs = m.C.q) ≈ exp.(T1) atol = 1e-8
    @test all(>(0), power(sol, m.R)[2:end])   # e·f > 0: power flows along the bond, out of R
end

@testset "Junction-to-junction bond through free ports" begin
    @named se = Se(1.0)
    @named R = Resistance(R = 1.0)
    @named C = Compliance(C = 1.0)
    @named j1 = Junction1(se => :pout, R, :out)   # bond leaves j1 ...
    @named j0 = Junction0(:in => :in, C)              # ... and enters j0
    @named m = System([connect(j1.out, j0.in)], t; systems = [j1, j0])
    @test simulate(m)(T1; idxs = m.j0.C.q) ≈ 1 .- exp.(-T1) atol = 1e-8
end

@testset "a custom @oneport matches the built-in Compliance" begin
    @oneport function MyC(; name, C = 1.0, q = 0.0)
        @parameters C = C
        @variables q(t) = q
        [e ~ q / C, D(q) ~ f]
    end
    rc(c) = (@named R = Resistance(R = 1.0); Junction0(R, c; name = :m))
    custom = rc(MyC(C = 2.0, q = 1.0; name = :C))
    builtin = rc(Compliance(C = 2.0, q = 1.0; name = :C))
    @test simulate(custom)(T1; idxs = custom.C.q) ≈ simulate(builtin)(T1; idxs = builtin.C.q) atol = 1e-10
end

# MTK gives every connector's flow at least a singleton connection set, so an
# unconnected PowerPort carries no flow (ModelingToolkitBase, src/systems/connectors.jl,
# generate_connection_set!). Oriented ports (ticket 12) rely on this for the unused port.
@testset "An unconnected PowerPort carries no flow" begin
    @named C = Compliance(C = 1.0, q = 1.0)
    @named m = System(Equation[], t; systems = [C])
    sol = simulate(m)
    @test sol(T1; idxs = m.C.q) ≈ ones(length(T1)) atol = 1e-12
    @test sol(T1; idxs = m.C.f) ≈ zeros(length(T1)) atol = 1e-12
end

@testset "A plain connect of element ports is a signed 0-junction" begin
    @named sf = Sf(1.0)
    @named R = Resistance(R = 2.0)
    @named C = Compliance(C = 1.0)
    @named plain = System([connect(sf.pout, C.power, R.power)], t; systems = [sf, C, R])
    @test simulate(plain)(T1; idxs = plain.C.q) ≈ 2 .* (1 .- exp.(-T1 ./ 2)) atol = 1e-8   # q̇ = F - q/RC
end

@testset "The removed [-1, x] member form throws" begin
    @named R = Resistance()
    @test_throws ArgumentError Junction1([-1, R]; name = :j)
    @test_throws ArgumentError Junction1(:a => :pout; name = :j)
    @test_throws ArgumentError Junction1(R => :in; name = :j)
end
