# Two-ports through the public API (Karnopp Eqs. 3.8-3.11): port 1 is the bond into the
# two-port, port 2 the bond out of it. A source reaches port 1 through a 1-junction
# (Se ⇀ 1 ⇀ TF), and port 2 bonds straight into a storage element (TF ⇀ I).

# Power into port 1 and out of port 2, in the bonds' own variables (f₂ = -p2.f).
power_in(sol, tp) = sol(TT; idxs = tp.p1.e) .* sol(TT; idxs = tp.p1.f)
power_out(sol, tp) = sol(TT; idxs = tp.p2.e) .* -sol(TT; idxs = tp.p2.f)

@testset "TF: e₁ = m e₂, f₂ = m f₁" begin
    @named tf = TF(m = 2.0)
    @named I = Inertia(I = 1.0)
    m = driven(tf, I; E = 3.0)
    sol = solve_tp(m)
    @test sol(TT; idxs = m.I.p) ≈ 1.5 .* TT atol = 1e-8        # ṗ = e₂ = E/m
    @test sol(TT; idxs = m.tf.p1.f) ≈ 0.75 .* TT atol = 1e-8   # f₁ = f₂/m = (p/I)/m
    @test maximum(abs, power_in(sol, m.tf) - power_out(sol, m.tf)) < 1e-8
end

@testset "GY: e₁ = r f₂, e₂ = r f₁" begin
    @named gy = GY(r = 4.0)
    @named C = Compliance(C = 1.0)
    m = driven(gy, C; E = 2.0)
    sol = solve_tp(m)
    @test sol(TT; idxs = m.C.q) ≈ 0.5 .* TT atol = 1e-8        # q̇ = f₂ = E/r
    @test sol(TT; idxs = m.gy.p1.f) ≈ 0.125 .* TT atol = 1e-8  # f₁ = e₂/r = (q/C)/r
    @test maximum(abs, power_in(sol, m.gy) - power_out(sol, m.gy)) < 1e-8
end

@testset "Modulated GY conserves power" begin
    @named gy = GY(r = 2 + t)
    @named C = Compliance(C = 1.0)
    m = driven(gy, C; E = 2.0)
    sol = solve_tp(m)
    @test sol(TT; idxs = m.C.q) ≈ 2 .* log.(1 .+ TT ./ 2) atol = 1e-8   # q̇ = E/(2+t)
    @test maximum(abs, power_in(sol, m.gy) - power_out(sol, m.gy)) < 1e-8
end

@testset "Modulated TF: the modulus is an expression" begin
    @named tf = TF(m = 1 + t)
    @named I = Inertia(I = 1.0)
    m = driven(tf, I; E = 3.0)
    sol = solve_tp(m)
    @test sol(TT; idxs = m.I.p) ≈ 3 .* log.(1 .+ TT) atol = 1e-8   # ṗ = E/(1+t)
    @test maximum(abs, power_in(sol, m.tf) - power_out(sol, m.tf)) < 1e-8
end

@testset "A custom @twoport matches the built-in TF" begin
    @twoport function MyTF(; name, m)
        @parameters m = m
        [e1 ~ m * e2, f2 ~ m * f1]
    end
    @named tf = MyTF(m = 2.0)
    @named I = Inertia(I = 1.0)
    custom = driven(tf, I; E = 3.0)
    @named tf = TF(m = 2.0)
    @named I = Inertia(I = 1.0)
    builtin = driven(tf, I; E = 3.0)
    @test solve_tp(custom)(TT; idxs = custom.I.p) ≈ solve_tp(builtin)(TT; idxs = builtin.I.p) atol = 1e-10
end
