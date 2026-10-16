# A junction-to-junction bond: textbook signs and power balance at both ends.
@testset "Named free ports preserve junction signs" begin
    @named se = Se(1.0)
    @named R = Resistance(R = 1.0)
    @named C = Compliance(C = 1.0)
    @named j1 = Junction1(se, [-1, R], [-1, :out])
    @named j2 = Junction1(:in, [-1, C])
    @named model = System([connect(j1.out, j2.in)], t; systems = [j1, j2])
    sol = simulate(model)
    effort(x) = sol(T1; idxs = x.e)
    flow(x) = sol(T1; idxs = x.f)
    @test sol(T1; idxs = model.j2.C.q) ≈ 1 .- exp.(-T1) atol = 1e-8
    @test maximum(abs, effort(model.j1.se.power) .* flow(model.j1.se.power) .-
                       effort(model.j1.R.power) .* flow(model.j1.R.power) .+
                       effort(model.j1.out) .* flow(model.j1.out)) < 1e-8
    @test maximum(abs, effort(model.j2.in) .* flow(model.j2.in) .-
                       effort(model.j2.C.power) .* flow(model.j2.C.power)) < 1e-8
    @test maximum(abs, flow(model.j1.out) .+ flow(model.j2.in)) < 1e-8
end
