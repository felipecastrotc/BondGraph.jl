# Readable equations: the model's equations in its states, with display names, and LaTeX.
using BondGraphToolkit.ModelingToolkit: t_nounits as t, D_nounits as D

function rc_model()
    @named R = Resistance(R = 1.0)
    @named C = Compliance(C = 1.0, q = 1.0)
    @named model = Junction0([-1, R], [-1, C])
end
# Names of the variables and parameters in an equation (derivative terms left out).
names_in(eq) = filter(!startswith("Differential"), string.(Symbolics.get_variables(eq.lhs - eq.rhs)))

@testset "equations in terms of the states" begin
    eqs = bg_equations(rc_model())
    @test length(eqs) == 1   # q̇ = -q/(RC): junction and port variables are gone
    @test sort(names_in(eqs[1])) == ["C₊C", "C₊q(t)", "R₊R"]
end

@testset "rename, also inside derivatives" begin
    @variables x(t) ẋ(t)
    @parameters τ
    model = rc_model()
    eqs = bg_equations(model; rename = Dict(model.C.q => x, model.R.R => τ))
    @test isequal(eqs[1].lhs, D(x))
    @test sort(names_in(eqs[1])) == ["C₊C", "x(t)", "τ"]
    eqs = bg_equations(model; rename = Dict(D(model.C.q) => ẋ, model.C.q => x))
    @test isequal(eqs[1].lhs, ẋ)
    # keys may also come from the compiled model
    sys = mtkcompile(model)
    @test isequal(bg_equations(sys; rename = Dict(sys.C.q => x))[1].lhs, D(x))
    @test_throws ArgumentError bg_equations(model; rename = Dict(model.R.power.e => x))
end

@testset "LaTeX" begin
    @variables x(t) ẋ(t)
    model = rc_model()
    tex = latexify(bg_equations(model; rename = Dict(D(model.C.q) => ẋ, model.C.q => x)))
    @test occursin("\\dot{x}", tex)
    @test !occursin("C.q", tex)
end
