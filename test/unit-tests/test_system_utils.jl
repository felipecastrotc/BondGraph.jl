using BondGraph
using Test
using Symbolics
using ModelingToolkit

import BondGraph: t, D
import BondGraph: substitute_dict
import ModelingToolkit: get_defaults

@testset "simplifysys" begin
    # Build a dummy system
    dflt = Dict(:mass₊I => exp(1.0), Symbol("mass₊power₊f(t)") => π/2, :damper₊R => π)
    @named mass = Mass(m=dflt[:mass₊I ], u=dflt[Symbol("mass₊power₊f(t)")])
    @named damper = Damper(c=dflt[:damper₊R])
    @named sys = Junction1(mass, damper)

    # Simplify system
    @named sys_simp = simplifysys(sys)

    # The system should have the following expression
    eq_truth = [D(mass.power.f) ~ (-damper.R*mass.power.f)/mass.I]
    @test equations(sys_simp) == eq_truth

    # Check parameter default values
    dflt_simp = ModelingToolkit.get_defaults(sys_simp)
    dflt_simp = Dict(Symbol(k) => v  for (k, v) in dflt_simp)
    @test dflt_simp == dflt
end

@testset "substitute_dict" begin

    @parameters a, b, c
    @variables x(t), y(t), z(t)

    # Simple expression test
    sub = Dict(a => c, b => z)
    expr = a + b + c
    @test Symbol(substitute_dict(expr, expr, sub)) == Symbol(2*c + z)
    # More complex expression test
    sub = Dict(a => c, b => z, x => y)
    expr = (a + b + c)/((x*y)^a)
    @test Symbol(substitute_dict(expr, expr, sub)) == Symbol((2*c + z)/(y^(2*c)))
end

# TODO: I need to find a system where the structural_simplify does not expand the observed