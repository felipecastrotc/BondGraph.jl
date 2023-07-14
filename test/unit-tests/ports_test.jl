using BondGraph
using BondGraph: t, D, bg
using ModelingToolkit
using ModelingToolkit: getdefault
using DifferentialEquations
using Test
using Symbolics

name = :test
split_str = "₊"

@testset "Se" begin
    # Default values
    cte = exp(1.0)
    @parameters ω, A
    expr = A*cos(2π*ω*t)
    # Function
    @named test_cte = Se(cte)
    @named test_expr = Se(expr)
    @test test.name == name
    test = Mass(name=name, m = m, u = u)
    @test test.name == name
    # Manual
    @named power = Power(flow=u)
    @parameters I = m
    # Test the power defaults
    @test getdefault(test.power.e) == 0.0
    @test getdefault(test.power.f) ≈ u
    @test getdefault(test.I) ≈ m
    # Test equation generated
    @test [D(power.f) ~ power.e/I] == equations(test)
    # Check the parameters
    ps = parameters(test)
    @test length(ps) == 1
    @test Symbol(ps[1]) == Symbol(I)
    # Check if the naming is corret
    naming = split(String(Symbolics.tosymbol(test.I)), "₊")
    @test Symbol(naming[1]) == name
    @test Symbol(naming[2]) == Symbolics.tosymbol(I)
end

@testset "Sf" begin
    # Default values
    u = π
    m = exp(1.0)
    # Function
    @named test = Mass(m = m, u = u)
    @test test.name == name
    test = Mass(name=name, m = m, u = u)
    @test test.name == name
    # Manual
    @named power = Power(flow=u)
    @parameters I = m
    # Test the power defaults
    @test getdefault(test.power.e) == 0.0
    @test getdefault(test.power.f) ≈ u
    @test getdefault(test.I) ≈ m
    # Test equation generated
    @test [D(power.f) ~ power.e/I] == equations(test)
    # Check the parameters
    ps = parameters(test)
    @test length(ps) == 1
    @test Symbol(ps[1]) == Symbol(I)
    # Check if the naming is corret
    naming = split(String(Symbolics.tosymbol(test.I)), "₊")
    @test Symbol(naming[1]) == name
    @test Symbol(naming[2]) == Symbolics.tosymbol(I)
end

@testset "Inertance" begin
    # Default values
    u = π
    m = exp(1.0)
    # Function
    @named test = Mass(m = m, u = u)
    @test test.name == name
    test = Mass(name=name, m = m, u = u)
    @test test.name == name
    # Manual
    @named power = Power(flow=u)
    @parameters I = m
    # Test the power defaults
    @test getdefault(test.power.e) == 0.0
    @test getdefault(test.power.f) ≈ u
    @test getdefault(test.I) ≈ m
    # Test equation generated
    @test [D(power.f) ~ power.e/I] == equations(test)
    # Check the parameters
    ps = parameters(test)
    @test length(ps) == 1
    @test Symbol(ps[1]) == Symbol(I)
    # Check if the naming is corret
    naming = split(String(Symbolics.tosymbol(test.I)), "₊")
    @test Symbol(naming[1]) == name
    @test Symbol(naming[2]) == Symbolics.tosymbol(I)
end


@testset "Compliance" begin
    # Default values
    x = π
    k = exp(1.0)
    # Function
    @named test = Spring(k=k, x=x)
    @test test.name == name
    test = Spring(name=name, k=k, x=x)
    @test test.name == name
    # Manual
    @named power = Power()
    @variables q(t) = x
    @parameters C = 1/k
    # Test the power defaults
    @test getdefault(test.power.e) == 0.0
    @test getdefault(test.power.f) == 0.0
    @test getdefault(test.C) ≈ 1/k
    # Test equation generated
    @test [power.e ~ q / C, D(q) ~ power.f] == equations(test)
    # Check the parameters
    ps = parameters(test)
    @test length(ps) == 1
    @test Symbol(ps[1]) == Symbol(C)
    # Check if the naming is corret
    naming = split(String(Symbolics.tosymbol(test.C)), "₊")
    @test Symbol(naming[1]) == name
    @test Symbol(naming[2]) == Symbolics.tosymbol(C)
end


@testset "Resistance" begin
    # Default values
    u = π
    c = exp(1.0)
    # Function
    @named test = Damper(c = c, u = u)
    @test test.name == name
    test = Damper(name=name, c = c, u = u)
    @test test.name == name
    # Manual
    @named power = Power(flow=u)
    @parameters R = c
    # Test the power defaults
    @test getdefault(test.power.e) == 0.0
    @test getdefault(test.power.f) ≈ u
    @test getdefault(test.R) ≈ c
    # Test equation generated
    @test [power.e ~ power.f * R]  == equations(test)
    # Check the parameters
    ps = parameters(test)
    @test length(ps) == 1
    @test Symbol(ps[1]) == Symbol(R)
    # Check if the naming is corret
    naming = split(String(Symbolics.tosymbol(test.R)), "₊")
    @test Symbol(naming[1]) == name
    @test Symbol(naming[2]) == Symbolics.tosymbol(R)
end

