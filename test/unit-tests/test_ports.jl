using BondGraph
using BondGraph: t, D, bg, j0, j1
import BondGraph: get_bg_junction
using ModelingToolkit
using ModelingToolkit: getdefault
using DifferentialEquations
using Test
using Symbolics

name = :test
split_str = "₊"

@testset "Se-cte" begin
    # Default values
    expr = exp(1.0)
    # Naming
    @named test = Se(expr)
    # Check naming
    @test test.name == name
    test = Se(expr, name=name)
    @test test.name == name
    # Manual
    @named power = Power()
    # Test the power defaults
    @test getdefault(test.power.e) == 0.0
    @test getdefault(test.power.f) == 0.0
    # Test equation generated
    @test [power.e ~ expr] == equations(test)
    # Check the parameters
    ps = parameters(test)
    @test length(ps) == 0
end

@testset "Se-expr" begin
    # Default values
    @parameters ω, A
    expr = A * cos(2π * ω * t)
    # Naming
    @named test = Se(expr)
    # Check naming
    @test test.name == name
    test = Se(expr, name=name)
    @test test.name == name
    # Manual
    @named power = Power()
    # Test the power defaults
    @test getdefault(test.power.e) == 0.0
    @test getdefault(test.power.f) == 0.0
    # Test equation generated
    @test [power.e ~ expr] == equations(test)
    # Check the parameters
    ps = parameters(test)
    @test length(ps) == 2
    # Check if the naming is corret
    naming = split(String(Symbolics.tosymbol(test.ω)), "₊")
    @test Symbol(naming[1]) == name
    @test Symbol(naming[2]) == Symbolics.tosymbol(ω)

    naming = split(String(Symbolics.tosymbol(test.A)), "₊")
    @test Symbol(naming[1]) == name
    @test Symbol(naming[2]) == Symbolics.tosymbol(A)
end

@testset "Sf-cte" begin
    # Default values
    expr = exp(1.0)
    # Naming
    @named test = Sf(expr)
    # Check naming
    @test test.name == name
    test = Sf(expr, name=name)
    @test test.name == name
    # Manual
    @named power = Power()
    # Test the power defaults
    @test getdefault(test.power.e) == 0.0
    @test getdefault(test.power.f) == 0.0
    # Test equation generated
    @test [power.f ~ expr] == equations(test)
    # Check the parameters
    ps = parameters(test)
    @test length(ps) == 0
end

@testset "Sf-expr" begin
    # Default values
    @parameters ω, A
    expr = A * cos(2π * ω * t)
    # Naming
    @named test = Sf(expr)
    # Check naming
    @test test.name == name
    test = Sf(expr, name=name)
    @test test.name == name
    # Manual
    @named power = Power()
    # Test the power defaults
    @test getdefault(test.power.e) == 0.0
    @test getdefault(test.power.f) == 0.0
    # Test equation generated
    @test [power.f ~ expr] == equations(test)
    # Check the parameters
    ps = parameters(test)
    @test length(ps) == 2
    # Check if the naming is corret
    naming = split(String(Symbolics.tosymbol(test.ω)), "₊")
    @test Symbol(naming[1]) == name
    @test Symbol(naming[2]) == Symbolics.tosymbol(ω)

    naming = split(String(Symbolics.tosymbol(test.A)), "₊")
    @test Symbol(naming[1]) == name
    @test Symbol(naming[2]) == Symbolics.tosymbol(A)
end

@testset "Inertance" begin
    # Default values
    u = π
    m = exp(1.0)
    # Function
    @named test = Mass(m=m, u=u)
    @test test.name == name
    test = Mass(name=name, m=m, u=u)
    @test test.name == name
    # Manual
    @named power = Power(flow=u)
    @parameters I = m
    # Test the power defaults
    @test getdefault(test.power.e) == 0.0
    @test getdefault(test.power.f) ≈ u
    @test getdefault(test.I) ≈ m
    # Test equation generated
    @test [D(power.f) ~ power.e / I] == equations(test)
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
    @parameters C = 1 / k
    # Test the power defaults
    @test getdefault(test.power.e) == 0.0
    @test getdefault(test.power.f) == 0.0
    @test getdefault(test.C) ≈ 1 / k
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
    @named test = Damper(c=c, u=u)
    @test test.name == name
    test = Damper(name=name, c=c, u=u)
    @test test.name == name
    # Manual
    @named power = Power(flow=u)
    @parameters R = c
    # Test the power defaults
    @test getdefault(test.power.e) == 0.0
    @test getdefault(test.power.f) ≈ u
    @test getdefault(test.R) ≈ c
    # Test equation generated
    @test [power.e ~ power.f * R] == equations(test)
    # Check the parameters
    ps = parameters(test)
    @test length(ps) == 1
    @test Symbol(ps[1]) == Symbol(R)
    # Check if the naming is corret
    naming = split(String(Symbolics.tosymbol(test.R)), "₊")
    @test Symbol(naming[1]) == name
    @test Symbol(naming[2]) == Symbolics.tosymbol(R)
end

for sgn in [1, -1] 
    @testset "Junction0 $sgn" begin
        # Default values
        u = π
        m = exp(1.0)
        # Function
        @named mass = Mass(m=m, u=u)
        @named damper = Damper(c=m, u=u)
        # Naming
        @named test = Junction0(mass, [sgn, damper])
        # Check naming
        @test test.name == name
        test = Junction0(mass, [sgn, damper], name=name)
        @test test.name == name
        # Check type junction
        @test get_bg_junction(test.power.e)[1] === j0
        @test get_bg_junction(test.power.f)[1] === j0
        # Manual
        @named power = Power(type=j0)
        # Test the power defaults
        @test getdefault(test.power.e) == 0.0
        @test getdefault(test.power.f) == 0.0
        # Test equation generated
        eqs_func = equations(test)
        if sgn == -1
            eqs_manual = [
                connect(mass.power, power),
                connect(power, damper.power),
            ]
        else
            eqs_manual = [
                connect(mass.power, power),
                connect(damper.power, power),
            ]
        end
        eqs_elm = []
        push!(eqs_elm, equations(compose(ODESystem(Equation[], t, [], []; name=:mass), mass))...)
        push!(eqs_elm, equations(compose(ODESystem(Equation[], t, [], []; name=:damper), damper))...)
        push!(eqs_manual, eqs_elm...)

        eqs_manual
        @test length(eqs_func) == length(eqs_manual)
        @test all([Symbol(eqs_manual[i]) == Symbol(eqs_func[i]) for i in 1:length(eqs_func)])

        # Check the parameters
        ps = parameters(test)
        @test length(ps) == (length(parameters(mass)) +  length(parameters(damper)))

        # Expand equations
        eqs_expd = equations(expand_connections(test))
        # Check if the elements equations were maintaned
        @test all([eqs_expd[i] == eqs_elm[i] for i in 1:2])
        eqs_expd
        @test (0 ~ mass.power.f + damper.power.f) in eqs_expd
        @test (mass.power.e ~ sgn*damper.power.e) in eqs_expd
    end
end

for sgn in [1, -1]
    @testset "Junction1 $sgn" begin
        # Default values
        u = π
        m = exp(1.0)
        # Function
        @named mass = Mass(m=m, u=u)
        @named damper = Damper(c=m, u=u)
        # Naming
        @named test = Junction1(mass, [sgn, damper])
        # Check naming
        @test test.name == name
        test = Junction1(mass, [sgn, damper], name=name)
        @test test.name == name
        # Check type junction
        @test get_bg_junction(test.power.e)[1] === j1
        @test get_bg_junction(test.power.f)[1] === j1
        # Manual
        @named power = Power(type=j1)
        # Test the power defaults
        @test getdefault(test.power.e) == 0.0
        @test getdefault(test.power.f) == 0.0
        # Test equation generated
        eqs_func = equations(test)
        if sgn == -1
            eqs_manual = [
                connect(mass.power, power),
                connect(power, damper.power),
            ]
        else
            eqs_manual = [
                connect(mass.power, power),
                connect(damper.power, power),
            ]
        end
        eqs_elm = []
        push!(eqs_elm, equations(compose(ODESystem(Equation[], t, [], []; name=:mass), mass))...)
        push!(eqs_elm, equations(compose(ODESystem(Equation[], t, [], []; name=:damper), damper))...)
        push!(eqs_manual, eqs_elm...)

        @test length(eqs_func) == length(eqs_manual)
        @test all([Symbol(eqs_manual[i]) == Symbol(eqs_func[i]) for i in 1:length(eqs_func)])

        # Check the parameters
        ps = parameters(test)
        @test length(ps) == (length(parameters(mass)) +  length(parameters(damper)))

        # Expand equations
        eqs_expd = equations(expand_connections(test))
        # Check if the elements equations were maintaned
        @test all([eqs_expd[i] == eqs_elm[i] for i in 1:2])
        eqs_expd
        @test (0 ~ mass.power.e + sgn*damper.power.e) in eqs_expd
        @test (mass.power.f ~ damper.power.f) in eqs_expd
    end
end


@named x, y = mGY(damper, mass)


# @testset "mGY" begin
    # Default values
    u = π
    m = exp(1.0)
    # Function
    # @named mass = Mass(m=m, u=u)
    @named damper = Damper(c=m, u=u)
    # Naming
    @named test = mGY(mass, [sgn, damper])
    # Check naming
    @test test.name == name
    test = Junction1(mass, [sgn, damper], name=name)
    @test test.name == name
    # Check type junction
    @test get_bg_junction(test.power.e)[1] === j1
    @test get_bg_junction(test.power.f)[1] === j1
    # Manual
    @named power = Power(type=j1)
    # Test the power defaults
    @test getdefault(test.power.e) == 0.0
    @test getdefault(test.power.f) == 0.0
    # Test equation generated
    eqs_func = equations(test)
    if sgn == -1
        eqs_manual = [
            connect(mass.power, power),
            connect(power, damper.power),
        ]
    else
        eqs_manual = [
            connect(mass.power, power),
            connect(damper.power, power),
        ]
    end
    eqs_elm = []
    push!(eqs_elm, equations(compose(ODESystem(Equation[], t, [], []; name=:mass), mass))...)
    push!(eqs_elm, equations(compose(ODESystem(Equation[], t, [], []; name=:damper), damper))...)
    push!(eqs_manual, eqs_elm...)

    @test length(eqs_func) == length(eqs_manual)
    @test all([Symbol(eqs_manual[i]) == Symbol(eqs_func[i]) for i in 1:length(eqs_func)])

    # Check the parameters
    ps = parameters(test)
    @test length(ps) == (length(parameters(mass)) +  length(parameters(damper)))

    # Expand equations
    eqs_expd = equations(expand_connections(test))
    # Check if the elements equations were maintaned
    @test all([eqs_expd[i] == eqs_elm[i] for i in 1:2])
    eqs_expd
    @test (0 ~ mass.power.e + sgn*damper.power.e) in eqs_expd
    @test (mass.power.f ~ damper.power.f) in eqs_expd
# end


# @testset "mGY" begin
#     # Default values
#     u = π
#     m = exp(1.0)
#     sgn = -1
#     # Function
#     @named mass = Mass(m=m, u=u)
#     @named damper = Damper(c=m, u=u)
#     # Naming
#     @named test = Junction1(mass, [sgn, damper])
#     # Check naming
#     @test test.name == name
#     test = Junction1(mass, [sgn, damper], name=name)
#     @test test.name == name
#     # Check type junction
#     @test get_bg_junction(test.power.e)[1] === j1
#     @test get_bg_junction(test.power.f)[1] === j1
#     # Manual
#     @named power = Power(type=j1)
#     # Test the power defaults
#     @test getdefault(test.power.e) == 0.0
#     @test getdefault(test.power.f) == 0.0
#     # Test equation generated
#     eqs_func = equations(test)
#     eqs_manual = [
#         connect(mass.power, power),
#         connect(power, damper.power),
#     ]
#     eqs_elm = []
#     push!(eqs_elm, equations(compose(ODESystem(Equation[], t, [], []; name=:mass), mass))...)
#     push!(eqs_elm, equations(compose(ODESystem(Equation[], t, [], []; name=:damper), damper))...)
#     push!(eqs_manual, eqs_elm...)

#     @test length(eqs_func) == length(eqs_manual)
#     @test all([Symbol(eqs_manual[i]) == Symbol(eqs_func[i]) for i in 1:length(eqs_func)])

#     # Check the parameters
#     ps = parameters(test)
#     @test length(ps) == (length(parameters(mass)) + length(parameters(damper)))

#     # Expand equations
#     eqs_expd = equations(expand_connections(test))
#     # Check if the elements equations were maintaned
#     @test all([eqs_expd[i] == eqs_elm[i] for i in 1:2])
#     @test (0 ~ mass.power.e + sgn * damper.power.e) in eqs_expd
#     @test (mass.power.f ~ sgn * damper.power.f) in eqs_expd
# end