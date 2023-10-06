using BondGraph
using BondGraph: t, D, bg, j0, j1
import BondGraph: get_bg_junction
using ModelingToolkit
using ModelingToolkit: getdefault
using DifferentialEquations
using Test
using Symbolics

@testset "oneport-macro" begin
    # Create the cubic spring model
    @oneport function Spring3(; name, k = 1.0, x = 0.0)

        ## Initialize the displacement state with its initial value `x`
        @variables q(t) = x
        ## Set the equation parameters
        @parameters C = 1 / k

        ## Define Spring3 element equations with the non-lineariyy of the cubic spring
        [
            power.e ~ q^3 / C
            D(q) ~ power.f
        ]
    end

    # Create the element
    @parameters a, b, c
    @named test = Spring3(k=a, x=b)

    # Define the system manually
    @variables q(t) = b
    @parameters C = 1/a
    @named power = Power()
    # Test the equations
    @test length(equations(test)) == 2
    # Test the equations
    @test (power.e ~ q^3/C)  == equations(test)[1]
    @test (D(q) ~ power.f)  == equations(test)[2]
    # Test variables for default
    @test getdefault(q) === getdefault(test.q)
    @test string(getdefault(C)) == string(getdefault(test.C))
end