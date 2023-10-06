using BondGraph
using Test
using Symbolics
using ModelingToolkit

import BondGraph: t
import BondGraph: equalityeqs, sumvar, flatinput, isindependent
import ModelingToolkit: unwrap, isvariable, istree

@testset "equalityeqs" begin

    @variables x(t), y(t)
    @parameters a, b

    vars = [a, b, x, y]
    eqs = equalityeqs(vars)
    # Check the numbers of equations
    @test length(eqs) == length(vars) - 1
    # Check if eqs is a subset of the complete combination set of equations
    # Get the combinations
    M = collect(Iterators.product(vars, vars))
    # Remove the diagonal elements and flatten the elements
    N = [M[i, j] for i in 1:size(M, 1) for j in 1:size(M, 2) if i != j]
    O = [x ~ y for (x, y) in N]
    neqs = sum([any([eq === c for eq in eqs]) for c in O])
    # The number neqs should be equal to the eqs length
    @test neqs == length(eqs)

    # Check the empty array case
    @test [] == equalityeqs([])

end

@testset "sumvar" begin

    @variables x(t), y(t)
    @parameters a, b

    vars = [a, b, x, y]
    eq = sumvar(vars)
    # Check the equation
    @test (0 ~ sum(vars)) == (0 ~ eq.rhs - eq.lhs)

    # Check the empty array case
    @test [] == equalityeqs([])
end

@testset "flatinput" begin

    function check_flatinput_ps(ps, truth_sys, truth_signs)
        subsys, signs = flatinput(ps)
        # It should not mute the subsystems
        for sub in subsys
            @test any([sub == truth for truth in truth_sys])
        end
        # The sign array should have the same size of subsys
        @test length(signs) == length(subsys)
        # Check the signs array
        @test isa(signs, Vector{Int})
        @test signs == truth_signs
    end

    @named mass = Mass()
    @named damper = Damper()

    # Check emty inputs
    subsys, signs = flatinput([])
    @test subsys == ModelingToolkit.AbstractSystem[]
    @test signs == Int[]

    # Check inputs without sign, the signs output it should default to +1
    ps = [mass, damper]
    check_flatinput_ps(ps, ps, ones(Int, length(ps)))

    # Check heterogeneous ps with different sign deffinition
    ps = [[1, mass], [2.0, damper], [-10, mass], [-1.0, damper], mass];
    truth_sys = [mass, damper, mass, damper, mass];
    truth_signs = [1, 1, -1, -1, 1]
    check_flatinput_ps(ps, truth_sys, truth_signs)
    # Check heterogeneous ps with different sign deffinition
    ps = [["1", mass], damper];
    @test_throws DomainError flatinput(ps)

end

@testset "isindependent" begin
    @variables x(t), b
    @parameters a, y(t)

    @named power = Power()

    @test isindependent(a)
    @test isindependent(b)

    @test !isindependent(x)
    @test !isindependent(y)
end
