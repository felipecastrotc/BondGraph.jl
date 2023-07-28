using BondGraph
using Test
using Symbolics
using ModelingToolkit

import BondGraph: t, bg
import BondGraph: check_bg_con, get_bg_connection_set!, get_var, add_idx
import BondGraph: gen_tp_con!

import ModelingToolkit: ConnectionSet, rename
import ModelingToolkit: namespaced_var, get_connection_type, getname
import ModelingToolkit: generate_connection_set!


@testset "check_bg_con + get_bg_connection_set!"

    @named mass = Mass()
    @named damper = Damper()

    # Manual Junction 1
    @named power = Power(type = j1)
    eqs = [
        connect(mass.power, power),
        connect(damper.power, power)
    ]
    sys = compose(ODESystem(eqs, t, [], [], name = :sys), power, mass, damper)

    # Based on the generate_connection_set function from MTK library
    connectionsets = ConnectionSet[]
    sys = generate_connection_set!(connectionsets, sys, nothing, nothing)

    # The number of connections must be the double of equations due to effort 
    # and flow variables
    @test length(connectionsets) == length(eqs)*2
    # Check if the connection all connection elements are identified as bg type
    for con in connectionsets
        @test check_bg_con(con)
    end

    # TODO: Increase the complexity of the system to add non-bg connections 
    @test length(get_bg_connection_set!(connectionsets)) == length(eqs)*2

end