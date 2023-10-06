using BondGraph
using Test
using Symbolics
using ModelingToolkit

import BondGraph: t, bg, tpgy
import BondGraph: check_bg_con, get_bg_connection_set!, add_idx, gen_tp_con!

import ModelingToolkit: ConnectionSet, rename
import ModelingToolkit: namespaced_var, get_connection_type, getname
import ModelingToolkit: generate_connection_set!


@testset "check_bg_con + get_bg_connection_set!" begin

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
    domain_csets = ConnectionSet[]
    sys = generate_connection_set!(connectionsets, domain_csets, sys, nothing, nothing)

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

@testset "add_idx" begin

    @parameters x=π
    @variables y(t)=exp(1)

    idx_x = 1
    idx_y = 2
    x_new = add_idx(x, idx_x)
    y_new = add_idx(y, idx_y)

    # Check the variable name
    @test getname(x_new) == Symbol("x" * string(idx_x))
    @test getname(y_new) == Symbol("y" * string(idx_y))

    # Check if the default value did not change
    @test Symbolics.getdefaultval(x_new) == π
    @test Symbolics.getdefaultval(y_new) == exp(1)

end

@testset "gen_tp_con!" begin
    # Dummy TP element
    # Generate the input and output connections
    @named pin = Power(type = tpgy)
    @named pout = Power(type = tpgy)

    sys_incomplete = ODESystem(Equation[], t, [], []; name = :s)
    sys_incomplete_pout = compose(sys_incomplete, pin)
    sys = compose(sys_incomplete_pout, pout)

    # Check if the system evaluates if it has pin and pout properties
    @test_throws AssertionError gen_tp_con!([], sys_incomplete, [])
    @test_throws AssertionError gen_tp_con!([], sys_incomplete_pout, [])

    # It should only generate connection equations if any subsys is passed
    subsys = []
    aux_eqs = []
    gen_tp_con!(aux_eqs, sys, subsys)
    @test length(aux_eqs) == 0

    # By default the first argument is left connect and the second argument 
    # is the connection on the right. 1stArg ---> TF ---> 2ndArg
    @named dummy = Mass()

    # Check the output for one ODE system passed
    subsys = dummy
    aux_eqs = []
    gen_tp_con!(aux_eqs, sys, subsys)
    @test aux_eqs == [connect(dummy.power, sys.pin)]

    # Check the output for one ODE system passed, but on the right connection
    subsys = [nothing, dummy]
    aux_eqs = []
    gen_tp_con!(aux_eqs, sys, subsys)
    @test aux_eqs == [connect(sys.pout, dummy.power,)]

    # Check the output for two ODE system passed
    subsys = [dummy, dummy]
    aux_eqs = []
    gen_tp_con!(aux_eqs, sys, subsys)
    @test aux_eqs == [connect(dummy.power, sys.pin), connect(sys.pout, dummy.power,)]
end