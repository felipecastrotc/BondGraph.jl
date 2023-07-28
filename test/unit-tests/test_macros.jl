using BondGraph
using Test
using Symbolics
using ModelingToolkit

using MacroTools
import BondGraph: t, oneport
import BondGraph: equalityeqs, sumvar, flatinput, isindependent
import ModelingToolkit: unwrap, isvariable, istree

@oneport

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