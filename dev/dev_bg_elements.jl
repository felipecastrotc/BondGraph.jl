using BondGraph, Plots, Symbolics.Latexify, DifferentialEquations
using ModelingToolkit

import BondGraph: t, D, Mass, Spring, Damper, Spring3
import ModelingToolkit: isvariable, istree, unwrap

# =============================================================================
# Elements


function Mass(; name, m = 1.0, u = 0.0)
    e, f = Power(flow = u)
    ps = @parameters I = m
    @variables p(t)

    eqs = [
        # It(e) ~ f*I
        # D(p) ~ e,
        # p ~ I * f,
        # f ~ p/I,
        D(f) ~ e / I,
    ]
    ODESystem(eqs, t, [e, f, p], ps; name = name)
    # ODESystem(eqs, t, [e, f], ps; name = name)
    # ODESystem(eqs, t, [], ps; name = name)
    # extend(ODESystem(eqs, t, [p], ps; name = name), power)
end

function Spring(; name, k = 10, x = 0.0)
    e, f = Power()
    @variables q(t) = x
    ps = @parameters C = 1 / k

    eqs = [
        e ~ q / C
        D(q) ~ f
        # D(e) ~ f/C
    ]
    ODESystem(eqs, t, [e, f, q], ps; name = name)
end

function Spring3(; name, k = 10, x = 0.0)
    e, f = Power()
    @variables q(t) = x
    ps = @parameters C = 1 / k

    eqs = [
        e ~ q^3 / C
        D(q) ~ f
        # D(e) ~ f/C
    ]
    ODESystem(eqs, t, [e, f, q], ps; name = name)
end

function Damper(; name, c = 10)
    e, f = Power()

    ps = @parameters R = c
    eqs = [e ~ f * R]
    ODESystem(eqs, t, [e, f], ps; name = name)
end
