"""
    BondGraphToolkit

A bond-graph dialect of ModelingToolkit: elements, junctions and two-ports are ordinary
MTK systems, bonds are `connect`s between [`PowerPort`](@ref)s, and models simulate with
`mtkcompile` → `ODEProblem` → `solve`.

Sign convention (Karnopp, Margolis & Rosenberg 2012, §§2.2, 3.1, 3.3): a bond's
half-arrow is the positive direction of power; element laws are written in the
variables of their bond as drawn; at a junction the shared variable is unsigned and the
sum is signed, + for bonds into the junction and − for bonds out of it.
"""
module BondGraphToolkit

using Reexport
@reexport using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SymbolicIndexingInterface: getname
using Symbolics: unwrap, get_variables, iscall, operation, arguments
using Latexify: latexify
using Graphviz_jll: Graphviz_jll
import REPL
using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput
using ModelingToolkitStandardLibrary.Electrical: Pin
using ModelingToolkitStandardLibrary.Mechanical: Translational, TranslationalPosition, Rotational
using ModelingToolkitStandardLibrary.Hydraulic: IsothermalCompressible

export PowerPort, @oneport, @twoport
export Se, Sf, Resistance, Compliance, Inertia, De, Df
export Mass, Spring, Damper
export Junction0, Junction1
export TF, GY
export bg_equations, latexify
export bondgraph, bgplot, savebg
export ElectricalAdapter, TranslationalAdapter, TranslationalPositionAdapter,
    RotationalAdapter, HydraulicAdapter

"""
    PowerPort(; name)

The connector through which an element exchanges power: effort `e` and flow `f`. MTK
counts `f` positive into the component that owns the port. A bond is a `connect` from
an outgoing port (an element's `pout`, a two-port's `p2`, a junction's outgoing free
port) to an incoming one (an element's `power`, a two-port's `p1`, a junction's
incoming free port); see [Concepts and Conventions](@ref concepts).
"""
@connector function PowerPort(; name)
    vars = @variables begin
        e(t)
        f(t), [connect = Flow]
    end
    System(Equation[], t, vars, []; name)
end

"""
    bond_ports()

The ports and bond variables of an element: `power`, the incoming port (its
bond points into the element), `pout`, the outgoing port (its bond points out of it),
and the bond's effort `e` and flow `f` along its half-arrow, whichever port carries it.
Returns `(power, pout, e, f, eqs)`. An element connects one of the two ports; MTK
zeroes the flow of the other, so `f = power.f - pout.f` is the bond flow either way.
MTK leaves the unused port's effort free, so it is tied to the bond's: both ports share
one effort, and only the flow carries the direction.
"""
function bond_ports()
    @named power = PowerPort()
    @named pout = PowerPort()
    @variables e(t) f(t)
    (power, pout, e, f, [e ~ power.e, pout.e ~ power.e, f ~ power.f - pout.f])
end

# Shared by @oneport and @twoport: make `fdef` an MTK component whose body ends with the
# vector of laws (long-form `function`). Its subsystems are `ports`, plus the signal
# connectors `input` (RealInput) and `output` (RealOutput) when the body defines them;
# `extra` are equations added to the laws.
function element(fdef, ports, preamble, source; extra = :(Equation[]))
    body = fdef.args[2]
    signals = [:(@isdefined($s) ? [$s] : []) for s in (:input, :output) if mentions(body, s)]
    body.args[end] = :($System(vcat($extra, $(body.args[end])), $t; name,
                               systems = [$(ports...); $(signals...)]))
    prepend!(body.args, preamble)
    esc(Expr(:macrocall, GlobalRef(ModelingToolkit, Symbol("@component")), source, fdef))
end
mentions(ex::Expr, s) = any(a -> mentions(a, s), ex.args)
mentions(ex, s) = ex === s

"""
    @oneport function Name(args...; name, kwargs...) ... [laws] end

Define a one-port element. The body can use `e` and `f`, the effort and flow of the
element's bond along its half-arrow, and must end with the vector of laws. The macro
adds the ports `power` (bond into the element) and `pout` (bond out of it, see
[`bond_ports`](@ref BondGraphToolkit.bond_ports)) and wraps the laws in a `System`, as an
MTK `@component`; MTK finds the unknowns and parameters itself. A signal connector the
body defines as `input` (`RealInput`) or `output` (`RealOutput`) becomes part of the
element.
"""
macro oneport(fdef)
    eqs = gensym(:bond)
    element(fdef, [:power, :pout], [:((power, pout, e, f, $eqs) = $bond_ports())],
            __source__; extra = eqs)
end

"""
    @twoport function Name(args...; name, kwargs...) ... [laws] end

Define a two-port element: port `p1` is the bond into it, port `p2` the bond out of it.
The body can use the bond variables `e1, f1` (port 1) and `e2, f2` (port 2), each flow
positive along its bond's half-arrow, and must end with the vector of laws, written as
in the books (e.g. Karnopp Eq. 3.9). The macro adds the ports and wraps the laws in a
`System`, as an MTK `@component`. Signal connectors as for [`@oneport`](@ref).

MTK counts every port's flow positive into the component, so `f2 = -p2.f`; write the
laws in `e1, f1, e2, f2` and this sign is handled for you.
"""
macro twoport(fdef)
    element(fdef, [:p1, :p2],
            [:(p1 = $PowerPort(name = :p1)), :(p2 = $PowerPort(name = :p2)),
             :(e1 = p1.e), :(f1 = p1.f), :(e2 = p2.e), :(f2 = -p2.f)], __source__)
end

include("elements.jl")
include("junctions.jl")
include("two_ports.jl")
include("equations.jl")
include("adapters.jl")
include("graph.jl")

end # module
