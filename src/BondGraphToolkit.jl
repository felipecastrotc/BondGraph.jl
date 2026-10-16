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
using ModelingToolkit: t_nounits as t, D_nounits as D, getname
using Symbolics: unwrap, get_variables, iscall, operation, arguments
using Latexify: latexify
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
export bondgraph, bondgraph_dot, bondgraph_svg
export ElectricalAdapter, TranslationalAdapter, TranslationalPositionAdapter,
    RotationalAdapter, HydraulicAdapter

"""
    PowerPort(; name)

The connector through which an element exchanges power: effort `e` and flow `f`. MTK
counts `f` positive into the component that owns the port.
"""
@connector function PowerPort(; name)
    vars = @variables begin
        e(t), [guess = 0.0]
        f(t), [connect = Flow, guess = 0.0]
    end
    System(Equation[], t, vars, []; name)
end

# Shared by @oneport and @twoport: make `fdef` an MTK component whose body ends with the
# vector of laws (long-form `function`). Its subsystems are `ports`, plus the signal
# connectors `input` (RealInput) and `output` (RealOutput) when the body defines them.
function element(fdef, ports, preamble, source)
    body = fdef.args[2]
    signals = [:(@isdefined($s) ? [$s] : []) for s in (:input, :output) if mentions(body, s)]
    body.args[end] = :($System(vcat($(body.args[end])), $t; name,
                               systems = [$(ports...); $(signals...)]))
    prepend!(body.args, preamble)
    esc(Expr(:macrocall, GlobalRef(ModelingToolkit, Symbol("@component")), source, fdef))
end
mentions(ex::Expr, s) = any(a -> mentions(a, s), ex.args)
mentions(ex, s) = ex === s

"""
    @oneport function Name(args...; name, kwargs...) ... [laws] end

Define a one-port element. The body can use `power.e` and `power.f`, the effort and
flow of the element's bond as drawn, and must end with the vector of laws. The macro
adds the `power` port and wraps the laws in a `System`, as an MTK `@component`; MTK
finds the unknowns and parameters itself. A signal connector the body defines as
`input` (`RealInput`) or `output` (`RealOutput`) becomes part of the element.
"""
macro oneport(fdef)
    element(fdef, [:power], [:(power = $PowerPort(name = :power))], __source__)
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
