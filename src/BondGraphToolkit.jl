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

export PowerPort, @oneport
export Se, Sf, Resistance, Compliance, Inertia
export Junction0, Junction1

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

"""
    @oneport function Name(args...; name, kwargs...) ... [laws] end

Define a one-port element. The body can use `power.e` and `power.f`, the effort and
flow of the element's bond as drawn, and must end with the vector of laws. The macro
adds the `power` port and wraps the laws in a `System`, as an MTK `@component`; MTK
finds the unknowns and parameters itself.
"""
macro oneport(fdef)
    body = fdef.args[2]
    body.args[end] = :($System(vcat($(body.args[end])), $t; name, systems = [power]))
    pushfirst!(body.args, :(power = $PowerPort(name = :power)))
    esc(Expr(:macrocall, GlobalRef(ModelingToolkit, Symbol("@component")), __source__, fdef))
end

include("elements.jl")
include("junctions.jl")

end # module
