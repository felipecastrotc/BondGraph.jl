# 0- and 1-junctions: one signed component (Karnopp Eqs. 3.16-3.21, ADRs 0001, 0002).
#
# A junction member is an element `x` (its bond points from the junction into x, through
# `x.power`), `x => :pout` (from x into the junction, through `x.pout`), or a free port
# that other junctions and two-ports bond to with `connect`: `:a` (the bond leaves the
# junction) or `:a => :in` (it enters). Each member's sign is its bond's half-arrow, +1
# into the junction. The core has one PowerPort per bond and MTK measures each core
# port's flow `c` into the core; for every kind of member the bond flow along its
# half-arrow is then `sign * c`.

const MEMBER_FORMS = "`x` (bond into element x), `x => :pout` (bond out of x), `:a` (free port, bond leaving the junction) or `:a => :in` (bond entering it)"

# A member as (sign, element or free-port name, port): port is :power/:pout, or :free.
function _member(a)
    a isa AbstractVector && throw(ArgumentError(
        "junction member `[±1, x]`: this form was removed; write $MEMBER_FORMS"))
    x, port = a isa Pair ? a : (a, a isa Symbol ? :out : :power)
    if x isa Symbol
        port in (:in, :out) || throw(ArgumentError("free port `$x => $port`: use `:$x` or `:$x => :in`"))
        return (port === :in ? 1 : -1, x, :free)
    end
    port in (:power, :pout) ||
        throw(ArgumentError("junction member `$(nameof(x)) => $port`: write $MEMBER_FORMS"))
    (port === :pout ? 1 : -1, x, port)
end

struct JunctionMembers end   # system metadata key

function junction_core(kind, signs; name)
    ps = [PowerPort(name = Symbol(:p, k)) for k in eachindex(signs)]
    bond_f = [s * p.f for (s, p) in zip(signs, ps)]
    eqs = kind == 0 ?
        [[p.e ~ ps[1].e for p in ps[2:end]]; 0 ~ sum(signs .* bond_f)] :
        [[f ~ bond_f[1] for f in bond_f[2:end]]; 0 ~ sum(s * p.e for (s, p) in zip(signs, ps))]
    System(eqs, t; name, systems = ps)
end

function junction(kind, args; name)
    members = _member.(args)
    signs = first.(members)
    core = junction_core(kind, signs; name = :core)
    ends = [port === :free ? PowerPort(name = x) : x for (_, x, port) in members]
    cons = [connect(port === :free ? sys : getproperty(sys, port), getproperty(core, Symbol(:p, k)))
            for (k, ((_, _, port), sys)) in enumerate(zip(members, ends))]
    # Recorded for `BondGraph`: each member's sign, name, and whether it is a free port.
    meta = [(s, nameof(e), port === :free) for ((s, _, port), e) in zip(members, ends)]
    System(cons, t; name, systems = [core; ends], metadata = Dict(JunctionMembers => (kind, meta)))
end

"""
    Junction0(members...; name)

0-junction: every bond shares one effort and `Σ ± fᵢ = 0`, + for bonds into the
junction. A member is an element `x` (its bond points from the junction into `x`),
`x => :pout` (its bond points out of `x`, into the junction), or a free port that other
junctions and two-ports bond to with `connect`: `:a` (the bond leaves the junction) or
`:a => :in` (it enters). A plain `connect` of element and free ports is itself a
0-junction.

```julia
@named rc = Junction0(sf => :pout, R, C)   # Sf ⇀ 0, 0 ⇀ R, 0 ⇀ C
```
"""
Junction0(args...; name) = junction(0, args; name)

"""
    Junction1(members...; name)

1-junction: every bond shares one flow and `Σ ± eᵢ = 0`, + for bonds into the junction.
Members as in [`Junction0`](@ref).

```julia
@named m = Junction1(se => :pout, R, C)   # Se ⇀ 1, 1 ⇀ R, 1 ⇀ C
```
"""
Junction1(args...; name) = junction(1, args; name)
