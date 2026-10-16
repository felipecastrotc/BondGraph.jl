# 0- and 1-junctions: one signed component (Karnopp Eqs. 3.16-3.21, ADR 0001).
#
# A junction member is `x` (bond into the junction, sign +1) or `[-1, x]` (bond out of
# it), where `x` is an element or a `Symbol` naming a free port. The core has one
# PowerPort per bond. MTK measures each core port's flow `c` into the core, so the bond
# flow along its half-arrow is `σ c`:
#   - element member: `connect(x.power, port)` makes x.power.f = -c, and the element's
#     variable is the bond flow whatever the sign, so σ = -1;
#   - free port: `c` is the flow into this junction, so σ = sign.

_member(a::AbstractVector) = length(a) == 2 && a[1] in (1, -1) ? (Int(a[1]), a[2]) :
    throw(ArgumentError("junction member must be `x` or `[-1, x]` (or `[1, x]`), got $a"))
_member(a) = (1, a)

function junction_core(kind, signs, σs; name)
    ps = [PowerPort(name = Symbol(:p, k)) for k in eachindex(signs)]
    bond_f = [σ * p.f for (σ, p) in zip(σs, ps)]
    eqs = kind == 0 ?
        [[p.e ~ ps[1].e for p in ps[2:end]]; 0 ~ sum(signs .* bond_f)] :
        [[f ~ bond_f[1] for f in bond_f[2:end]]; 0 ~ sum(s * p.e for (s, p) in zip(signs, ps))]
    System(eqs, t; name, systems = ps)
end

function junction(kind, args; name)
    signs, xs = first.(_member.(args)), last.(_member.(args))
    σs = [x isa Symbol ? s : -1 for (s, x) in zip(signs, xs)]
    core = junction_core(kind, signs, σs; name = :core)
    ends = [x isa Symbol ? PowerPort(name = x) : x for x in xs]
    cons = [connect(x isa Symbol ? sys : sys.power, getproperty(core, Symbol(:p, k)))
            for (k, (x, sys)) in enumerate(zip(xs, ends))]
    System(cons, t; name, systems = [core; ends])
end

"""
    Junction0(members...; name)

0-junction: every bond shares one effort and `Σ ± fᵢ = 0`, + for bonds into the
junction. A member is an element `x` (its bond points into the junction), `[-1, x]` (its
bond points out of the junction, towards `x`), or a free port `:name` / `[-1, :name]`
that other junctions and two-ports bond to with `connect`. A plain `connect` of free
ports is itself a 0-junction.

```julia
@named rc = Junction0([-1, R], [-1, C])   # 0 ⇀ R, 0 ⇀ C
```
"""
Junction0(args...; name) = junction(0, args; name)

"""
    Junction1(members...; name)

1-junction: every bond shares one flow and `Σ ± eᵢ = 0`, + for bonds into the junction.
Members as in [`Junction0`](@ref).

```julia
@named m = Junction1(se, [-1, R], [-1, C])   # Se ⇀ 1, 1 ⇀ R, 1 ⇀ C
```
"""
Junction1(args...; name) = junction(1, args; name)
