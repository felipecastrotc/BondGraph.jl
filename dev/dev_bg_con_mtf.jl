using ModelingToolkit
using Symbolics

@variables t
D = Differential(t)

struct bg end
struct bgeffort end
struct bgflow end
struct op end
struct tpgy end
struct tptf end
struct j0 end
struct j1 end

get_connection_type(s) = getmetadata(
    Symbolics.unwrap(s),
    ModelingToolkit.VariableConnectType,
    ModelingToolkit.Equality,
)

# How the connect is parsed
# https://github.com/JuliaSymbolics/Symbolics.jl/blob/b7cae533e0a2aae8d7e227b2e0db4e25e6cb8d09/src/variable.jl
# _parse_vars
@macroexpand @variables e(t) [connect = :a, unit = :b]
# It adds metadata to the variable being created the connect unit and options 

# One idea is to set the metadata after like this
@variables e(t) [connect = bg]

e = Symbolics.wrap(setmetadata(Symbolics.unwrap(e), bg, j1))

getmetadata(e, bg)
get_connection_type(e)

# BG functions

set_bg_metadata(s, type) = Symbolics.wrap(setmetadata(Symbolics.unwrap(s), bg, type))
get_bg_junction(s) = getmetadata(Symbolics.unwrap(s), bg)
update_mtk_con(s, type) = Symbolics.wrap(
    setmetadata(Symbolics.unwrap(s), ModelingToolkit.VariableConnectType, type),
)

# Convert the bg metadata to comply to the modelingtoolkit connection types
function convert_bg_metadata(cset)

    i = findall(x -> x.isouter, cset.set)
    i = length(i) > 0 ? i[1] : 1

    cfg = get_bg_junction(cset.set[i].v)

    # get_bg_junction(cset.set[1].sys.sys.e)

    # get_bg_junction(mj.power.e)

    # cset.set[2].isouter

    T = ModelingToolkit.ConnectionElement
    newconnectionset = T[]

    vartype = nothing
    if (cfg[1] === j1 && cfg[2] === bgeffort) || (cfg[1] === j0 && cfg[2] === bgflow)
        vartype = Flow
    elseif (cfg[1] === j1 && cfg[2] === bgflow) || (cfg[1] === j0 && cfg[2] === bgeffort)
        vartype = :equality
    end

    for ele in cset.set
        namespace = ele.sys

        v = deepcopy(ele.v)
        v = update_mtk_con(v, vartype)

        isouter = ele.isouter

        push!(newconnectionset, T(namespace, v, isouter))
    end

    return ModelingToolkit.ConnectionSet(newconnectionset)
end

@connector function Power(; name, effort = 0.0, flow = 0.0, type = op)
    sts = @variables e(t) = effort [connect = bg] f(t) = flow [connect = bg]
    sts = set_bg_metadata.(sts, [[type, bgeffort], [type, bgflow]])
    ODESystem(Equation[], t, sts, []; name = name)
end

function Mass(; name, m = 1.0, u = 0.0)
    @named power = Power(flow = u)
    ps = @parameters I = m

    eqs = [D(power.f) ~ power.e / I]
    compose(ODESystem(eqs, t, [], ps; name = name), power)
end

function Spring(; name, k = 10, x = 0.0)
    @named power = Power()

    @variables q(t) = x
    ps = @parameters C = 1 / k

    eqs = [
        power.e ~ q / C
        D(q) ~ power.f
    ]
    compose(ODESystem(eqs, t, [q], ps; name = name), power)
end

function Damper(; name, c = 10, u = 1.0)
    @named power = Power(flow = u)

    ps = @parameters R = c
    eqs = [power.e ~ power.f * R]

    compose(ODESystem(eqs, t, [], ps; name = name), power)
end

function Junction0(ps...; name, couple = true)

    @named power = Power(type = j0)
    # Get connections
    ps = collect(Base.Flatten([ps]))

    # Get signs and subsystems 
    subsys, signs = flatinput(ps)

    eqs = Equation[]
    for (s, sign) in zip(subsys, signs)
        if sign == 1
            push!(eqs, connect(s.power, power))
        elseif sign == -1
            push!(eqs, connect(power, s.power))
        end
    end

    # Build subsystem
    sys = ODESystem(eqs, t, [], [], name = name)
    compose(sys, power, subsys...)
end

function Junction1(ps...; name, couple = true)

    @named power = Power(type = j1)
    # Get connections
    ps = collect(Base.Flatten([ps]))

    # Get signs and subsystems 
    subsys, signs = flatinput(ps)

    eqs = Equation[]
    for (s, sign) in zip(subsys, signs)
        if sign == 1
            push!(eqs, connect(s.power, power))
        elseif sign == -1
            push!(eqs, connect(power, s.power))
        end
    end

    # Build subsystem
    sys = ODESystem(eqs, t, [], [], name = name)
    compose(sys, power, subsys...)
end

function Se(expr; name)

    @named power = Power(type = op)

    eqs = [power.e ~ expr]

    sts = []

    if ModelingToolkit.isvariable(expr)
        ps = [expr]
    elseif ModelingToolkit.istree(Symbolics.unwrap(expr))
        vars = collect(Set(ModelingToolkit.get_variables(g)))
        sts = vcat(sts, filter(x -> ~isindependent(Num(x)), vars))
        ps = filter(x -> isindependent(Num(x)), vars)
    else
        ps = []
    end

    compose(ODESystem(eqs, t, sts, ps; name = name), power)
end

# -----------------------------------------------------------------------------
# Setting equations by functions

@named rm = Mass(m = 1.0)
@named rk = Spring(k = 10.0, x = 0.1)
@named rf = Damper(c = 1.0)

@named pj = Mass(m = 1.0)
@named pf = Damper(c = 0.01)
@variables T
@named pT = Se(3)

# -----------------------------------------------------------------------------
# Using Junctions and mTR
# Bond graph example 

# Gear radius
@variables R
@variables g(t)
R = GlobalScope(R)

@named psub = Junction1(pT, pj, pf)
@named rsub = Junction1(rm, rk, rf)

subsys = [psub, rsub]
c = subsys isa ODESystem ? [subsys, nothing] : collect(subsys)

pos = .!isnothing.(c)
c = c[pos]

r = R + g + 10
# Generate the in and out connection
@named pin = Power(type = tptf)
@named pout = Power(type = tptf)

# Alias for simpler view about the mGY port
e₁, f₁ = pin.e, pin.f
e₂, f₂ = pout.e, pout.f

eqs = [f₁ * r ~ f₂, e₂ * r ~ e₁]

if ModelingToolkit.isvariable(r)
    sts, ps = [], [r]
elseif ModelingToolkit.istree(ModelingToolkit.unwrap(r))
    sts = []
    ps = collect(Set(ModelingToolkit.get_variables(r)))
else
    sts, ps = [], []
end

push!(eqs, connect(c[1].power, pin), connect(pout, c[2].power))
eqs

name = :tf
# sys = ODESystem(eqs, t, sts, ps; name = name)
tf = compose(ODESystem(eqs, t, sts, ps; name = name), pin, pout)

tf = compose(sys, pin, pout)
equations(tf)

@named psys = ODESystem([connect(psub.power, tf.pin), connect(tf.pout, rsub.power)], t)
mdl = compose(psys, tf, psub, rsub)
mdl = compose(tf, psub, rsub)
equations(mdl)

generate_graph(mdl)

equations(mdl)
emdl = expand_connections(mdl)

equations(emdl)
generate_graph(mdl)

equations(alias_elimination(emdl))
equations(structural_simplify(emdl))

@named red = reducedobs(structural_simplify(emdl))

equations(red)
