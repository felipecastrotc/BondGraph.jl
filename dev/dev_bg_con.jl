using ModelingToolkit
using Symbolics
using GraphRecipes
using Plots

import Symbolics: unwrap, wrap

@variables t
D = Differential(t)

struct bg end
struct bgeffort end
struct bgflow end
struct op end
struct j0 end
struct j1 end

get_connection_type(s) = getmetadata(Symbolics.unwrap(s), ModelingToolkit.VariableConnectType, ModelingToolkit.Equality)

# How the connect is parsed
# https://github.com/JuliaSymbolics/Symbolics.jl/blob/b7cae533e0a2aae8d7e227b2e0db4e25e6cb8d09/src/variable.jl
# _parse_vars
@macroexpand @variables e(t) [connect = :a, unit = :b]
# It adds metadata to the variable being created the connect unit and options 

# One idea is to set the metadata after like this
@variables e(t) [connect = bg]
e = Symbolics.wrap(setmetadata(Symbolics.unwrap(e), bg, j1))

# BG functions

set_bg_metadata(s, type) = wrap(setmetadata(unwrap(s), bg, type))
get_bg_junction(s) = getmetadata(unwrap(s), bg)
update_mtk_con(s, type) = wrap(setmetadata(unwrap(s), ModelingToolkit.VariableConnectType, type))

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
        
        isouter  = ele.isouter

        push!(newconnectionset, T(namespace, v, isouter))
    end

    return ModelingToolkit.ConnectionSet(newconnectionset)
end

@connector function Power(; name, effort=0.0, flow=0.0, type=op)
    sts = @variables e(t) = effort [connect = bg] f(t) = flow [connect = bg]
    sts = set_bg_metadata.(sts, [[type, bgeffort], [type, bgflow]])
    ODESystem(Equation[], t, sts, []; name=name)
end

function Mass(; name, m = 1.0, u = 0.0)
    @named power = Power(flow=u)    
    ps = @parameters I = m

    eqs = [
        D(power.f) ~ power.e / I,
    ]
    compose(ODESystem(eqs, t, [], ps; name = name), power)
end

function Damper(; name, c = 10, u=1.0)
    @named power = Power(flow=u)
    
    ps = @parameters R = c
    eqs = [power.e ~ power.f * R]

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


function Junction0(ps...; name, couple=true)
    @named power = Power(type=j0)

    # couple = true
    # ps = [mj, sd]
    # name = :oi
    # Get connections
    ps = collect(Base.Flatten([ps]))
    
    # Split connections
    oneport = getoneport(ps)
    multiport = getmultiport(ps)

    eqs = Equation[]
    if length(oneport) > 0
        sgns, e = couple ? ([1], power.e) : ([], [])
        f = couple ? power.f : 0.0
        # Σ efforts
        push!(eqs, 0 ~ sumvar(oneport, :f) + f)
        # f₁ = f₂ = f₃
        push!(eqs, equalityeqs(oneport, :e, sgns, e)...)
    end
    if length(multiport) > 0
        for m in multiport
            push!(eqs, connect(m.power, power))
        end
    end

    # Build subsystem
    sys = ODESystem(eqs, t, [], [], name = name)
    compose(sys, power, oneport..., multiport...)
end

function Junction1(ps...; name, couple=true)

    @named power = Power(type=j1)

    # couple = true
    # ps = [mj, sd]
    # name = :oi
    # Get connections
    ps = collect(Base.Flatten([ps]))
    
    # Split connections
    oneport = getoneport(ps)
    multiport = getmultiport(ps)

    eqs = Equation[]
    if length(oneport) > 0
        sgns, f = couple ? ([1], power.f) : ([], [])
        e = couple ? power.e : 0.0
        # Σ efforts
        push!(eqs, 0 ~ sumvar(oneport, :e) + e)
        # f₁ = f₂ = f₃
        push!(eqs, equalityeqs(oneport, :f, sgns, f)...)
    end
    if length(multiport) > 0
        for m in multiport
            push!(eqs, connect(m.power, power))
        end
    end

    # Build subsystem
    sys = ODESystem(eqs, t, [], [], name = name)
    compose(sys, power, oneport..., multiport...)
end

# ndof
@named m = Mass()
@named d = Damper()
@named d2 = Damper()
@named s = Spring()

@named sd = Junction1(d, s)
@named mj = Junction1(m)

@named b1 = Junction0(mj, sd)
@named b2 = Junction0(mj, sd)
equations(b2)

@named b0 = Junction1(m, d)

@named psys = ODESystem([connect(b1.mj.power, b2.power), connect(b1.power, b0.power)], t)
mdl = compose(psys, b1, b0, b2)
# @named psys = ODESystem([connect(b2.power, b1.mj.power)], t)
# @named psys = ODESystem([connect(b1.mj.power, b2.power)], t)
# mdl = compose(psys, b1, b2)
generate_graph(mdl)

equations(mdl)
emdl = expand_connections(mdl)

equations(emdl)
equations(alias_elimination(emdl))
sys = structural_simplify(emdl)
equations(sys)

@named sys = reducedobs(sys)
equations(sys)

# ==========================================================
# 1dof

@named b = Junction1(m, s, d)
emdl = expand_connections(b)
sys = structural_simplify(emdl)
@named sys = reducedobs(sys)
equations(sys)
generate_graph(b)

# ==========================================================
# 2dof

@named sd = Junction1(s, d)
@named mj = Junction1(m)
@named b1 = Junction0(mj, sd)

@named b0 = Junction1(m, s, d)

@named sys = ODESystem([connect(b0.power, b1.power)], t)
mdl = compose(sys, b1, b0)
generate_graph(mdl)

emdl = expand_connections(mdl)
equations(emdl)
@named emdl = reducedobs(structural_simplify(emdl))
equations(emdl)


# ==========================================================
# New test with multiple outputs for pump leakage

@named lek = Junction1(d)

@named suc = Junction1(d)
@named pm = Junction0()

@named imp = Junction1(m)

@named val = Junction1(d)
@named pj = Junction0()


cons = [connect(suc.power, pm.power), connect(lek.power, pm.power), connect(pm.power, imp.power), connect(imp.power, pj.power), connect(pj.power, lek.power), connect(pj.power, val.power)]
@named psys = ODESystem(cons, t)
mdl = compose(psys, lek, suc, pm, imp, pj, val)
generate_graph(mdl)

equations(mdl)
emdl = expand_connections(mdl)

equations(emdl)
equations(alias_elimination(emdl))
sys = structural_simplify(emdl)
equations(sys)

@named sys = reducedobs(sys)
equations(sys)
