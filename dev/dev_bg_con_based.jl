using ModelingToolkit
using Symbolics
import Symbolics: unwrap, wrap

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

set_bg_metadata(s, type) = wrap(setmetadata(unwrap(s), bg, type))
get_bg_junction(s) = getmetadata(unwrap(s), bg)
update_mtk_con(s, type) =
    wrap(setmetadata(unwrap(s), ModelingToolkit.VariableConnectType, type))

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

# Update this function to handle the bond graph type
function ModelingToolkit.generate_connection_equations_and_stream_connections(
    csets::AbstractVector{<:ModelingToolkit.ConnectionSet},
)
    eqs = Equation[]
    stream_connections = ModelingToolkit.ConnectionSet[]

    for cset in csets
        v = cset.set[1].v
        if hasmetadata(v, Symbolics.GetindexParent)
            v = getparent(v)
        end
        vtype = get_connection_type(v)
        if vtype === Stream
            push!(stream_connections, cset)
            continue
        elseif vtype === Flow
            rhs = 0
            for ele in cset.set
                v = ModelingToolkit.namespaced_var(ele)
                rhs += ele.isouter ? -v : v
            end
            push!(eqs, 0 ~ rhs)
        elseif vtype === bg
            ncset = convert_bg_metadata(cset)
            neqs, _ = ModelingToolkit.generate_connection_equations_and_stream_connections([
                ncset,
            ])
            push!(eqs, neqs...)
        else # Equality
            base = ModelingToolkit.namespaced_var(cset.set[1])
            for i = 2:length(cset.set)
                v = ModelingToolkit.namespaced_var(cset.set[i])
                push!(eqs, base ~ v)
            end
        end
    end
    eqs, stream_connections
end

function Base.merge(csets::AbstractVector{<:ModelingToolkit.ConnectionSet})
    mcsets = ModelingToolkit.ConnectionSet[]
    ele2idx = Dict{ModelingToolkit.ConnectionElement,Int}()
    cacheset = Set{ModelingToolkit.ConnectionElement}()

    for cset in csets
        if get_connection_type(cset.set[1].v) != bg
            idx = nothing
            for e in cset.set
                idx = get(ele2idx, e, nothing)
                idx !== nothing && break
            end
            if idx === nothing
                push!(mcsets, cset)
                for e in cset.set
                    ele2idx[e] = length(mcsets)
                end
            else
                for e in mcsets[idx].set
                    push!(cacheset, e)
                end
                for e in cset.set
                    push!(cacheset, e)
                end
                empty!(mcsets[idx].set)
                for e in cacheset
                    ele2idx[e] = idx
                    push!(mcsets[idx].set, e)
                end
                empty!(cacheset)
            end
        end
    end

    # bg
    str2con = name2sys(csets)
    if length(str2con) > 1
        al = csets2adjlist(csets)
        push!(mcsets, adjlist2csets(al, str2con)...)
    end
    mcsets
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

function Damper(; name, c = 10, u = 1.0)
    @named power = Power(flow = u)

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

function Se(expr; name)

    @named power = Power(type = op)

    eqs = [power.e ~ expr]

    sts = []
    if ModelingToolkit.isvariable(expr)
        ps = [expr]
    elseif ModelingToolkit.istree(ModelingToolkit.unwrap(expr))
        vars = collect(Set(ModelingToolkit.get_variables(g)))
        sts = vcat(sts, filter(x -> ~isindependent(Num(x)), vars))
        ps = filter(x -> isindependent(Num(x)), vars)
    else
        ps = []
    end

    compose(ODESystem(eqs, t, sts, ps; name = name), power)
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

@named m = Mass()
@named d = Damper()
@named d2 = Damper()
@named s = Spring()

@named b1 = Junction1(m, d)
@named b2 = Junction1(m, d)
@named b3 = Junction1(m, d)

# Huge information on how it is organized
# https://github.com/SciML/ModelingToolkit.jl/blob/28b36c8f9932acfc6c54737c51437d35b192b0d7/test/stream_connectors.jl
@named sys = ODESystem([connect(b1.power, b2.power)], t)
sys = expand_connections(sys)

nsys, csets = ModelingToolkit.generate_connection_set(sys)

ceqs, instream_csets =
    ModelingToolkit.generate_connection_equations_and_stream_connections(csets)

equations(sys)
mdl = compose(sys, b1, b2)
mdl = alias_elimination(mdl)
equations(mdl)

@named mdl = reducedobs(structural_simplify(mdl))
equations(mdl)

# Simple MSD
@named sys = Junction1(m, s, d, couple = false)
equations(sys)
equations(alias_elimination(sys))


# 2DOF
@named sd = Junction1(d)
@named mj = Junction1(m)
@named b1 = Junction0(mj, sd)

@named b0 = Junction1(m, d)

@named sys = ODESystem([connect(b0.power, b1.power)], t)
mdl = compose(sys, b1, b0)
mdl = expand_connections(mdl)
equations(mdl)
equations(alias_elimination(mdl))
@named mdl = reducedobs(structural_simplify(mdl))
eqs = equations(mdl)

# nDOF
@named sd = Junction1(d, d2)
@named mj = Junction1(m)

@named b1 = Junction0(mj, sd)
@named b2 = Junction0(mj, sd)
equations(b2)

@named b0 = Junction1(m, d)

# @named psys = ODESystem([connect(b2.power, b1.mj.power), connect(b1.power, b0.power)], t)
# mdl = compose(psys, b1, b0, b2)
# @named psys = ODESystem([connect(b2.power, b1.mj.power)], t)
@named psys = ODESystem([connect(b1.mj.power, b2.power)], t)
mdl = compose(psys, b1, b2)

equations(mdl)
emdl = expand_connections(mdl)

equations(emdl)

equations(alias_elimination(emdl))
equations(structural_simplify(emdl))


@named psys = ODESystem([connect(b0.power, b1.power)], t)
# sys = compose(psys, b1, b0,b2)
sys = compose(psys, b1, b0)

# ==========================================================
# MSD test for signal

@named F = Se(10)

@named b = Junction1([1, m], [1, s], [1, d], F)

mdl = b


# ==========================================================
# New test with multiple outputs for pump leakage

@named lek = Junction1(d)

@named suc = Junction1(m)
@named pm = Junction0()

@named imp = Junction1(m)

@named val = Junction1(m)
@named pj = Junction0()


cons = [
    connect(suc.power, pm.power),
    connect(lek.power, pm.power),
    connect(pm.power, imp.power),
    connect(imp.power, pj.power),
    connect(pj.power, lek.power),
    connect(pj.power, val.power),
]
@named psys = ODESystem(cons, t)
mdl = compose(psys, lek, suc, pm, imp, pj, val)
