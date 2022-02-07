using Plots
using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Symbolics

@variables t
D = Differential(t)

?@connector


@connector function Power(;name, effort=0.0, flow=0.0)
    sts = @variables e(t)=effort f(t)=flow
    ODESystem(Equation[], t, sts, []; name=name)
end

function Mass(; name, m = 1.0, u = 0.)
    @named power = Power(flow=u)
    @unpack e, f = power
    ps = @parameters I=m
    eqs = [
            D(f) ~ e/I
        ]
    extend(ODESystem(eqs, t, [], ps; name=name), power)
end

function Spring(; name, k = 10, x = 0.)
    @named power = Power()
    @unpack e, f = power

    @variables q(t)=x
    ps = @parameters C=1/k
    eqs = [
            e ~ q/C
            D(q) ~ f
            # D(e) ~ f/C
        ]
    # extend(ODESystem(eqs, t, [], ps; name=name), power)
    extend(ODESystem(eqs, t, [q], ps; name=name), power)
    # V = e
    # I = f
    # ODESystem(Equation[], t, [q], ps; name)
end

function Damper(; name, c = 10)
    @named power = Power()
    @unpack e, f = power

    ps = @parameters R=c
    eqs = [
            e ~ f*R
        ]
    extend(ODESystem(eqs, t, [], ps; name=name), power)
end

function Junction1(ps...; name, subsys=[], couple=true)
    # Get connections
    con = collect(Base.Flatten([ps, subsys]))
    # println(typeof(con), "  ", typeof(ps))
    if couple
        @named power = Power()
        @unpack e, f = power
    
        # Σ efforts
        eqs = [e ~ sum(c -> c.e, con)]
        # f₁ = f₂ = f₃
        if length(con) > 1
            for i in 1:(length(con) - 1)
                push!(eqs, con[i].f ~ con[i + 1].f)
            end
        end
        push!(eqs, con[1].f ~ f)
        # Build subsystem
        compose(extend(ODESystem(eqs, t, [], []; name=name), power), ps...)
    else
        # Σ efforts
        eqs = [0 ~ sum(p -> p.e, con)]
        # f₁ = f₂ = f₃
        if length(con) > 1
            for i in 1:(length(con) - 1)
                push!(eqs, con[i].f ~ con[i + 1].f)
            end
        end
        # Build system
        compose(ODESystem(eqs, t, [], []; name=name), ps...)
    end
end

function Junction0(ps...; name, subsys=[], couple=true)
    # Get connections
    con = collect(Base.Flatten([ps, subsys]))
    println(typeof(con), "  ", typeof(ps))
    println(length(con))
    if couple
        @named power = Power()
        @unpack e, f = power
    
        # Σ flows
        eqs = [f ~ sum(c -> c.f, con)]
        # e₁ = e₂ = e₃
        if length(con) > 1
            for i in 1:(length(con) - 1)
                push!(eqs, con[i].e ~ con[i + 1].e)
            end
        end
        push!(eqs, con[1].e ~ e)
        # Build subsystem
        compose(extend(ODESystem(eqs, t, [], []; name=name), power), ps...)
    else
        # Σ flows
        eqs = [0 ~ sum(p -> p.f, con)]
        # e₁ = e₂ = e₃
        if length(con) > 1
            for i in 1:(length(con) - 1)
                push!(eqs, con[i].e ~ con[i + 1].e)
            end
        end
        # Build system
        compose(ODESystem(eqs, t, [], []; name=name), ps...)
    end
end

function Junction1(ps...; name, subsys=[], couple=true)
    # Get connections
    con = collect(Base.Flatten([ps, subsys]))
    # println(typeof(con), "  ", typeof(ps))
    if couple
        @named power = Power()
        @unpack e, f = power
    
        # Σ efforts
        eqs = [e ~ sum(c -> c.e, con)]
        # f₁ = f₂ = f₃
        if length(con) > 1
            for i in 1:(length(con) - 1)
                push!(eqs, con[i].f ~ con[i + 1].f)
            end
        end
        push!(eqs, con[1].f ~ f)
        # Build subsystem
        compose(extend(ODESystem(eqs, t, [], []; name=name), power), ps...)
    else
        # Σ efforts
        eqs = [0 ~ sum(p -> p.e, con)]
        # f₁ = f₂ = f₃
        if length(con) > 1
            for i in 1:(length(con) - 1)
                push!(eqs, con[i].f ~ con[i + 1].f)
            end
        end
        # Build system
        compose(ODESystem(eqs, t, [], []; name=name), ps...)
    end
end

function Junction0(ps...; name, subsys=[], couple=true)
    # Get connections
    con = collect(Base.Flatten([ps, subsys]))
    println(typeof(con), "  ", typeof(ps))
    println(length(con))
    if couple
        @named power = Power()
        @unpack e, f = power
    
        # Σ flows
        eqs = [f ~ sum(c -> c.f, con)]
        # e₁ = e₂ = e₃
        if length(con) > 1
            for i in 1:(length(con) - 1)
                push!(eqs, con[i].e ~ con[i + 1].e)
            end
        end
        push!(eqs, con[1].e ~ e)
        # Build subsystem
        compose(extend(ODESystem(eqs, t, [], []; name=name), power), ps...)
    else
        # Σ flows
        eqs = [0 ~ sum(p -> p.f, con)]
        # e₁ = e₂ = e₃
        if length(con) > 1
            for i in 1:(length(con) - 1)
                push!(eqs, con[i].e ~ con[i + 1].e)
            end
        end
        # Build system
        compose(ODESystem(eqs, t, [], []; name=name), ps...)
    end
end

function junction0(ps...)
    vcat([
        0 ~ sum(p -> p.f, ps)
        ],[
        ps[i].e ~ ps[i+1].e for i in 1:(length(ps) - 1)
    ])
end

function junction1(ps...)
    vcat([
        0 ~ sum(p -> p.e, ps)
        ],[
        ps[i].f ~ ps[i+1].f for i in 1:(length(ps) - 1)
    ])
end

# =============================================================================
# Definitions

m = 1.0
X = 1.
k = 10.
c = 1
center = 0.
@named mass = Mass(m=m)
@named spring = Spring(k=k, x=X)
@named damper = Damper(c=c)

# =============================================================================
# 1DOF

@named j1 = Junction1(spring, damper)

eqs = junction0(mass, j1)
eqs = junction1(mass, j1)

@named _model = ODESystem(eqs, t)
@named model = compose(_model, mass, j1)

sys = structural_simplify(model)
equations(sys)

prob = ODEProblem(sys, [], (0., 20.))
sol = solve(prob)

plot(sol)

# =============================================================================
# 2DOF

@named sd = Junction1(spring, damper)
@named bb = Junction1(mass, spring, damper)

eqs = junction0(mass, sd, bb)
_model = Junction0(mass, sd; name=:bt, subsys=[bb], couple=false)


@named _model = ODESystem(eqs, t)
@named model = compose(_model, spring, damper, mass, sd, bb)
@named model = compose(_model, sd, bb, mass)
@named model = compose(_model, bb)

sys = structural_simplify(model)
sys = structural_simplify(_model)
equations(sys)
states(sys)

prob = ODEProblem(sys, [1, 0, 0, 0], (0., 20.))
sol = solve(prob)

plot(sol)

# =============================================================================
# 4DOF
@named sd = Junction1(spring, damper)

S = [Junction0(sd, mass; name=:b3j0)];
push!(S, Junction1(mass; name=:b2j1, subsys=[S[1]]));
# push!(S, Junction0(sd; name=:b2j0, subsys=[S[2]]));
# push!(S, Junction1(mass; name=:b1j1, subsys=[S[3]]));
# push!(S, Junction0(sd; name=:b1j0, subsys=[S[4]]));
# _model = Junction1(mass, spring, damper; name=:model, subsys=[S[5]], couple=false)
_model = Junction0(sd; name=:b2j0, subsys=[S[2]], couple=false)
model = compose(_model, S...; name=:model)

model




?ModelingToolkit.connect
?@connector

equations(_model)
states(_model)

_model.systems[1].systems[1]
model.systems[1].systems[2]

typeof(model)


ModelingToolkit.get_systems(S[1])

ModelingToolkit.namespace_equations(S[1])
ModelingToolkit.namespace_equations(_model)

?ModelingToolkit.get_eqs
?ModelingToolkit.namespace_equations


equations(_model)
equations(model)
ModelingToolkit.get_eqs(_model)
[i.name for i in collect(_model.systems)]

S[2].eqs

[i.name for i in S]

ModelingToolkit.get_eqs(S[3])
ModelingToolkit.get_eqs(_model)
ModelingToolkit.namespace_equations(S[2])
ModelingToolkit.independent_variables(S[2])

eqs = ModelingToolkit.get_eqs(S[2])
O = ModelingToolkit.unwrap(eqs[1].rhs)

ModelingToolkit.istree(O)
args = ModelingToolkit.arguments(O)
x = ModelingToolkit.unwrap(args[1])

ModelingToolkit.getmetadata(x, SymScope, LocalScope())


ModelingToolkit.get_eqs(_model)
ModelingToolkit.get_eqs(S[2])
ModelingToolkit.get_eqs(S[1])

ModelingToolkit.get_eqs(sd)
ModelingToolkit.get_eqs(mass)
equations(S[1])
?equations



?ModelingToolkit.@set!
ModelingToolkit.get_systems(_model)

ModelingToolkit.get_systems(S[1])[2].name
ModelingToolkit.get_systems(S[1])[1].name
ModelingToolkit.get_systems(S[2])[1].name
states(model)

equations(model)


structural_simplify(S[1])


b2j1₊b3j0₊e(t) + b2j1₊mass₊e(t) ~ sd₊damper₊R*sd₊damper₊f(t) + sd₊spring₊q(t)*(sd₊spring₊C^-1)
b3j0₊mass₊e(t) ~ b3j0₊sd₊damper₊R*b3j0₊sd₊damper₊f(t) + b3j0₊sd₊spring₊q(t)*(b3j0₊sd₊spring₊C^-1)
?compose

b3j0₊f(t) ~ b3j0₊mass₊f(t) + b3j0₊sd₊damper₊f(t)

b2j1₊mass₊f(t) ~ b2j1₊b3j0₊f(t)


Differential(t)(b3j0₊sd₊spring₊q(t)) ~ b3j0₊sd₊damper₊f(t)
Differential(t)(sd₊spring₊q(t)) ~ -b2j1₊mass₊f(t)
Differential(t)(b3j0₊mass₊f(t)) ~ b3j0₊mass₊e(t)*(b3j0₊mass₊I^-1)
Differential(t)(b2j1₊mass₊f(t)) ~ b2j1₊mass₊e(t)*(b2j1₊mass₊I^-1)




sys = structural_simplify(model)
sys = structural_simplify(_model)
states(sys)
equations(sys)


sys = initialize_system_structure(alias_elimination(_model))
sys = initialize_system_structure(alias_elimination(model))
sys = tearing(sys, simplify=true)
sys = tearing(sys)
states(sys)
equations(sys)


# =============================================================================

@named power = Power()
@unpack e, f = power

ps = [mass, S[1]];
ps = [mass, s1];
# Σ efforts
eqs = [e ~ sum(p -> p.e, ps)]
# f₁ = f₂ = f₃
if length(ps) > 1
    push!(eqs, [ps[i].f ~ ps[i+1].f for i in 1:(length(ps) - 1)]...)
end
push!(eqs, ps[1].f ~ f)
s1= compose(extend(ODESystem(eqs, t, [], []; name=:s1), power), ps...)
s2 = compose(extend(ODESystem(eqs, t, [], []; name=:s2), power), ps[1])
s2 = compose(extend(ODESystem(eqs, t, [], []; name=:s2), power), ps...)

equations(s2)
states(s2)
states(s1)[1:10]

s3 = extend(ODESystem(eqs, t, [], []; name=:s3), S[1]; name=:t1)

states(s3)

# =============================================================================

# s1 = compose(extend(ODESystem(eqs, t, [], []; name=:tst), power), ps...)
# compose(extend(extend(ODESystem(eqs, t, [], []; name=:tst), S[1]), power), ps...)
# compose(extend(extend(ODESystem(eqs, t, [], []; name=:tst), S[1]), power), ps[1])


sys = initialize_system_structure(alias_elimination(s1))
sys = tearing(sys)
states(sys)


sys = structural_simplify(model)
sys = initialize_system_structure(alias_elimination(model))
sys = tearing(sys)
equations(sys)


S = [Junction0(sd, mass; name=:b3j0)];
push!(S, Junction1(S[1], mass; name=:b2j1));
push!(S, Junction0(S[2], sd; name=:b2j0));
push!(S, Junction1(S[3], mass; name=:b1j1));
push!(S, Junction0(S[4], sd; name=:b1j0));
eqs = junction0(S[end], sd, mass)

@named _model = ODESystem(eqs, t)
@named model = compose(_model, mass, sd, S...)
# @named model = compose(_model, mass, sd, S...)

sys = structural_simplify(model)
equations(sys)

sys = initialize_system_structure(alias_elimination(model))
equations(sys)

sys = tearing(sys)
equations(sys)

prob = ODEProblem(sys, [1, 0, 0, 0], (0., 20.))
prob = ODEProblem(sys, [], (0., 20.))
sol = solve(prob)

plot(sol)


# function Junction0(ps...; name, couple=true)
#     if couple
#         @named power = Power()
#         @unpack e, f = power

#         # Σ flows
#         eqs = [f ~ sum(p -> p.f, ps)]
#         # e₁ = e₂ = e₃
#         if length(ps) > 1
#             push!(eqs, [ps[i].e ~ ps[i+1].e for i in 1:(length(ps) - 1)]...)
#         end
#         push!(eqs, ps[1].e ~ e)
#         compose(extend(ODESystem(eqs, t, [], []; name=name), power), ps...)
#     else
#         # Σ flows
#         eqs = [0 ~ sum(p -> p.f, ps)]
#         # e₁ = e₂ = e₃
#         if length(ps) > 1
#             push!(eqs, [ps[i].e ~ ps[i+1].e for i in 1:(length(ps) - 1)]...)
#         end
#         compose(ODESystem(eqs, t, [], []; name=name), ps...)
#     end
# end

