using Plots
using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra

# @parameters σ ρ β
# @variables t x(t) y(t) z(t)
# D = Differential(t)

# eqs = [D(x) ~ σ*(y-x),
#         D(y) ~ x*(ρ-z)-y,
#         D(z) ~ x*y - β*z]

# de = ODESystem(eqs,t,[x,y,z],[σ,ρ,β])

@variables t
D = Differential(t)

function Mass(; name, m = 1.0, x = 0. , u = 0.)
    ps = @parameters m=m
    sts = @variables x(t)=x v(t)=u
    eqs = D(x) ~ v
    ODESystem(eqs, t, [x, v], ps; name)
end

function Spring(; name, k = 10)
    ps = @parameters k=k
    @variables x(t)
    ODESystem(Equation[], t, [x], ps; name)
end

function Damper(; name, c = 0.1)
    ps = @parameters c=c
    @variables v(t)
    ODESystem(Equation[], t, [v], ps; name)
end

function connect_spring(spring, a, b)
    [
        # spring.x ~ norm(collect(a .- b))
        spring.x ~ a - b
    ]
end

function connect_damper(damper, a, b)
    [
        # spring.x ~ norm(collect(a .- b))
        damper.v ~ a - b
    ]
end

spring_force(spring) = -spring.k*spring.x
damper_force(damper) = -damper.c*damper.v

m = 1.0
X = 1.
k = 10.
c = 1
center = 0.
@named mass = Mass(m=m, x=X)
@named spring = Spring(k=k)
@named damper = Damper(c=c)

eqs = [
    connect_spring(spring, mass.x, center)
    connect_damper(damper, mass.v, center)
    D(mass.v) ~ spring_force(spring) / mass.m + damper_force(damper) / mass.m
]

@named _model = ODESystem(eqs, t)
@named model = compose(_model, mass, spring, damper)

sys = structural_simplify(model)
sys.eqs

prob = ODEProblem(sys, [], (0., 10.))
sol = solve(prob)

plot(sol)

Δt = 0.01
t = 0:Δt:20;
t = 0:Δt:4;
x̄ = sol(t);

plot(x̄)

# =============================================================================

@variables t
D = Differential(t)


function Mass(; name, m = 1.0, x = 0. , u = 0.)
    ps = @parameters I=m
    sts = @variables q(t)=x f(t)=u e(t)
    eqs = [
            D(q) ~ f
            D(f) ~ e/I
        ]
    ODESystem(eqs, t, [q, f, e], ps; name)
end

function Damper(; name, c = 10, u = 0.)
    ps = @parameters R=c
    @variables f(t)=u e(t)
    eqs = [
            e ~ f*R
        ]
    ODESystem(eqs, t, [e, f], ps; name)
    # ODESystem(Equation[], t, [q], ps; name)
end

function Spring(; name, k = 10, x = 0. , u = 0.)
    ps = @parameters C=1/k
    @variables q(t)=x# 2DOF

@named sd = Junction1(spring, damper)
@named bb = Junction1(mass, spring, damper)

eqs = junction0(mass, sd, bb)

@named _model = ODESystem(eqs, t)
@named model = compose(_model, spring, damper, mass, sd, bb)
@named model = compose(_model, sd, bb, mass)

sys = structural_simplify(model)
equations(sys)

sys = initialize_system_structure(alias_elimination(model))
equations(sys)

sys = tearing(sys)
equations(sys)

prob = ODEProblem(sys, [1, 0, 0, 0], (0., 20.))
sol = solve(prob)

plot(sol)
        ]
    ODESystem(eqs, t, [q, e, f], ps; name)
    # ODESystem(Equation[], t, [q], ps; name)
end

function junction1(m, s, d)
    [
        m.e ~ -d.e  - s.e
        m.f ~ d.f
        m.f ~ s.f
    ]
end

m = 1.0
X = 1.
k = 10.
c = 1
center = 0.
@named mass = Mass(m=m, x=X)
@named spring = Spring(k=k, x=X)
@named damper = Damper(c=c)

eqs = junction1(mass, spring, damper)

@named _model = ODESystem(eqs, t)
@named model = compose(_model, mass, spring, damper)

sys = structural_simplify(model)
sys.eqs

prob = ODEProblem(sys, [], (0., 20.))
sol = solve(prob)

plot(sol)

# =============================================================================

@variables t
D = Differential(t)

function Mass(; name, m = 1.0, x = 0. , u = 0.)
    ps = @parameters I=m
    sts = @variables q(t)=x f(t)=u e(t)
    eqs = [
            D(q) ~ f
            D(f) ~ e/I
        ]
    ODESystem(eqs, t, [q, f, e], ps; name)
end

function Damper(; name, c = 10, u = 0.)
    ps = @parameters R=c
    @variables f(t)=u e(t)
    eqs = [
            e ~ f*R
        ]
    ODESystem(eqs, t, [e, f], ps; name)
    # ODESystem(Equation[], t, [q], ps; name)
end

function Spring(; name, k = 10, x = 0. , u = 0.)
    ps = @parameters C=1/k
    @variables q(t)=x f(t)=u e(t)
    eqs = [
            e ~ q/C
            f ~ D(q)
        ]
    ODESystem(eqs, t, [q, e, f], ps; name)
    # ODESystem(Equation[], t, [q], ps; name)
end

function junction1(objs; name)
    @variables f(t) e(t)
    eqs = vcat([
        e ~ sum([o.e for o in objs])
        ],[
        f ~ o.f for o in objs
    ])
    compose(ODESystem(eqs, t, [e, f], []; name), objs...)
    # ODESystem(Equation[], t, [q], ps; name)
end

function junction0(m, j1)
    [
        m.e ~ j1.e
        m.f ~ j1.f
    ]
end


m = 1.0
X = 1.
k = 10.
c = 1
center = 0.
@named mass = Mass(m=m, x=X)
@named spring = Spring(k=k, x=X)
@named damper = Damper(c=c)
@named j1 = junction1([spring, damper])

eqs = junction0(mass, j1)

@named _model = ODESystem(eqs, t)
@named model = compose(_model,j1)
model.eqs

sys = structural_simplify(model)


prob = ODEProblem(sys, [], (0., 20.))
sol = solve(prob)

plot(sol)


function junction0(m, s, d)
    [
        m.e ~ -d.e  - s.e
        m.f ~ d.f
        m.f ~ s.f
    ]
end

# =============================================================================

@variables t
D = Differential(t)

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

function Junction1(ps...; name)
    @named power = Power()
    @unpack e, f = power

    # Σ efforts
    eqs = [e ~ sum(p -> p.e, ps)]
    # f₁ = f₂ = f₃
    if length(ps) > 1
        push!(eqs, [ps[i].f ~ ps[i+1].f for i in 1:(length(ps) - 1)]...)
    end
    push!(eqs, ps[1].f ~ f)
    compose(extend(ODESystem(eqs, t, [], []; name=name), power), ps...)
end

function Junction0(ps...; name)
    @named power = Power()
    @unpack e, f = power

    # Σ flows
    eqs = [f ~ sum(p -> p.f, ps)]
    # e₁ = e₂ = e₃
    if length(ps) > 1
        push!(eqs, [ps[i].e ~ ps[i+1].e for i in 1:(length(ps) - 1)]...)
    end
    push!(eqs, ps[1].e ~ e)
    compose(extend(ODESystem(eqs, t, [], []; name=name), power), ps...)
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


m = 1.0
X = 1.
k = 10.
c = 1
center = 0.
@named mass = Mass(m=m)
@named spring = Spring(k=k, x=X)
@named damper = Damper(c=c)
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

# 2DOF

@named sd = Junction1(spring, damper)
@named bb = Junction1(mass, spring, damper)

eqs = junction0(mass, sd, bb)

@named _model = ODESystem(eqs, t)
@named model = compose(_model, spring, damper, mass, sd, bb)
@named model = compose(_model, sd, bb, mass)

sys = structural_simplify(model)
equations(sys)

sys = initialize_system_structure(alias_elimination(model))
equations(sys)

sys = tearing(sys)
equations(sys)

prob = ODEProblem(sys, [1, 0, 0, 0], (0., 20.))
sol = solve(prob)

plot(sol)

# 4DOF

@named sd = Junction1(spring, damper)
@named bb = Junction1(mass, spring, damper)

eqs = junction0(mass, sd, bb)

@named _model = ODESystem(eqs, t)
@named model = compose(_model, spring, damper, mass, sd, bb)
@named model = compose(_model, sd, bb, mass)

sys = structural_simplify(model)
equations(sys)

sys = initialize_system_structure(alias_elimination(model))
equations(sys)

sys = tearing(sys)
equations(sys)

prob = ODEProblem(sys, [1, 0, 0, 0], (0., 20.))
sol = solve(prob)

plot(sol)

# =============================================================================

# DEBUG on the model

# sys = initialize_system_structure(alias_elimination(model))
# sys = tearing(sys)
# deqs = equations(sys)

# function move_diffs(eq::Equation,r)
#     # Do not modify `D(x) ~ ...`, already correct
#     # Ignore `δ(x) ~ ...` for now
#     if !(eq.lhs isa Term && ModelingToolkit.operation(eq.lhs) isa Differential) && !(eq.lhs isa Term && ModelingToolkit.operation(eq.lhs) isa Difference)
#        _eq = eq.rhs-eq.lhs
#        rhs = r(_eq)
#        if rhs === nothing
#            eq
#        else
#            lhs = _eq - rhs
#            if !(lhs isa Number) && ModelingToolkit.operation(lhs) isa Differential
#                lhs ~ -rhs
#            else
#                -lhs ~ rhs
#            end
#        end
#    else
#        eq
#    end
# end

# r1 = SymbolicUtils.@rule *(~~a, D(~~b), ~~c) => 0
# r2 = SymbolicUtils.@rule D(~~b) => 0
# remove_diffs = SymbolicUtils.Postwalk(SymbolicUtils.Chain([r1,r2]))
# deqs = [move_diffs(eq,remove_diffs) for eq in deqs]

# @named nm = ODESystem(deqs, t)
# # @named model = compose(_model, mass, j1)
# equations(nm)

