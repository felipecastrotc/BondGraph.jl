using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Symbolics, Symbolics.Latexify

include("lib_bg.jl")


@variables t
D = Differential(t)


@connector function Power(; name, effort = 0.0, flow = 0.0)
    sts = @variables e(t) = effort f(t) = flow
    ODESystem(Equation[], t, sts, []; name = name)
end

function Mass(; name, m = 1.0, u = 0.0)
    @named power = Power(flow = u)
    @unpack e, f = power
    ps = @parameters I = m
    eqs = [D(f) ~ e / I]
    extend(ODESystem(eqs, t, [], ps; name = name), power)
end

function Spring(; name, k = 10, x = 0.0)
    @named power = Power()
    @unpack e, f = power

    @variables q(t) = x
    ps = @parameters C = 1 / k

    eqs = [
        e ~ q / C
        D(q) ~ f
        # D(e) ~ f/C
    ]
    extend(ODESystem(eqs, t, [q], ps; name = name), power)
end

function Damper(; name, c = 10)
    @named power = Power()
    @unpack e, f = power

    ps = @parameters R = c
    eqs = [e ~ f * R]
    extend(ODESystem(eqs, t, [], ps; name = name), power)
end

function Junction1(ps...; name, subsys = [], couple = true)
    # Check subsys type
    sys = subsys isa ODESystem ? [subsys] : subsys

    # Get connections
    con = collect(Base.Flatten([ps, sys]))
    if couple
        @named power = Power()
        @unpack e, f = power

        # Σ efforts
        eqs = [e ~ sum(c -> c.e, con)]
        # f₁ = f₂ = f₃
        if length(con) > 1
            for i = 1:(length(con)-1)
                push!(eqs, con[i].f ~ con[i+1].f)
            end
        end
        push!(eqs, con[1].f ~ f)
        # Build subsystem
        compose(extend(ODESystem(eqs, t, [], []; name = name), power), ps...)
    else
        # Σ efforts
        eqs = [0 ~ sum(p -> p.e, con)]
        # f₁ = f₂ = f₃
        if length(con) > 1
            for i = 1:(length(con)-1)
                push!(eqs, con[i].f ~ con[i+1].f)
            end
        end
        # Build system
        compose(ODESystem(eqs, t, [], []; name = name), ps...)
    end
end

function Junction0(ps...; name, subsys = [], couple = true)
    # Check subsys type
    sys = subsys isa ODESystem ? [subsys] : subsys

    # Get connections
    con = collect(Base.Flatten([ps, sys]))
    if couple
        @named power = Power()
        @unpack e, f = power

        # Σ flows
        eqs = [f ~ sum(c -> c.f, con)]
        # e₁ = e₂ = e₃
        if length(con) > 1
            for i = 1:(length(con)-1)
                push!(eqs, con[i].e ~ con[i+1].e)
            end
        end
        push!(eqs, con[1].e ~ e)
        # Build subsystem
        compose(extend(ODESystem(eqs, t, [], []; name = name), power), ps...)
    else
        # Σ flows
        eqs = [0 ~ sum(p -> p.f, con)]
        # e₁ = e₂ = e₃
        if length(con) > 1
            for i = 1:(length(con)-1)
                push!(eqs, con[i].e ~ con[i+1].e)
            end
        end
        # Build system
        compose(ODESystem(eqs, t, [], []; name = name), ps...)
    end
end


# =============================================================================
# Definitions

m = 1.0
X = 1.0
k = 10.0
c = 1.0
@named mass = Mass(m = m)
@named spring = Spring(k = k, x = X)
@named damper = Damper(c = c)

# =============================================================================
# 1DOF

@named _model = Junction1(spring, mass, damper, couple = false)

sys = structural_simplify(_model)
equations(sys)

latexify(equations(sys)[1])

latexify(equations(sys)[1], convert_unicode = false)

prob = ODEProblem(sys, [], (0.0, 20.0))
sol = solve(prob)

plot(sol)

# =============================================================================
# 2DOF

@named sd = Junction1(spring, damper)

M = [Junction1(spring, damper, mass; name = :b1)];
push!(M, Junction0(sd; name = :b2j0, subsys = M[1]));
push!(M, Junction1(mass; name = :b2j1, subsys = M[2], couple = false));

model = ODESystem(equations(M), t; name = :tst)
sys = structural_simplify(model)
equations(sys)

prob = ODEProblem(sys, [], (0.0, 20.0))
sol = solve(prob)
plot(sol)

# Alternative scheme - REVISAR PQ ESTÀ DANDO ERRADO!!!
M2 = [Junction1(spring, damper, mass; name = :b1)];
push!(M2, Junction1(mass; name = :b2j1));
push!(M2, Junction0(sd, M2[1], M2[2]; name = :b2j0, couple = false));

model2 = ODESystem(equations(M2[3]), t; name = :tst)
sys2 = structural_simplify(model2)
equations(sys2)

equations(sys)

prob2 = ODEProblem(sys2, [], (0.0, 20.0))
prob2
sol2 = solve(prob2)
plot(sol2)

sys2.states
sys.states

Array(sol2)[1, :] - Array(sol)[3, :]
Array(sol2)[2, :] - Array(sol)[1, :]
Array(sol2)[3, :] - Array(sol)[2, :]
Array(sol2)[4, :] + Array(sol)[4, :]

# =============================================================================
# 4DOF

@named sd = Junction1(spring, damper)

M = [Junction0(sd, mass; name = :b3j0)];
push!(M, Junction1(mass; name = :b2j1, subsys = M[1]));
push!(M, Junction0(sd; name = :b2j0, subsys = M[2]));
push!(M, Junction1(mass; name = :b1j1, subsys = M[3]));
push!(M, Junction0(sd; name = :b1j0, subsys = M[4]));
push!(M, Junction1(mass, spring, damper; name = :b0, subsys = M[5], couple = false));

model = ODESystem(equations(M), t; name = :tst)
sys = structural_simplify(model)
equations(sys)

prob = ODEProblem(sys, [], (0.0, 20.0))
sol = solve(prob)

plot(sol)

# =============================================================================
# 4DOF

@named sd = Junction1(spring, damper)

M = [Junction0(sd, mass; name = :b3j0)];
push!(M, Junction1(mass; name = :b2j1, subsys = M[1]));
push!(M, Junction0(sd; name = :b2j0, subsys = M[2]));
push!(M, Junction1(mass; name = :b1j1, subsys = M[3]));
push!(M, Junction0(sd; name = :b1j0, subsys = M[4]));
push!(M, Junction1(mass, spring, damper; name = :b0, subsys = M[5], couple = false));

model = ODESystem(equations(M), t; name = :tst)
sys = structural_simplify(model)
equations(sys)

prob = ODEProblem(sys, [], (0.0, 20.0))
sol = solve(prob)

plot(sol)

# =============================================================================
# Test 1DOF

using Flux, Statistics
using DiffEqSensitivity
using Flux.Optimise: ADAM, update!
using Printf

h(x, y) = x^2 + y
h(x, p) = x * (p' * p)
h(x, p) = x * p
p = 5

m = Chain(Dense(1, 30), Dense(30, 1))
m = Flux.f64(m)
p, re = Flux.destructure(m);
h(x, p) = dot(re(p)([x]), [1])
@register h(x, p)


function SpringU(; name, x = 0.0, P = 1.0)
    @named power = Power()
    @unpack e, f = power

    @variables q(t) = x
    # @variables s2(q)

    ps = @parameters p = P

    eqs = [
        e ~ h(q, p)
        # e ~ s2
        D(q) ~ f
        # D(e) ~ f/C
    ]
    extend(ODESystem(eqs, t, [q], ps; name = name), power)
end

@named springU = SpringU(x = 1.0)
@named modelU = Junction1(mass, springU, damper, couple = false)

sysU = structural_simplify(modelU)
equations(sysU)
parameters(sysU)

probU = ODEProblem(sysU, [], (0.0, 20.0), [springU.p => p])
solU = solve(probU)
plot(solU)

# Reference for further reading
# https://github.com/SciML/ModelingToolkit.jl/issues/174
# ode_func = ODEFunction(sysU, dvs=states(sysU), parameters(sysU))

function pred_ode(p)
    Array(solve(probU, Tsit5(), u0 = u0, p = [1.0, p, 1], saveat = t_rng))
    # Array(solve(probU, saveat=t_rng))
end

plot(sol)
plot!(solU)
plot(pred_ode(p)')

t_rng = 0:0.1:20
y = Array(sol(t_rng))

ps = Flux.params(p);

# Initialise optimiser
# opt = AMSGrad();
opt = ADAM(0.01);

loss(x, y) = mean((circshift(x, 1) .- y) .^ 2);

# Iterate
# it = 19;
it = 300;
hist = zeros(it);
for i = 1:it
    g = gradient((p) -> loss(pred_ode(p), y), p)
    # ∇loss = gradient(() -> loss(pred_ode(p), y), ps);
    # collect(∇loss)
    # update!(opt, ps, ∇loss);
    # update!(opt, ps, ∇loss);
    p += -g[1]
    hist[i] = loss(pred_ode(p), y)
    @printf("It: %d - loss: %.3e \n", i, hist[i])
end

plot(pred_ode(p)')
plot!(y')

# =============================================================================
# Trying 1DOF manually

function pODE(du, u, p, t)
    m = 1.0
    c = 1.0
    du[1] = (m^-1) * (-h(u[2], p) - (c * u[1]))
    du[2] = u[1]
end

m = Chain(Dense(1, 30, relu), Dense(30, 1))
# m = Flux.f64(m)
p, re = Flux.destructure(m);
h(x, p) = re(p)([x])[1]
# h(x, p) = x*p[1]
# p = [2.0]
@register h(x, p)

t_f = 5.0
u0 = [0.0, 1.0];
tspan = (0.0, 20.0);
tspan = (0.0, t_f);
t_rng = 0:0.1:t_f

prob = ODEProblem(pODE, u0, tspan)

function predict_n_ode()
    Array(solve(prob, Tsit5(), u0 = u0, p = p, saveat = t_rng))
end

plot(predict_n_ode()')

y = Array(sol(t_rng));
y = circshift(y, 1);

ps = Flux.params(p);

# Initialise optimiser
# opt = AMSGrad();
opt = ADAM(0.05);

# loss(x, y) = mean((x .- y).^2);
loss(x, y) = sum(abs2, y .- x);

# Iterate
# it = 19;
it = 1000;
hist = zeros(it);
for i = 1:it
    ∇loss = gradient(() -> loss(predict_n_ode(), y), ps)
    update!(opt, ps, ∇loss)
    hist[i] = loss(predict_n_ode(), y)
    cb()
    @printf("It: %d - loss: %.3e \n", i, hist[i])
end


cb = function (; doplot = false) # callback function to observe training
    pred = predict_n_ode()
    #   display(loss(pred, y))
    # plot current prediction against data
    pl = scatter(t_rng, y', label = "data")
    scatter!(pl, t_rng, pred', label = "prediction")
    display(plot(pl))
end
cb()


# =============================================================================
# Test 4DOF

using Flux, Statistics
using DiffEqSensitivity
using Flux.Optimise: ADAM, update!
using Printf

function SpringU(; name, x = 0.0, P = 1.0)
    @named power = Power()
    @unpack e, f = power

    @variables q(t) = x
    # @variables s2(q)

    ps = @parameters p = P

    eqs = [
        e ~ h(q, p)
        # e ~ s2
        D(q) ~ f
        # D(e) ~ f/C
    ]
    extend(ODESystem(eqs, t, [q], ps; name = name), power)
end


h(x, p) = x * p
p = [1.0, 1.0]
@register h(x, p)

@named springU = SpringU(x = 1.0)
@named sd = Junction1(springU, damper)

MU = [Junction0(sd, mass; name = :b3j0)];
push!(MU, Junction1(mass; name = :b2j1, subsys = MU[1]));
push!(MU, Junction0(sd; name = :b2j0, subsys = MU[2]));
push!(MU, Junction1(mass; name = :b1j1, subsys = MU[3]));
push!(MU, Junction0(sd; name = :b1j0, subsys = MU[4]));
push!(MU, Junction1(mass, springU, damper; name = :b0, subsys = MU[5], couple = false));

modelU = ODESystem(equations(MU), t; name = :tst)
sysU = structural_simplify(modelU)
equations(sysU)
parameters(sysU)

probU = ODEProblem(sysU, [], (0.0, 20.0))
solU = solve(probU)
plot(solU)

# u0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
u0 = probU.u0

function pred_ode(p)
    # ps = [p[1], 1.0, 1.0, 1.0, p[2], 1.0, 1.0, p[3], 1.0, 1.0, p[4], 1.0]
    ps = [p, 1.0, 1.0, 1.0, p, 1.0, 1.0, p, 1.0, 1.0, p, 1.0]
    Array(solve(probU, Tsit5(), u0 = u0, p = ps, saveat = t_rng))
    # Array(solve(probU, saveat=t_rng))
end

p = [1, 5, 11, 9]
p = [10, 10, 10, 10]

t_rng = 0:0.1:20
y = Array(sol(t_rng))


plot(sol)
plot!(solU)
plot(pred_ode(p)')

ps = Flux.params(p);

# Initialise optimiser
# opt = AMSGrad();
opt = ADAM(0.01);

loss(x, y) = mean((x .- y) .^ 2);

p = 30
# Iterate
# it = 19;
it = 3000;
hist = zeros(it);
for i = 1:it
    g = gradient((p) -> loss(pred_ode(p), y), p)
    # ∇loss = gradient(() -> loss(pred_ode(p), y), ps);
    # collect(∇loss)
    # update!(opt, ps, ∇loss);
    # update!(opt, ps, ∇loss);
    p += -g[1]
    hist[i] = loss(pred_ode(p), y)
    @printf("It: %d - loss: %.3e \n", i, hist[i])
end

plot(pred_ode(p)')
plot!(y')
