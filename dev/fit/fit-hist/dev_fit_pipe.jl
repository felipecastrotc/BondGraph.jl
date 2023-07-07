using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D
using DifferentialEquations, JLD2

using Flux, Statistics
using DiffEqSensitivity
using Flux.Optimise: ADAM, update!
using Printf

using JLD2
using LinearAlgebra

# -----------------------------------------------------------------------------
# Functions

function MassU(h; name, P = 1.0, Θ=1.0, u = 0.0)
    e, f = Power(flow = u)
    
    ps = @parameters p = P θ = Θ

    eqs = [
        # D(f) ~ e / h(e, f, p, θ),
        D(f) ~ e / h(p, θ),
    ]
    ODESystem(eqs, t, [e, f], ps; name = name)
end

function SpringU(h ; name, x = 0.0, P = 1.0, Θ=1.0)
    e, f = Power()

    @variables q(t) = x

    ps = @parameters p = P θ = Θ

    eqs = [
        # e ~ q*h(q, e, f, p, θ)
        e ~ q*h(p, θ)
        D(q) ~ f
    ]
    ODESystem(eqs, t, [e, f, q], ps; name = name)
end

function DamperU(h; name, P = 1.0, Θ=1.0)
    e, f = Power()

    ps = @parameters p = P θ = Θ

    # eqs = [e ~ f*h(e, f, p, θ)]
    eqs = [e ~ f*h(p, θ)]

    ODESystem(eqs, t, [e, f], ps; name = name)
end

# -----------------------------------------------------------------------------
# Laminar without fluid capacitance

Hin = 0.01        # m - ΔH

g = 9.81        # m/s^2 - Gravity

l = 1.0         # m - Pipe segment length
d = 0.01        # m - Pipe segment diameter
A = π * d^2 / 4     # m^2 - Pipe section area
μ = 1e-3        # Pa*s - Viscosity
ρ = 998         # kg/m³ - Water
B = 2.2 * 1e9       # Pa -> Fluid bulk modulus

Rᵥ = 128 * μ * l / (π * d^4)
Cᵥ = A * l / B
Iᵥ = ρ * l / A

@named Pin = Se(Hin * ρ * g)

@named C = Spring(k = Cᵥ)
@named R = Damper(c = Rᵥ)
@named I = Mass(m = Iᵥ)

@named j1 = Junction1(Pin, -I, -R, -C, sgn = -1, couple=false)
@named sys = reducedobs(structural_simplify(j1))
equations(sys)

ts = (0.0, 40)
prob = ODEProblem(sys, [], ts)
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)

plot(sol.t, sol[I.f]./ A)

# -----------------------------------------------------------------------------
# Load data experiment 1

file = jldopen("../numerical/data/sim_doe1.jld2")

ex = file["1"]
ρ, μ, d, L, ϵ, a = ex["ρ"], ex["μ"], ex["d"], ex["L"], ex["ϵ"], ex["a"]

A = π*d^2/4
Rᵥ = 128 * μ * L / (π * d^4)
Iᵥ = ρ * L / A

Hin = ex["Hr"] - ex["Hf"]

g = 9.81        # m/s^2 - Gravity
# -----------------------------------------------------------------------------
# Only Mass
# -----------------------------------------------------------------------------
# Linear elements as reference 
@named Pin = Se(Hin * ρ * g)
@named R = Damper(c = Rᵥ)
@named M = Mass(m = 2*Iᵥ)

# Build system
@named j1 = Junction1(Pin, -M, -R, sgn = -1, couple=false)
@named sys = reducedobs(structural_simplify(j1))
equations(sys)

# Solve
prob = ODEProblem(sys, [], (0.0, 20.0))
sol = solve(prob)
plot(sol)

# -----------------------------------------------------------------------------
# Simple function to debug
# M => ρ, L, d, A
# NN
mnn = Chain(Dense(4, 30), Dense(30, 1));
mnn = Flux.f64(mnn);
mθ, mre = Flux.destructure(mnn);
hm(p, θ) = dot(mre(θ)([p[1], p[2], p[3], p[4]]), [1])
# Linear
hm(p, θ) = p*θ
@register hm(p, θ)

# Set elements
@named Mu = MassU(hm)

# Build system
@named j1 = Junction1(Pin, -Mu, -R, sgn = -1, couple=false)
@named sysu = reducedobs(structural_simplify(j1))
equations(sysu)

# Solve
probu = ODEProblem(sysu, [], (0.0, 20.0))
solu = solve(probu)
plot!(solu)

parameters(sysu)
u0 = prob.u0;

# Linear
ps = [2*Iᵥ, Rᵥ, 1]
# Dot product
ps = [Rᵥ, [2, 0,0 ], [Iᵥ, 1, 2]']
# NN
pm = [ρ, L, d, A]
ps = [Rᵥ, mθ, pm]

# Using a different approach to solve the system
t_rng = 0:0.1:20
out = Array(solve(probu, Tsit5(), u0 = u0, p = ps, saveat = t_rng))

# Working
plot!(t_rng, out')

# -----------------------------------------------------------------------------
# Linear without state variables
# -----------------------------------------------------------------------------
# Mass + resistance
# R => μ, ϵ, L, d, A
# M => ρ, L, d, A

# NN mass
mnn = Chain(Dense(4, 30), Dense(30, 1));
mnn = Flux.f64(mnn);
mθ, mre = Flux.destructure(mnn);
hm(p, θ) = dot(mre(θ)([p[1], p[2], p[3], p[4]]), [1])
@register hm(p, θ)

# NN damper
dnn = Chain(Dense(5, 30), Dense(30, 1))
dnn = Flux.f64(dnn)
dθ, dre = Flux.destructure(dnn);
hd(p, θ) = dot(dre(θ)([p[1], p[2], p[3], p[4], p[5]]), [1])
@register hd(p, θ)

# Set elements
@named Pin = Se(Hin * ρ * g)
@named Mu = MassU(hm)
@named Ru = DamperU(hd)

# Build system
@named j1 = Junction1(Pin, -Mu, -Ru, sgn = -1, couple=false)
@named sysu = reducedobs(structural_simplify(j1))
equations(sysu)

# Parameters system
parameters(sysu)
pm = [ρ, L, d, A]
pd = [μ, ϵ, L, d, A]
ps = [dθ, mθ]
ps = mθ
ps = [pd, dθ, mθ, pm]

# Convert the system to a function
f = ODEFunction(sysu);
parameters(sysu)
# dummy(du, u, p, t) = f(du, u, [pd, p[1], p[2], pm], t);
dummy(du, u, p, t) = f(du, u, [pd, dθ, p, pm], t);
dummy(du, u, p, t) = f(du, u, p, t);

# Test dummy function
du = [0.0]
u0 = [0.0]
dummy(du, u0, ps, 0.0)

# Create a ODE problem from the system function
probu = ODEProblem(dummy, [0.0], (0.0, 3.0))

# Initial condition
states(sysu)
u0 = [0.0];
y = ex["V"]
t_rng = ex["t"]
out = Array(solve(probu, Tsit5(), u0 = u0, p = ps, saveat=t_rng))[1, :]
plot(t_rng, out)


# Now trying to fit
function pred_ode(p)
    ps = [pd, p[1], p[2], pm]
    Array(solve(probu, Tsit5(), u0 = u0, p = ps, saveat = t_rng))[1, :]
end

loss(x, y) = mean((x .- y) .^ 2);
loss(pred_ode([dθ, mθ]), y)

y = ex["V"][1:100]
t_rng = ex["t"][1:100]


gradient((p) -> loss(pred_ode(p), y), dθ, mθ)


g = gradient((p) -> loss(pred_ode(p), y), ps)

grads = gradient(() -> loss(pred_ode(dθ, mθ), y), dθ, mθ)



loss(y) = mean((pred_ode(ps) .- y) .^ 2);
loss(y)


gradient(() -> loss(y), params([dθ, mθ]))



s




















# opt = AMSGrad();
opt = ADAM(0.01);

loss(x, y) = mean((x .- y) .^ 2);


pred_ode(θ)
g = gradient((θ) -> loss(pred_ode(θ), y), θ)

gradient(θ -> loss(pred_ode(θ), y), θ)

# Iterate
it = 10;
hist = zeros(it);
for i = 1:it
    g = gradient((p) -> loss(pred_ode(θ), y), θ)
    # ∇loss = gradient(() -> loss(pred_ode(p), y), ps);
    # collect(∇loss)
    # update!(opt, ps, ∇loss);
    # update!(opt, ps, ∇loss);
    p += -g[1]
    hist[i] = loss(pred_ode(p), y)
    @printf("It: %d - loss: %.3e \n", i, hist[i])
end



# -----------------------------------------------------------------------------
# Non-linear with state variables
# -----------------------------------------------------------------------------
# Mass + resistance

# R => μ, ϵ, L, d, A
# M => ρ, L, d, A
# K => ρ, L, A, a 


# NN mass
mnn = Chain(Dense(6, 30), Dense(30, 1));
mnn = Flux.f64(mnn);
mθ, mre = Flux.destructure(mnn);
hm(e, f, p, θ) = dot(re(θ)([e, f, p[1], p[2], p[3], p[4]]), [1])
@register hm(e, f, p, θ)

# NN spring
sm = Chain(Dense(7, 30), Dense(30, 1))
sm = Flux.f64(sm)
sθ, sre = Flux.destructure(sm);
hs(q, e, f, p, θ) = dot(re(θ)([q, e, f, p[1], p[2], p[3], p[4]]), [1])
@register hs(q, e, f, p, θ)

# NN damper
dm = Chain(Dense(7, 30), Dense(30, 1))
dm = Flux.f64(dm)
dθ, dre = Flux.destructure(dm);
hd(e, f, p, θ) = dot(re(θ)([e, f, p[1], p[2], p[3], p[4], p[5]]), [1])
@register hd(e, f, p, θ)