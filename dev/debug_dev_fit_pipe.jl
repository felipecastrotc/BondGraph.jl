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
        e ~ q*h(q, e, f, p, θ)
        D(q) ~ f
    ]
    ODESystem(eqs, t, [e, f, q], ps; name = name)
end

function DamperU(h; name, P = [1.0], Θ=[1.0])
    e, f = Power()

    # ps = @parameters p = P θ = Θ
    ps = @parameters p θ

    # eqs = [e ~ f*h(e, f, p, θ)]
    eqs = [e ~ f*h(p, θ)]

    ODESystem(eqs, t, [e, f], ps; name = name)
end

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
# hd(p, θ) = dre(θ)([p[1], p[2], p[3], p[4], p[5]])
# hd(p, θ) = dre(θ)([p[1], p[2], p[3], p[4], p[5]])[1]
@register hd(p, θ)

# Set elements
@named Pin = Se(Hin * ρ * g)
@named Mu = MassU(hm)
@named Ru = DamperU(hd)

@named M = Mass(m = 2*Iᵥ)

# Build system
# @named j1 = Junction1(Pin, -Mu, -Ru, sgn = -1, couple=false)
@named j1 = Junction1(Pin, -M, -Ru, sgn = -1, couple=false)
@named sysu = reducedobs(structural_simplify(j1))
equations(sysu)

# Solve
# probu = ODEProblem(sysu, [], (0.0, 10.0), [Ru.p => dp, Ru.θ => dθ])
probu = ODEProblem(sysu, [], (0.0, 10.0), [Ru.p => 1, Ru.θ => 2])

# Parameters system
parameters(sysu)
mp = [ρ, L, d, A]
dp = [μ, ϵ, L, d, A]
# ps = [dp, dθ, mθ, mp]
ps = [Iᵥ, dp, dθ]

# Initial condition
u0 = probu.u0;
# Using a different approach to solve the system
t_rng = 0:0.1:20
out = Array(solve(probu, Tsit5(), u0 = u0, p = ps, saveat = t_rng))
plot(out')

# Working
# plot(t_rng, out[1, :])

# Now trying to fit

function pred_ode(θ::Vector{Float64})
    # ps = [dp, θ[1:length(dθ)], θ[length(dθ)+1:end], mp]
    # ps = [dp, θ[1:211], θ[212:end], mp]
    # ps = [Iᵥ, dp, dθ]
    # ps = [Iᵥ, dp, θ]
    # ps = [Iᵥ, dp, θ]
    Array(solve(probu, Tsit5(), u0 = u0, p =  [Iᵥ, dp, θ], saveat = t_rng))[1, :]
    # solve(probu, Tsit5(), u0 = u0, p = ps, saveat = t_rng).u
end

# Flat weights
# θ = vcat(dθ, mθ)
θ = dθ

# Θ = Flux.params(θ);

# Initialise optimiser
# opt = AMSGrad();
opt = ADAM(0.01);

# loss(x, y) = mean((x .- y) .^ 2);
loss(x, y) = mean((x - y) .^ 2);
loss(x, y) = sum((x - y));

y = ex["V"]
t_rng = ex["t"]

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

# ----------------------------------------------------------------------------
# Testing how to generate and visualize the function generated

using MacroTools

myode_oop = generate_function(sysu)[2]; # first one is the out-of-place function
MacroTools.striplines(myode_oop) # print without line numbers

nlsys_func = generate_function(sysu)[2]
f = eval(nlsys_func)

# f = ODEFunction(sysu)     # Generate using mtk library

dp = [μ, ϵ, L, d, A]
ps = [Iᵥ, dp, dθ]
du = [0.0]
u0 = [0.0]
f(du, u0, ps, 0.0)

probu = ODEProblem(f, [0.0], (0.0, 20.0))

function predict_ode(dθ)
    ps = [Iᵥ, dp, dθ]
    Array(solve(probu, Tsit5(), u0 = u0, p = ps, saveat = t_rng))[1, :]
end

y = ex["V"]
t_rng = ex["t"]

predict_ode(dθ)

loss(x, y) = sum((x - y));

g = gradient((p) -> loss(predict_ode(p), y), θ)

# ----------------------------------------------------------------------------
# Manual - Working

function myf(du, u, p, t)
    du[1] = (/)((+)(6263.702421289803, (*)((*)(-1//1, u[1]), (hd)(dp, p))), Iᵥ)
end;

function predict_ode(p)
    Array(solve(probu, Tsit5(), u0 = u0, p = p, saveat = t_rng))[1, :]
end

du = [0.0];
u0 = [0.0]
myf(du, [0], dθ, 1.0)
y = ex["V"]
t_rng = ex["t"]


probu = ODEProblem(myf, [0.0], (0.0, 20.0))

predict_ode(dθ)

loss(x, y) = sum((x - y));
loss(predict_ode(dθ), y)

g = gradient((p) -> loss(predict_ode(p), y), dθ)

# ----------------------------------------------------------------------------
# Manual - extending

myode_oop = generate_function(sysu)[2]; # first one is the out-of-place function
MacroTools.striplines(myode_oop) # print without line numbers

# f = ODEFunction(sysu)
nlsys_func = generate_function(sysu)[2]
f = eval(nlsys_func)

function myf(du, u, p, t)
    f(du, u, [Iᵥ, dp, p], t)
end

function predict_ode(p)
    Array(solve(probu, Tsit5(), u0 = u0, p = p, saveat = t_rng))[1, :]
end

du = [0.0];
u0 = [0.0]
myf(du, [0], dθ, 1.0)

y = ex["V"][1:1200]
t_rng = ex["t"]

probu = ODEProblem(myf, [0.0], (0.0, 20.0))
predict_ode(dθ)

loss(predict_ode(dθ), y)

g = gradient((p) -> loss(predict_ode(p), y), dθ)

# ----------------------------------------------------------------------------
# Manual - ODEFunction

f = ODEFunction(sysu)

function myf(du, u, p, t)
    f(du, u, [Iᵥ, dp, p], t)
end

function predict_ode(p)
    Array(solve(probu, Tsit5(), u0 = u0, p = p, saveat = t_rng))[1, :]
end

du = [0.0];
u0 = [0.0]
myf(du, [0], dθ, 1.0)

y = ex["V"][1:100]
t_rng = ex["t"][1:100]

probu = ODEProblem(myf, [0.0], (0.0, 20.0))
predict_ode(dθ)

loss(predict_ode(dθ), y)

g = gradient((p) -> loss(predict_ode(p), y), dθ)

# ----------------------------------------------------------------------------
# Auto Failing

probu = ODEProblem(sysu, [], (0.0, 10.0), [Ru.p => 1, Ru.θ => 2])

function predict_ode(p)
    Array(solve(probu, Tsit5(), u0 = u0, p = [Iᵥ, dp, p], saveat = t_rng))[1, :]
end

du = [0.0];
u0 = [0.0]

y = ex["V"][1:100]
t_rng = ex["t"][1:100]

predict_ode(dθ)

loss(predict_ode(dθ), y)

g = gradient((p) -> loss(predict_ode(p), y), dθ)


