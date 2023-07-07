using DiffEqFlux, DifferentialEquations, Plots, GalacticOptim
# using DifferentialEquations, Plots, GalacticOptim
using JLD2, LinearAlgebra

using Flux, Statistics
using DiffEqSensitivity
using Flux.Optimise: ADAM, update!
using Printf


# -----------------------------------------------------------------------------
# Load data experiment 1

function calcRe(V, ρ, μ, d)
    return ρ .* V .* d ./ μ
end

g = 9.81        # m/s^2 - Gravity

file = jldopen("../numerical/data/sim_doe1.jld2")

for i in keys(file)

    ex = file[i]
    ρ, μ, d, L, ϵ, a = ex["ρ"], ex["μ"], ex["d"], ex["L"], ex["ϵ"], ex["a"]

    A = π*d^2/4
    Rᵥ = 128 * μ * L / (π * d^4)
    Iᵥ = ρ * L / A

    Hin = ex["Hr"] - ex["Hf"]

    Re = calcRe(ex["V"][end], ρ, μ, d)
    if (Re < 10000) & (Re > 4000)
        println(Re,"   ", i)
    end
end

ex = file["116"]
ex = file["162"]
ex = file["163"]
ρ, μ, d, L, ϵ, a = ex["ρ"], ex["μ"], ex["d"], ex["L"], ex["ϵ"], ex["a"]
Δh = ex["Hr"] - ex["Hf"]
Δp = Δh*ρ*g

Re = calcRe(ex["V"][end], ρ, μ, d)
A = π*d^2/4
Rᵥ = 128 * μ * L / (π * d^4)
Iᵥ = ρ * L / A

# -----------------------------------------------------------------------------
# Define the NNs

m̂ = FastChain(FastDense(4, 30), FastDense(30, 1))
d̂ = FastChain(FastDense(5, 30, relu), FastDense(30, 1))
m̂!(u, p) = dot(m̂(u, p), [1])
d̂!(u, p) = dot(d̂(u, p), [1])

mθ = initial_params(m̂)
dθ = initial_params(d̂)

mp = [ρ, L, d, A]
dp = [μ, ϵ, L, d, A]

function diffeq!(du, u, p, t)
    du[1] = (/)((+)(Δp, (*)((*)(-1//1, u[1]), (d̂!)(dp, p[1]))), (m̂!)(mp, p[2]))
end

du = [0.0];
u0 = [0.0]
diffeq!(du, [0], [dθ, mθ], 1.0)

rng = 1:100
y = ex["V"][rng]
t = ex["t"][rng]
trng = 0:step(t):step(t)*(length(rng)-1)
u0 = [y[1].*A]
prob = ODEProblem(diffeq!, u0, 20.0, [dθ, mθ])
sol = solve(prob, Tsit5(), saveat = trng)
plot(sol)

function predict_neuralode(p)
    tmp_prob = remake(prob, p = p)
    Array(solve(tmp_prob, Tsit5(), saveat = trng))
end

function loss_neuralode(p)
    pred = predict_neuralode(p)
    loss = sum(abs2, y .- pred)
    return loss, pred
end

function loss_neuralodel(p)
    pred = predict_neuralode(p)
    loss = sum(abs2, y .- pred)
    return loss
end

callback = function (p, l, pred; doplot = true)
  display(l)
  # plot current prediction against data
  plt = scatter(trng, ode_data[1,:], label = "data")
  scatter!(plt, trng, pred[1,:], label = "prediction")
  if doplot
    display(plot(plt))
  end
  return false
end

result_neuralode = DiffEqFlux.sciml_train(loss_neuralode, prob.p, cb = callback)

Zygote.gradient(loss_neuralodel, prob.p)

# -----------------------------------------------------------------------------
# Define the NNs

loss(x, y) = sum(abs2, y .- x)
loss(x, y) = mean((y .- x).^2)

m̂ = FastChain(FastDense(4, 30, relu), FastDense(30, 1))
d̂ = FastChain(FastDense(5, 30, relu), FastDense(30, 1))
m̂!(u, p) = dot(m̂(u, p), [1])
d̂!(u, p) = dot(d̂(u, p), [1])
# m̂!(u, p) = abs(dot(m̂(u, p), [1]))
# d̂!(u, p) = abs(dot(d̂(u, p), [1]))

mθ = initial_params(m̂)
dθ = initial_params(d̂)

mp = [ρ, L, d, A]
dp = [μ, ϵ, L, d, A]

m̂!(mp, mθ)
d̂!(dp, dθ)

function diffeq!(du, u, p, t)
    du[1] = (/)((+)(Δp, (*)((*)(-1//1, u[1]), (d̂!)(p[1], p[2]))), (m̂!)(p[3], p[4]))
end

diffeqd!(du, u, p, t) = diffeq!(du, u, [dp, p, mp, mθ], t)
diffeqm!(du, u, p, t) = diffeq!(du, u, [dp, dθ, mp, p], t)

rng = 300:100:6000
rng = 300:200:8000
rng = 1:1:50
rng = 300:50:3000
rng = 1000:50:3000
rng = 9000:100:12890
rng = 1:500:12890
y = ex["V"][rng]
t = ex["t"][rng]
trng = 0:step(t):step(t)*(length(rng)-1)
u0 = [y[1].*A]
prob_d = ODEProblem(diffeqd!, u0, 20.0)
prob_m = ODEProblem(diffeqm!, u0, 20.0)

predict_ode_d(p) = Array(solve(prob_d, Tsit5(), u0=u0, p=p, saveat=trng)) ./ A
predict_ode_m(p) = Array(solve(prob_m, Tsit5(), u0=u0, p=p, saveat=trng)) ./ A

predict_ode_d(dθ)
predict_ode_m(mθ)

loss(predict_ode_d(dθ), y)

∇d = gradient(p -> loss(predict_ode_d(p), y), dθ)
∇m = gradient(p -> loss(predict_ode_m(p), y), mθ)

mθ = initial_params(m̂)
dθ = initial_params(d̂)

plot(t, y)
plot!(t, predict_ode_d(dθ)')

# Training
optd = ADAM(0.1);
optm = ADAM(0.1);

optd = ADAM(1);
optm = ADAM(1);

m̂!(mp, mθ)
d̂!(dp, dθ)
Iᵥ
Rᵥ
m̂!(mp, mθ)/d̂!(dp, dθ)
Iᵥ/Rᵥ

α = diff(collect(extrema(y)))[1]
it = 3000;
hist = zeros(it);
for i = 1:it
    # i = 2
    ∇d = gradient(p -> loss(predict_ode_d(p), y), dθ)
    ∇m = gradient(p -> loss(predict_ode_m(p), y), mθ)

    # if loss(predict_ode_d(dθ), y) < 42
    update!(optd, dθ, ∇d[1])
    update!(optm, mθ, ∇m[1]);
    # else
    # dθ += -(0.01*norm(dθ)/norm(∇d[1]))*∇d[1]
    # mθ += -(0.01*norm(mθ)/norm(∇m[1]))*∇m[1]
    # end

    # loss(predict_ode_d(dθ), y)
    hist[i] = loss(predict_ode_d(dθ), y)

    gloss = mean(diff(hist[max(i - 4, 1):i]))
    # @printf("It: %d - loss: %.5e - dloss %.3e\n", i, hist[i], gloss)
    @printf("It: %d - loss: %.5e m: %.5e d: %.5e\n", i, hist[i], m̂!(mp, mθ), d̂!(dp, dθ))

    if i == 1
        plot(t, y)
    elseif (i % 30) == 0
        p = plot!(t, predict_ode_d(dθ)')
        display(p)
        # plot!(t, y)
    end

    if sum(diff(hist[max(i - 6, 1):i]) .> 0) > 10
        break
    end
    # if hist[]
    # i += 1
end

GC.gc()

i = findfirst(hist .==0.0)
plot(hist[1:i-1], yaxis=:log10)
plot(t, predict_ode_d(dθ)')
plot!(t, y)


# -----------------------------------------------------------------------------
# Laurent with custom flux layer

# custom join layer
struct Join{T, F}
    combine::F
    paths::T
end

# custom split layer
struct Split{T}
  paths::T
end

Join(combine, paths...) = Join(combine, paths)
Split(paths...) = Split(paths)

Flux.@functor Join
Flux.@functor Split

(m::Join)(xs::Tuple) = m.combine(map((f, x) -> f(x), m.paths, xs)...)
(m::Join)(xs...) = m(xs)
(m::Split)(x::AbstractArray) = tuple(map(f -> f(x), m.paths))[1]

nm = 4
mm = Chain(Dense(nm, nm), Split(Dense(nm, 1, x -> x.^-2), Dense(nm, 1, x -> x.^-1), Dense(nm, 1), Dense(nm, 1, x -> x.^2)), Join(vcat, Dense(1, 1), Dense(1, 1), Dense(1, 1), Dense(1, 1)), Dense(4, 1))

nd = 5
dm = Chain(Dense(nd, nd), Split(Dense(nd, 1, x -> x.^-2), Dense(nd, 1, x -> x.^-1), Dense(nd, 1), Dense(nd, 1, x -> x.^2)), Join(vcat, Dense(1, 1), Dense(1, 1), Dense(1, 1), Dense(1, 1)), Dense(4, 1))

mθ, m̂ = Flux.destructure(mm);
m̂!(x, p) = m̂(p)(x)[1]

dθ, d̂ = Flux.destructure(dm);
d̂!(x, p) = d̂(p)(x)[1]

# loss(x, y) = mean((y .- x).^2)
loss(x, y) = sum(abs2, y .- x)

mp = [ρ, L, d, A]
dp = [μ, ϵ, L, d, A]

m̂!(mp, mθ)
d̂!(dp, dθ)

function diffeq!(du, u, p, t)
    du[1] = (/)((+)(6263.702421289803, (*)((*)(-1//1, u[1]), (d̂!)(p[1], p[2]))), (m̂!)(p[3], p[4]))
end

diffeqd!(du, u, p, t) = diffeq!(du, u, [dp, p, mp, mθ], t)
diffeqm!(du, u, p, t) = diffeq!(du, u, [dp, dθ, mp, p], t)

rng = 300:50:5500
y = ex["V"][rng]
t = ex["t"][rng]
u0 = [y[1].*A]
prob_d = ODEProblem(diffeqd!, u0, 20.0)
prob_m = ODEProblem(diffeqm!, u0, 20.0)

predict_ode_d(p) = Array(solve(prob_d, Tsit5(), u0=u0, p=p, saveat=t)) ./ A
predict_ode_m(p) = Array(solve(prob_m, Tsit5(), u0=u0, p=p, saveat=t)) ./ A

predict_ode_d(dθ)
predict_ode_m(mθ)

gradient(p -> loss(predict_ode_d(p), y), dθ)

∇d = gradient(p -> loss(predict_ode_d(p), y), dθ)
∇m = gradient(p -> loss(predict_ode_m(p), y), mθ)

mθ = randn(nvarm, length(pol))[:]
dθ = randn(nvard, length(pol))[:]

plot(t, y)
plot!(t, predict_ode_d(dθ)')

# Training
optd = ADAM(0.01);
optm = ADAM(0.01);

it = 3000;
hist = zeros(it);
# for i = 1:it
    # i = 2
    ∇d = gradient(p -> loss(predict_ode_d(p), y), dθ)
    ∇m = gradient(p -> loss(predict_ode_m(p), y), mθ)

    # if loss(predict_ode_d(dθ), y) < 42
        update!(optd, dθ, ∇d[1]);
        update!(optm, mθ, ∇m[1]);
    # else
        # dθ += -(0.01*norm(dθ)/norm(∇d[1]))*∇d[1]
        # mθ += -(0.01*norm(mθ)/norm(∇m[1]))*∇m[1]
    # end
    
    # loss(predict_ode_d(dθ), y)
    hist[i] = loss(predict_ode_d(dθ), y)

    @printf("It: %d - loss: %.5e \n", i, hist[i])

    if (i % 30) == 0
        p = plot!(t, predict_ode_d(dθ)')
        display(p)
        # plot!(t, y)
    end

    if sum(diff(hist[max(i-6, 1):i]) .> 0) > 3
        break
    end
    # if hist[]
    # i += 1
end



i = findfirst(hist .==0.0)
plot(hist[1:i-1], yaxis=:log10)
plot(t, predict_ode_d(dθ)')
plot!(t, y)


# -----------------------------------------------------------------------------
# Laurent 

# Laurent Layer
λ(X, W, n) = sum((X' * W[:, i]) .^ c for (i, c) in enumerate(n))
λ(X, W) = (X' * W[:, 1]) .^ -2 + (X' * W[:, 2]) .^ -1 + (X.* 0)'*W[:, 3] + (X' * W[:, 4]) .^ 1 + (X' * W[:, 5]) .^ 2

# Layer settings
pol = -2:2
nvarm = 4
nvard = 5
mθ = randn(nvarm, length(pol))[:]
dθ = randn(nvard, length(pol))[:]

# Create a dummy function
m̂!(u, p) = λ(u, reshape(p, nvarm, length(pol)), pol)
d̂!(u, p) = λ(u, reshape(p, nvard, length(pol)), pol)
m̂!(u, p) = λ(u, reshape(p, nvarm, 5))
d̂!(u, p) = λ(u, reshape(p, nvard, 5))

# m̂!(u, p) = Iᵥ*p[1]
# d̂!(u, p) = Rᵥ*p[1]

# loss(x, y) = sum(abs2, y .- x) + 1e5/m̂!(mp, mθ)
# loss(x, y) = sum(abs2, y .- x)
loss(x, y) = mean((y' .- x).^2)

mp = [ρ, L, d, A]
dp = [μ, ϵ, L, d, A]

m̂!(mp, mθ)
d̂!(dp, dθ)

function diffeq!(du, u, p, t)
    du[1] = (/)((+)(Δp, (*)((*)(-1//1, u[1]), (d̂!)(p[1], p[2]))), (m̂!)(p[3], p[4]))
end

diffeqd!(du, u, p, t) = diffeq!(du, u, [dp, p, mp, mθ], t)
diffeqm!(du, u, p, t) = diffeq!(du, u, [dp, dθ, mp, p], t)

rng = 1:1:50
rng = 6000:2:14000
rng = 6000:100:10000
rng = 1:200:4000
rng = 5000:100:6000
rng = 300:100:3000
rng = 9000:100:12890
rng = 1:10:1000
rng = 1:300:3000
rng = 1:500:12890
rng = 1:500:8000
stp = Int(round(length(ex["V"]) / 10))
rng = 1:stp:length(ex["V"])
y = ex["V"][rng]
t = ex["t"][rng]
trng = 0:step(t):step(t)*(length(rng)-1)
u0 = [y[1].*A]
prob_d = ODEProblem(diffeqd!, u0, 20.0)
prob_m = ODEProblem(diffeqm!, u0, 20.0)

predict_ode_d(p) = Array(solve(prob_d, Tsit5(), u0=u0, p=p, saveat=trng, reltol=1e-8, abstol=1e-8)) ./ A
predict_ode_m(p) = Array(solve(prob_m, Tsit5(), u0=u0, p=p, saveat=trng, reltol=1e-8, abstol=1e-8)) ./ A

predict_ode_d(dθ)
predict_ode_m(mθ)

∇d = gradient(p -> loss(predict_ode_d(p), y), dθ)
∇m = gradient(p -> loss(predict_ode_m(p), y), mθ)

#mθ = randn(nvarm, length(pol))[:]
#dθ = randn(nvard, length(pol))[:]
mθ = rand(nvarm, length(pol))[:]
dθ = rand(nvard, length(pol))[:]

#plot(sims[1]["t"], predict(p, sims[1], prob))
plot(t, y)
plot!(t, predict_ode_d(dθ)')

m̂!(mp, mθ)
d̂!(dp, dθ)
Iᵥ
Rᵥ
m̂!(mp, mθ)/d̂!(dp, dθ)
Iᵥ/Rᵥ

loss(predict_ode_m(mθ), y)

mean((predict_ode_m(mθ)' .- y).^2)


# Training

optd = ADAM(0.005);
optm = ADAM(0.005);

optd = ADAM(10);
optm = ADAM(10);

optd = ADAM(1);
optm = ADAM(1);

it = 20000;
hist = zeros(it);
for i = 1:it
    # i = 2
    ∇d = gradient(p -> loss(predict_ode_d(p), y), dθ)
    ∇m = gradient(p -> loss(predict_ode_m(p), y), mθ)

    # if loss(predict_ode_d(dθ), y) < 42
    update!(optm, mθ, ∇m[1])
    update!(optd, dθ, ∇d[1])
    # else
    #dθ += -(5e-5*norm(dθ)/norm(∇d[1]))*∇d[1]
    #dθ += -(1e-6*norm(dθ)/norm(∇d[1]))*∇d[1]
    #dθ += -(1e-10*norm(dθ)/norm(∇d[1]))*∇d[1]
    #mθ += -(1e-4*norm(mθ)/norm(∇m[1]))*∇m[1]
    #mθ += -(1e-6*norm(mθ)/norm(∇m[1]))*∇m[1]
    #mθ += -(1e-7*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-2*norm(mθ)/norm(∇m[1]))*∇m[1]
    # end

    # loss(predict_ode_d(dθ), y)
    # hist[i] = loss(predict_ode_d(dθ), y)
    hist[i] = loss(predict_ode_m(mθ), y)

    @printf("It: %d - loss: %.5e m: %.5e d: %.5e\n", i, hist[i], m̂!(mp, mθ), d̂!(dp, dθ))

    if i == 1
        plot(t, y)
    elseif ((i % 500) == 0)
        p = plot!(t, predict_ode_d(dθ)')
        display(p)
        # plot!(t, y)
    end

    if sum(diff(hist[max(i - 6, 1):i]) .> 0) > 50
        break
    end
    # if hist[]
    # i += 1
end

plot(t, y)
plot!(t, predict_ode_d(dθ)')
plot(t, ((predict_ode_d(dθ)' .- y).^2))

i = findfirst(hist .==0.0)
# i = length(hist)
plot(hist[1:i-1], yaxis=:log10)
# plot(hist, yaxis=:log10)
plot(t, predict_ode_d(dθ)')
plot!(t, y)