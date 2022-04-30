using DiffEqFlux, DifferentialEquations, Plots, GalacticOptim
using JLD2, LinearAlgebra

using Flux, Statistics
using DiffEqSensitivity
using Flux.Optimise: ADAM, update!
using Printf


# -----------------------------------------------------------------------------
# Support functions

function calcRe(V, ρ, μ, d)
    return ρ .* V .* d ./ μ
end

function getre(file, minre=1000, maxre=10000)

    out = Dict()
    for i in keys(file)

        ex = file[i]
        ρ, μ, d, L, ϵ, a = ex["ρ"], ex["μ"], ex["d"], ex["L"], ex["ϵ"], ex["a"]

        A = π * d^2 / 4
        # Rᵥ = 128 * μ * L / (π * d^4)
        # Iᵥ = ρ * L / A

        Hin = ex["Hr"] - ex["Hf"]

        Re = calcRe(ex["V"][end], ρ, μ, d)
        if (Re > minre) & (Re < maxre)
            # println(Re, "   ", i)
            out[i] = Re
        end
    end

    return out
end

# -----------------------------------------------------------------------------
# Laurent 

# Laurent Layer
λ(X, W, n) = sum((X' * W[:, i]) .^ c for (i, c) in enumerate(n))
# λ(X, W) = (X' * W[:, 1]) .^ -2 + (X' * W[:, 2]) .^ -1 + (X .* 0)' * W[:, 3] + (X' * W[:, 4]) .^ 1 + (X' * W[:, 5]) .^ 2

# Layer settings
pol = -2:2
# pol = -4:4
# nvarm = 4
# nvard = 5
nvarm = 3
nvard = 3

mθ = rand(nvarm, length(pol))[:]
dθ = rand(nvard, length(pol))[:]

# Create a dummy function
m̂!(u, p) = λ(u, reshape(p, nvarm, length(pol)), pol)
d̂!(u, p) = λ(u, reshape(p, nvard, length(pol)), pol)
# m̂!(u, p) = λ(u, reshape(p, nvarm, 5))
# d̂!(u, p) = λ(u, reshape(p, nvard, 5))

# -----------------------------------------------------------------------------
# m̂ = FastChain(FastDense(4, 5, relu), FastDense(5, 1))
# d̂ = FastChain(FastDense(5, 5, relu), FastDense(5, 1))
m̂ = FastChain(FastDense(3, 50, relu), FastDense(50, 1))
d̂ = FastChain(FastDense(3, 50, relu), FastDense(50, 1))
# # m̂!(u, p) = dot(m̂(u, p), [1])
# # d̂!(u, p) = dot(d̂(u, p), [1])
m̂!(u, p) = abs(dot(m̂(u, p), [1]))
d̂!(u, p) = abs(dot(d̂(u, p), [1]))

mθ = initial_params(m̂)
dθ = initial_params(d̂)

# -----------------------------------------------------------------------------

# Generate data
g = 9.81        # m/s^2 - Gravity

# file = jldopen("../numerical/data/sim_doe1.jld2")
file = jldopen("../numerical/data/sim_doe2.jld2")


function gendata(file, n=10, minre=1000, maxre=10000)
    Res = getre(file, minre, maxre)
    # Res = ["162", "163"]
    # Res = ["162"]
    sims = []
    for k in keys(Res)
    # for k in Res
        data = Dict()

        ex = file[k]
        ρ, μ, d, L, ϵ, a = ex["ρ"], ex["μ"], ex["d"], ex["L"], ex["ϵ"], ex["a"]
        A = π * d^2 / 4

        # Neural networks inputs
        data["A"] = A
        # data["mp"] = [ρ, L, d, A]
        # data["dp"] = [μ, ϵ, L, d, A]
        data["mp"] = [ρ, L, A]
        data["dp"] = [μ, L, d^2]
        data["Δp"] = (ex["Hr"] - ex["Hf"]) * ρ * g
        data["Re"] = calcRe(ex["V"][end], ρ, μ, d)

        data["Rᵥ"] = 128 * μ * L / (π * d^4)
        data["Iᵥ"] = ρ * L / A

        # Simulation settings
        tsim = ex["t"] .< 2
        stp = Int(round(length(ex["t"][tsim]) / n))
        rng = 1:stp:length(ex["t"][tsim])
        y = ex["V"][rng]
        t = ex["t"][rng]
        trng = 0:step(t):step(t)*(length(rng)-1)

        data["u0"] = [y[1] .* A]
        data["y"] = y
        data["t"] = t
        data["trng"] = trng

        push!(sims, data)
    end
    return sims
end

sims = gendata(file, 10, 1000, 2000)
sims2 = gendata(file, 10, 1000, 2000)
sims = sims2[[1]]

# Scaling
# for p in ["mp", "dp"]Δp
Δp
#     for i in 1:length(sims[1][p])
#         x = [s[p][i] for s in sims]
#         x̄ = mean(x)
#         ẋ = std(x)
#         for s in sims
#             s[p][i] = (s[p][i] - x̄)/ẋ
#         end
#     end
# end

# -----------------------------------------------------------------------------
# ODE

function diffeq!(du, u, p, t)
    du[1] = (/)((+)(p[5], (*)((*)(-1//1, u[1]), (d̂!)(p[1], p[2]))), (m̂!)(p[3], p[4]))
end

# -----------------------------------------------------------------------------
# Set problem

data = sims[1]
mp = data["mp"]
dp = data["dp"]
Δp = data["Δp"]
u0 = sims[1]["u0"]

diffeqd!(du, u, p, t) = diffeq!(du, u, [dp, p, mp, mθ, Δp], t)
diffeqm!(du, u, p, t) = diffeq!(du, u, [dp, dθ, mp, p, Δp], t)
# diffeqe!(du, u, p, t) = diffeq!(du, u, [dp, p[1], mp, p[2], Δp], t)

probm = ODEProblem(diffeqm!, u0, (0.0, 20.0), mθ)
probd = ODEProblem(diffeqd!, u0, (0.0, 20.0), dθ)
# probe = ODEProblem(diffeqe!, u0, (0.0, 20.0), [dθ, mθ])

function predict(p, data, prob)
    # global mp, dp, Δp
    tmp_prob = remake(prob, u0=data["u0"], p=p)
    # mp = data["mp"]
    # dp = data["dp"]
    # Δp = data["Δp"]
    diffeqd!(du, u, p, t) = diffeq!(du, u, [data["dp"], p, data["mp"], mθ, data["Δp"]], t)
    diffeqm!(du, u, p, t) = diffeq!(du, u, [data["dp"], dθ, data["mp"], p, data["Δp"]], t)
    vec(solve(tmp_prob, Tsit5(), saveat=data["trng"]))./data["A"]
end


# -----------------------------------------------------------------------------
# Define loss functions

function loss(p, prob)
    n = length(sims[1]["y"]) * length(sims)
    # return sum((s["y"] .- predict_neuralode(p, s)) .^ 2)/n
    return sqrt(sum([sum((s["y"] .- predict(p, s, prob)) .^ 2) for s in sims])/n)
end

# p = dθ
# prob = probd
# n = length(sims[1]["y"]) * length(sims)
# sum([sum((s["y"] .- predict(p, s, prob)) .^ 2) for s in sims]) / n
# mean(vcat([(s["y"] .- predict(p, s, prob)) .^ 2 for s in sims]...))

# -----------------------------------------------------------------------------
# Check gradient functions

∇d = gradient((p) -> loss(p, probd), dθ)
∇m = gradient((p) -> loss(p, probm), mθ)

# Generate random weights
mθ = randn(nvarm, length(pol))[:]
dθ = randn(nvard, length(pol))[:]

# mθ = initial_params(m̂)
# dθ = initial_params(d̂)

# -----------------------------------------------------------------------------
# Plot initial
ns = 1
data = sims[ns]
plot(data["t"], data["y"])
plot!(data["t"], predict(dθ, sims[ns], probd))
# plot(data["t"], predict(dθ, sims[1], probd))

loss(dθ, probd)

# -----------------------------------------------------------------------------
# Training
optd = ADAM(1);
optm = ADAM(1);

# optd = AdaBelief(1);
# optm = AdaBelief(1);

sims2 = gendata(file, 10, 1000, 1600)
sims = sims2[[1]]
sims = sims2[[1, 2]]

it = 24000;
hist = zeros(it);
for i = 1:it
    # i = 2
    ∇d = gradient((p) -> loss(p, probd), dθ)
    ∇m = gradient((p) -> loss(p, probm), mθ)

    update!(optm, mθ, ∇m[1])
    update!(optd, dθ, ∇d[1])
    # dθ += -(5e-2*norm(dθ)/norm(∇d[1]))*∇d[1]
    # dθ += -(5e-5*norm(dθ)/norm(∇d[1]))*∇d[1]
    # dθ += -(1e-6*norm(dθ)/norm(∇d[1]))*∇d[1]
    # dθ += -(1e-10*norm(dθ)/norm(∇d[1]))*∇d[1]
    # mθ += -(1e-1*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-4*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-6*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-7*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-8*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-2*norm(mθ)/norm(∇m[1]))*∇m[1]

    hist[i] = loss(dθ, probd)
    mp = sims[1]["mp"]
    dp = sims[1]["dp"]
    Iᵥ = sims[1]["Iᵥ"]
    Rᵥ = sims[1]["Rᵥ"]
    @printf("It: %d - loss: %.9e m: %.5e/%.5e d: %.5e/%.5e\n", i, hist[i], m̂!(mp, mθ), Iᵥ,d̂!(dp, dθ), Rᵥ)
    print

    if i == 1
        plot(sims[ns]["t"], sims[ns]["y"])
    elseif ((i % 100) == 0)
        p = plot!(sims[1]["t"], predict(dθ, sims[1], probd))
        # p = plot()
        # for (i, s) in enumerate(sims)
        #     c = palette(:default)[i%16]
        #     plot!(s["t"], s["y"], ls=:dash, linecolor=c)sum
        #     p = plot!(s["t"], predict(dθ, s, probd), linecolor=c)
        # end
        display(p)
        
    end

    if sum(diff(hist[max(i - 6, 1):i]) .> 0) > 100
        break
    end

end
sims

plot()
for (i, s) in enumerate(sims)
    c = palette(:default)[i%16]
    plot!(s["t"], s["y"], ls=:dash, linecolor=c)
    plot!(s["t"], predict(dθ, s, probd), linecolor=c)
end
plot!()


i = findfirst(hist .== 0.0)
# i = length(hist)
plot(hist[1:i-1], yaxis=:log10)
# plot(hist, yaxis=:log10)
plot(t, predict_ode_d(dθ)')
plot!(t, y)
