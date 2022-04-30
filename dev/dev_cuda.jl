using DiffEqFlux, DifferentialEquations, Plots, GalacticOptim
using JLD2, LinearAlgebra

using Flux, Statistics
using DiffEqSensitivity
using Flux.Optimise: ADAM, update!
using Printf
using CUDA


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
#λ(X, W, n) = sum((X' * W[:, i]) .^ c for (i, c) in enumerate(n))
λ(X, W) = (X' * W[:, 1]) .^ -2 + (X' * W[:, 2]) .^ -1 + (X .* 0)' * W[:, 3] + (X' * W[:, 4]) .^ 1 + (X' * W[:, 5]) .^ 2

# Layer settings
pol = -2:2
nvarm = 4
nvard = 5

mθ = rand(nvarm, length(pol))[:] |> gpu
dθ = rand(nvard, length(pol))[:] |> gpu

# Create a dummy function
#m̂!(u, p) = λ(u, reshape(p, nvarm, length(pol)), pol)
#d̂!(u, p) = λ(u, reshape(p, nvard, length(pol)), pol)
m̂!(u, p) = λ(u, reshape(p, nvarm, 5))
d̂!(u, p) = λ(u, reshape(p, nvard, 5))

# -----------------------------------------------------------------------------

# Generate data
g = 9.81        # m/s^2 - Gravity

file = jldopen("../numerical/data/sim_doe1.jld2")

function gendata(file, n=10, minre=1000, maxre=10000)
    Res = getre(file, minre, maxre)

    sims = []
    for k in keys(Res)
        data = Dict()

        ex = file[k]
        ρ, μ, d, L, ϵ, a = ex["ρ"], ex["μ"], ex["d"], ex["L"], ex["ϵ"], ex["a"]
        A = π * d^2 / 4

        # Neural networks inputs
        data["A"] = A |> gpu |> gpu
        data["mp"] = [ρ, L, d, A] |> gpu
        data["dp"] = [μ, ϵ, L, d, A] |> gpu
        data["Δp"] = (ex["Hr"] - ex["Hf"]) * ρ * g |> gpu

        # Simulation settings
        stp = Int(round(length(ex["V"]) / n))
        rng = 1:stp:length(ex["V"])
        y = ex["V"][rng]
        t = ex["t"][rng]
        trng = 0:step(t):step(t)*(length(rng)-1) |> gpu

        data["u0"] = [y[1] .* A] |> gpu
        data["y"] = y |> gpu
        data["t"] = collect(t) |> gpu
        data["trng"] = collect(trng) |> gpu

        push!(sims, data)
    end
    return sims
end

sims = gendata(file, 10, 1000, 10000);

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
    global mp, dp, Δp
    tmp_prob = remake(prob, u0=data["u0"], p=p)
    mp = data["mp"]
    dp = data["dp"]
    Δp = data["Δp"]
    vec(solve(tmp_prob, Tsit5(), saveat=data["trng"]))./data["A"]
end

# -----------------------------------------------------------------------------
# Define loss functions

function loss(p, prob)
    s = sims[1]
    n = length(s["y"]) * length(sims)
    # return sum((s["y"] .- predict_neuralode(p, s)) .^ 2)/n
    return sum([sum((s["y"] .- predict(p, s, prob)) .^ 2) for s in sims])/n
end

# -----------------------------------------------------------------------------
# Check gradient functions

∇d = gradient((p) -> loss(p, probd), dθ)
∇m = gradient((p) -> loss(p, probm), mθ)

# Generate random weights
mθ = rand(nvarm, length(pol))[:] |> gpu
dθ = rand(nvard, length(pol))[:] |> gpu

# -----------------------------------------------------------------------------
# Plot initial
plot(data["t"], data["y"])
plot!(data["t"], predict(dθ, sims[1], probd))

loss(dθ, probd)

# -----------------------------------------------------------------------------
# Training
optd = ADAM(1);
optm = ADAM(1);

it = 20000;
hist = zeros(it);
for i = 1:it
    # i = 2
    ∇d = gradient((p) -> loss(p, probd), dθ)
    ∇m = gradient((p) -> loss(p, probm), mθ)

    update!(optm, mθ, ∇m[1])
    update!(optd, dθ, ∇d[1])

    hist[i] = loss(dθ, probd)

    Rᵥ = 128 * dp[1] * dp[3] / (π * dp[4]^4)
    Iᵥ = mp[1] * mp[2] / mp[4]
    @printf("It: %d - loss: %.5e m: %.3e/%.3e d: %.3e/%.3e\n", i, hist[i], m̂!(mp, mθ), Iᵥ,d̂!(dp, dθ), Rᵥ)

    if i == 1
        plot(data["t"], data["y"])
    elseif ((i % 500) == 0)
        p = plot!(data["t"], predict(dθ, sims[1], probd))
        display(p)
        # plot!(t, y)
    end

    if sum(diff(hist[max(i - 6, 1):i]) .> 0) > 100
        break
    end

end

plot(t, y)
plot!(t, predict_ode_d(dθ)')
plot(t, ((predict_ode_d(dθ)' .- y) .^ 2))

i = findfirst(hist .== 0.0)
# i = length(hist)
plot(hist[1:i-1], yaxis=:log10)
# plot(hist, yaxis=:log10)
plot(t, predict_ode_d(dθ)')
plot!(t, y)


