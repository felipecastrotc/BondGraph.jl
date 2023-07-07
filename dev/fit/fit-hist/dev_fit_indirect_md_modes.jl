include("lib_dev_fit.jl")

# -----------------------------------------------------------------------------
# Get data

file = jldopen("../numerical/data/sim_doe5.jld2", "r")

sims_all = gendata(file, 30, 0, 2200, remove_slow=0.9);
sims = sims_all[[1, 2]]

# -----------------------------------------------------------------------------
# Approx functions

k̂!(u, p) = p[1] * u[2];
d̂!(u, p) = p[1] * u[1];
m̂!(u, p) = p[1];

θ = [1e6, 1e6, 1e6];

# -----------------------------------------------------------------------------
# Differential equation setup

function diffeq!(du, u, p, t)
    return du[1] = (p[1] - d̂!(u[1], p[2])) / m̂!(u[1], p[3])
end

function diffeqmsd!(du, u, p, t)
    du[1] = (p[1] - k̂!(u, p[2]) - d̂!(u, p[3])) / m̂!(u, p[4])
    du[2] = u[1]
end

function diffeqmsd!(du, u, p, t)
    du .= nn(vcat(u, p[1]), p[2:end])
end

#prob = ODEProblem(diffeqmsd!, [0.0, 0.0], (0.0, 20.0), [1e6, 1e6, 1e6, 1e6])
prob = ODEProblem(diffeqmsd!, [0.0, 0.0], (0.0, 200.0), [])

# -----------------------------------------------------------------------------
# Train setup

function predict(p, nsim, prob; t=nothing)
    saveat = isnothing(t) ? nsim["trng"] : t
    tmp_prob = remake(prob; u0=vcat(nsim["u0"], 0.0), p=vcat(nsim["Δp"], p))
    out = solve(tmp_prob, Tsit5(); saveat=saveat, reltol=1e-8, abstol=1e-8)
    return vec(out[1, :]) ./ nsim["A"]
end

function loss(p, prob)
    Ŷ = [predict(p, s, prob) for s in sims]
    ρ = mean([-cor(s["y"], ŷ) for (ŷ, s) in zip(Ŷ, sims)])
    rmse = sqrt(
        mean([mean((s["y"] .- ŷ) .^ 2) for (ŷ, s) in zip(Ŷ, sims)])
    )
    return ρ + rmse
end

# function lossres(p, prob, res, t)
#     Ŷ = [predict(p, s, prob, t=t) for s in sims]
#     ρ = mean([-cor(s, ŷ) for (ŷ, s) in zip(Ŷ, res)])
#     rmse = sqrt(
#         mean([mean((s .- ŷ) .^ 2) for (ŷ, s) in zip(Ŷ, res)])
#     )
#     return ρ + rmse
# end
function lossres(p, prob, res, t)
    Ŷ = [predict(p, s, prob, t=t) for s in sims]
    ρ = mean([-cor(s, ŷ) for (ŷ, s) in zip(Ŷ, res)])
    rmse = sqrt(
        mean([mean((s .- ŷ) .^ 2) for (ŷ, s) in zip(Ŷ, res)])
    )
    return ρ + rmse + 10/p[3]
end

function lossout(p, prob)
    Ŷ = [predict(p, s, prob) for s in sims]
    ρ = mean([-cor(s["y"], ŷ) for (ŷ, s) in zip(Ŷ, sims)])
    rmse = sqrt(
        mean([mean((s["y"] .- ŷ) .^ 2) for (ŷ, s) in zip(Ŷ, sims)])
    )
    # θ = sum(log10.(abs.(dθ))) + sum(log10.(abs.(mθ)))
    θ = 0.0
    return ρ, rmse, θ
end

function lossoutres(p, prob, res, t)
    Ŷ = [predict(p, s, prob, t=t) for s in sims]
    ρ = mean([-cor(s, ŷ) for (ŷ, s) in zip(Ŷ, res)])
    # println(size(Ŷ[1]), "   ", size(res[1]))
    rmse = sqrt(
        mean([mean((s .- ŷ) .^ 2) for (ŷ, s) in zip(Ŷ, res)])
    )
    return ρ, rmse
end


# Generate gradient functions
∇θ = gradient((p) -> loss(p, prob), [1e6, 1e6, 1e6])
# ∇θ = gradient((p) -> loss(p, probs), [1e6, 1e6])

# -----------------------------------------------------------------------------
# Train

tik = 0.0
c = 0.0

info = Dict(
    :s => 0,
    :ρ => 0.0,
    :θ => 0.0,
    :rmse => 0.0,
    :tik => tik,
    :eta => 0.0,
)

sims_train = sims_all

# θ = [1e6, 1e6]
θ = [1e6, 1e6, 1e6]

s = 10
sims = sims_train[[s]]
# info[:s] = k

it = 400000
α = 1000
hist = zeros(it)

for i = 1:it
    # ∇θ = gradient((p) -> loss(p, probs), θ)
    ∇θ = gradient((p) -> loss(p, prob), θ)

    λ = 5e-3
    θ += -(λ * norm(θ) / norm(∇θ[1])) * ∇θ[1]

    ρ, rmse, pru = lossout(θ, prob)
    hist[i] = ρ + rmse
    info[:ρ], info[:rmse], info[:θ] = ρ, rmse, pru

    # Stop criteria
    if stopcrit!(i, hist, α)
        break
    elseif (ρ < -(1 - 1e-4)) & (rmse < 1e-2)
        break
    end

    println(printit(i, hist[i], info))
end

nsim = sims[1]
plot(nsim["t"], predict(θ, nsim, prob))
plot!(nsim["t"], nsim["y"])

k = nsim["key"]
ts = file[k]["t"]
ŷs = predict(θ, nsim, prob, t=ts)

# using LsqFit
# model(x, p) = @. p[1] - p[2]*exp(-p[3]*x)

n = length(ts)
n = Int(round(length(ts) / 2))
n = 300
y = file[k]["V"]

# mdl = curve_fit(model, ts[1:n], y[1:n], [1.0, 1.0, 1.0])
# ŷ = model(ts, mdl.param)

n = 100
W = hcat(ts[1:n])
a = W \ y[1:n]

plot(ts[1:n], y[1:n], legend=false)
plot!(ts[1:n], W * a)
# plot!(ts[1:n], ŷ[1:n])
plot!(ts[1:n], ŷs[1:n])


ϕ = Dict(:s => 1e4, :v => 1e3, :pa => 1e4)
ϕ[:m] = ϕ[:pa]/(ϕ[:s]/ϕ[:v])
ϕ[:k] = (ϕ[:s]^2)*ϕ[:m]
ϕ[:c] = (ϕ[:s])*ϕ[:m]

sims = deepcopy(sims_all[[s]])
sims[1]["Δp"] /= ϕ[:pa]
# nn = FastChain(FastDense(3, 10, tanh), FastDense(10, 2))
# γ = initial_params(nn)

tres = ts[1:n]*ϕ[:s]
res = [-(W * a .- y[1:n])*ϕ[:v]]
plot(tres, res[1])

γ = [10.0, 1.0, 100.0]
it = 400000
α = 300
hist = zeros(it)

opt = ADAM(0.1)

for i = 1:it
    ∇γ = gradient((p) -> lossres(p, prob, res, tres), γ)

    update!(opt, γ, ∇γ[1])

    ρ, rmse = lossoutres(γ, prob, res, tres)
    hist[i] = ρ + rmse
    info[:ρ], info[:rmse], info[:θ] = ρ, rmse, 0.0

    # Stop criteria
    if stopcrit!(i, hist, α)
        break
    elseif (ρ < -(1 - 1e-4)) & (rmse < 1e-2)
        break
    end

    println(printit(i, hist[i], info))
end

# γ = [100, 0.01, 10]
nsim = sims[1]
plot(tres, predict(γ, nsim, prob, t=tres), legend=false)
#plot(tres, predict([2e12, 2e6, 3.96e4], nsim, prob, t=tres))
# plot(tres, predict(γ, nsim, prob, t=tres))
plot!(tres, res)

γ
nsim = deepcopy(sims_all[s])
Γ = [ϕ[:k]*γ[1], ϕ[:v]*γ[2], ϕ[:m]*γ[3]]
plot(tres/ϕ[:s], predict(Γ, nsim, prob, t=tres/ϕ[:s]))
plot!(tres/10000, res/1000)

plot(predict(Γ, nsim, prob, t=tres/ϕ[:s]) .- res[1]./1000)








using Interpolations

x = 2.0.*randn(10001)
t = 0:0.01:100
f = LinearInterpolation(t, x)

#f = x -> 1.0 + x*0

function diffeq2(du, u, p, t)
    du[1] = (f(t) -p[1]*u[2] -p[2]*u[1])/p[3]
    du[2] = u[1]
end

prob = ODEProblem(diffeq2, [0.0, 0.0], [0.0, 20.0], [2, 2.0, 1.0])
sol = solve(prob)

plot(sol)
