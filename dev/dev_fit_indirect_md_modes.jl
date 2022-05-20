include("lib_dev_fit.jl")

# -----------------------------------------------------------------------------
# Get data

file = jldopen("../numerical/data/sim_doe5.jld2", "r")

sims_all = gendata(file, 30, 0, 2200, remove_slow = 0.9);
sims = sims_all[[1, 2]]

# -----------------------------------------------------------------------------
# Approx functions

k̂!(u, p) = p[1]*u[2];
d̂!(u, p) = p[1]*u[1];
m̂!(u, p) = p[1];

kθ = [1e6];
mθ = [1e6];
dθ = [1e6];

# -----------------------------------------------------------------------------
# Differential equation setup

function diffeq!(du, u, p, t)
    return du[1] = (p[1] - d̂!(u[1], p[2])) / m̂!(u[1], p[3])
end

function diffeqmsd!(du, u, p, t)
    du[1] = (p[1] - k̂!(u, p[2]) - d̂!(u, p[3]) ) / m̂!(u, p[4])
    du[2] = u[1]
end

prob = ODEProblem(diffeqmsd!, [0.0, 0.0], (0.0, 20.0), [1e6, 1e6, 1e6, 1e6])
probs = ODEProblem(diffeq!, 0.0, (0.0, 20.0), [1e6, 1e6, 1e6])

# -----------------------------------------------------------------------------
# Train setup

function predict(p, nsim, prob; t=nothing)
    saveat = isnothing(t) ? nsim["trng"] : t
    tmp_prob = remake(prob; u0 = vcat(nsim["u0"], 0.0), p = vcat(nsim["Δp"], p))
    out = solve(tmp_prob, Tsit5(); saveat = saveat)
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

function lossout(p, prob)
    Ŷ = [predict(p, s, prob) for s in sims]
    ρ = mean([-cor(s["y"], ŷ) for (ŷ, s) in zip(Ŷ, sims)])
    rmse = sqrt(
        mean([mean((s["y"] .- ŷ) .^ 2) for (ŷ, s) in zip(Ŷ, sims)])
        )
    θ = sum(log10.(abs.(dθ))) + sum(log10.(abs.(mθ)))
    return ρ, rmse, θ
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

    ρ, rmse, pru  = lossout(θ, prob)
    hist[i] = ρ + rmse
    info[:ρ], info[:rmse], info[:θ] = ρ, rmse, pru

    # Stop criteria
    if stopcrit!(i, hist, α)
        break
    elseif (ρ < -(1-1e-4)) & (rmse < 1e-2)
        break
    end

    println(printit(i, hist[i], info))
end

nsim = sims[1]
plot(nsim["t"], predict(θ, nsim, prob))
plot!(nsim["t"], nsim["y"])

k = nsim["key"]
ts = file[k]["t"]
sim = predict(θ, nsim, prob, t=ts)

n = 600
plot(ts[1:n], sim[1:n])
plot!(ts[1:n], file[k]["V"][1:n])

plot(ts[1:n], file[k]["V"][1:n])

# using LsqFit

n = length(ts)
n = 4000
n = Int(round(length(ts)/2))

model(x, p) = @. p[1] - p[2]*exp(-p[3]*x)

tdata = collect(ts[1:n])

mdl = curve_fit(model, tdata, y, [1.0, 1.0, 1.0])
mdl.param

plot(tdata, model(tdata, mdl.param))
plot!(tdata, y)

n = length(ts)
n = 20000
plot(ts[1:n], y[1:n] - sim[1:n])
plot!(ts[1:n], model(ts[1:n], mdl.param) .- y[1:n])

# W = hcat(ts[1:n], ones(n, 1))
# y = file[k]["V"][1:n]
# a = W\y

# plot(ts[1:n], W*a)
# plot!(ts[1:n], y)


# plot!(ts[1:n], y - W*a)


# Fit residue

