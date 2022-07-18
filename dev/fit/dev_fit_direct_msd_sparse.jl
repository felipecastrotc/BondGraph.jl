using Random
using CUDA
include("lib_dev_fit.jl")

# For file history go to the fit-hist/dev_fit_indirect_md_ode.jl. 
# This is a cleaned file version of it, with some solved previous issues.

sim = Dict()

# -----------------------------------------------------------------------------
# Get data

# sim[:name] = "sim_doe3";
sim[:name] = "sim_doe5";
sim[:name] = "sim_doe6";
# sim[:npoitns], sim[:minRe], sim[:maxRe] = 30, 0, Inf
sim[:n], sim[:minRe], sim[:maxRe] = 30, 0, 2200
sim[:npoitns], sim[:minRe], sim[:maxRe] = 10, 0, Inf
sim[:n], sim[:minRe], sim[:maxRe] = 10, 3000, 10000

file = jldopen("../numerical/data/" * sim[:name] * ".jld2", "r")

sims_all = gendata(file, sim[:n], sim[:minRe], sim[:maxRe], remove_slow = 0.9);
sims = sims_all;

for s in sims_all
    # Reynolds without velocity
    s["sRe"] = s["all"][1]*s["all"][3]/s["all"][2]
    s["θ"] = [0.0, s["Iᵥ"], 0.0]
end

# -----------------------------------------------------------------------------
# Approx functions

include("lib_sparse_units.jl")


# -----------------------------------------------------------------------------
# Differential equation setup

function diffeq!(du, u, p, t)
    du[1] = (p[1] - d̂!(u, p[6:end])) / (m̂!(u, p[4])*p[2])
    du[2] = u[1]
end


prob = ODEProblem(diffeq!, [0.0, 0.0], (0.0, 20.0), γ)

# -----------------------------------------------------------------------------
# Train setup

function predict(p, nsim, prob, t=nothing)
    saveat = isnothing(t) ? nsim["trng"] : t
    tmp_prob = remake(prob; u0 = vcat(nsim["u0"], 0.0), p = vcat(nsim["Δp"],nsim["A"], nsim["θ"], nsim["all"][[1,3,4,5]], nsim["sRe"], p))
    out = solve(tmp_prob, Tsit5(); saveat=saveat, reltol=1e-8, abstol=1e-8)
    return vec(out[1, :])
end

# Define loss functions

function loss(p, prob)
    Ŷ = [predict(p, s, prob) for s in sims]
    ρ = mean([-cor(s["y"], ŷ) for (ŷ, s) in zip(Ŷ, sims)])
    rmse = sqrt(
        mean([mean(((s["y"] .- ŷ)./s["y"][end]) .^ 2) for (ŷ, s) in zip(Ŷ, sims)])
        )
    return ρ + rmse
end

function lossout(p, prob)
    Ŷ = [predict(p, s, prob) for s in sims]
    ρ = mean([-cor(s["y"], ŷ) for (ŷ, s) in zip(Ŷ, sims)])
    rmse = sqrt(
        mean([mean((s["y"] .- ŷ) .^ 2) for (ŷ, s) in zip(Ŷ, sims)])
        )
    return ρ, rmse
end

# Generate gradient functions
sims = sims_all[[1]]
∇γ = gradient((p) -> loss(p, prob), γ)

# -----------------------------------------------------------------------------
# Train

tik = 0.0

info = Dict(
    :s => 0,
    :ρ => 0.0,
    :rmse => 0.0,
    :tik => tik,
    :eta => 0.0,
)

idxs = 1:length(sims_all)
idxs = shuffle(idxs)

# i = [1]
i = 10
sims = deepcopy(sims_all[idxs[[i...]]])
sims[1]["Re"]

opt = ADAM(0.01)

it = 400000
α = 100
hist = zeros(it)
tik = @elapsed for i = 1:it
    
    ∇γ = gradient((p) -> loss(p, prob), γ)

    update!(opt, γ, ∇γ[1])

    ρc, rmse  = lossout(γ, prob)
    hist[i] = ρc + rmse
    info[:ρ], info[:rmse] = ρc, rmse

    # Stop criteria
    if stopcrit!(i, hist, α)
        break
    elseif (ρc < -(1-1e-4)) & (rmse < 1e-2)
        break
    elseif (rmse < 3e-2)
        break
    end
    # Print info it
    println(printit(i, hist[i], info))
end

i = 1
nsim = sims[i]
nsim["Re"]
plot(nsim["t"], nsim["y"])
# plot!(nsim["t"], predict(θ, nsim, prob), legend=false)
plot!(nsim["t"], predict(γ, nsim, prob), legend=false)

