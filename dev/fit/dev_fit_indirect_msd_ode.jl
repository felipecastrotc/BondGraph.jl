using Random

include("lib_dev_fit.jl")

# For file history go to the fit-hist/dev_fit_indirect_md_ode.jl. 
# This is a cleaned file version of it, with some solved previous issues.

sim = Dict()

# -----------------------------------------------------------------------------
# Get data

# sim[:name] = "sim_doe3";
sim[:name] = "sim_doe5";
# sim[:npoitns], sim[:minRe], sim[:maxRe] = 30, 0, Inf
sim[:n], sim[:minRe], sim[:maxRe] = 30, 4000, 6000
sim[:n], sim[:minRe], sim[:maxRe] = 30, 0, 2200

file = jldopen("../numerical/data/" * sim[:name] * ".jld2", "r")

sims_all = gendata(file, sim[:n], sim[:minRe], sim[:maxRe], remove_slow = 0.9);
sims = sims_all;

# -----------------------------------------------------------------------------
# Approx functions

# k̂!(u, p) = p[1] * u[2]^2;
k̂!(u, p) = p[1] * u[2];
k̂!(u, p) = 0 * u[2];
d̂!(u, p) = p[1] * u[1];
m̂!(u, p) = p[1];

θ = [1e6, 1e6, 1e6];
cvt = Dict(:m => x -> x[3], :s => x -> x[1], :d => x -> x[2])

# -----------------------------------------------------------------------------
# Differential equation setup

function diffeq!(du, u, p, t)
    du[1] = (p[1] - k̂!(u, p[2]) - d̂!(u, p[3])) / m̂!(u, p[4])
    du[2] = u[1]
end

prob = ODEProblem(diffeq!, [0.0, 0.0], (0.0, 20.0), θ)

# -----------------------------------------------------------------------------
# Train setup

function predict(p, nsim, prob, t=nothing)
    saveat = isnothing(t) ? nsim["trng"] : t
    tmp_prob = remake(prob; u0 = vcat(nsim["u0"], 0.0), p = vcat(nsim["Δp"], p))
    out = solve(tmp_prob, Tsit5(); saveat=saveat, reltol=1e-8, abstol=1e-8)
    return vec(out[1, :]) ./ nsim["A"]
end

# Define loss functions

function loss(p, prob)
    return sqrt(
        mean([mean((s["y"] .- predict(p, s, prob)) .^ 2) for s in sims]),
    )
end

function loss(p, prob)
    Ŷ = [predict(p, s, prob) for s in sims]
    ρ = mean([-cor(s["y"], ŷ) for (ŷ, s) in zip(Ŷ, sims)])
    rmse = sqrt(
        mean([mean((s["y"] .- ŷ) .^ 2) for (ŷ, s) in zip(Ŷ, sims)])
        )
    # θ = sum(log10.(abs.(dθ))) + sum(log10.(abs.(mθ)))
    # return ρ + rmse + 0.1*θ
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
θ = [1e6, 1e6, 1e6];
∇θ = gradient((p) -> loss(p, prob), θ)

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
θlist, losslist, timelist, keylist, xlist = [], [], [], [], []
idxs = 1:1
for (j, s) in enumerate(idxs)
    global sims, θ

    sims = sims_all[[s]]
    info[:s] = j

    θ = [1e6, 1e6, 1e6]

    it = 400000
    α = 1000
    hist = zeros(it)
    tik = @elapsed for i = 1:it
        
        ∇θ = gradient((p) -> loss(p, prob), θ)

        λ = 5e-3
        θ += -(λ * norm(θ) / norm(∇θ[1])) * ∇θ[1]

        ρ, rmse  = lossout(θ, prob)
        hist[i] = ρ + rmse
        info[:ρ], info[:rmse] = ρ, rmse

        # Stop criteria
        if stopcrit!(i, hist, α)
            break
        elseif (ρ < -(1-1e-4)) & (rmse < 1e-2)
            break
        end

        # Print info it
        println(printit(i, hist[i], info))
    end
    push!(θlist, θ)
    push!(xlist, sims[1]["all"])
    push!(keylist, sims[1]["key"])

    i = findfirst(hist .== 0.0)
    push!(losslist, hist[i-1])
    push!(timelist, tik)

    info[:eta] = (length(idxs) - j) * (sum(timelist)/length(timelist))
    info[:tik] = tik
end

θ
# Visualize fit
i = 1
θ = θlist[i]
nsim = sims_all[i]
nsim["Re"]
plot(nsim["t"], nsim["y"])
plot!(nsim["t"], predict(θ, nsim, prob), legend=false)


function predict(s, nsim, prob, t=nothing)
    saveat = isnothing(t) ? nsim["trng"] : t
    if s == :bg
        p = vcat(nsim["Δp"], nsim["Rᵥ"], nsim["Iᵥ"])
        tmp_prob = remake(prob; u0=nsim["u0"], p=p)
    else
        p = [nsim["Δp"], s]
        println(p)
        tmp_prob = remake(prob; u0 = vcat(nsim["u0"], 0.0), p = p)
    end
    out = solve(tmp_prob, Tsit5(); saveat=saveat, reltol=1e-8, abstol=1e-8)
    return vec(out[1, :]) ./ nsim["A"]
end


plt = plot()
# Color palette selection
c = palette(:default)[i % 15]
# Plot the reference
plt = plot!(nsim["t"], nsim["y"]; ls=:dash, lc=c, label="Ref. PDE", xaxis="Time (s)", yaxis="Velocity (m/s)")
ls = [:solid :dashdot :dashdotdot][1]
tmp_prob = remake(prob; u0 = vcat(nsim["u0"], 0.0), p =vcat(nsim["Δp"], θ))
out = solve(tmp_prob, Tsit5(); saveat=nsim["trng"], reltol=1e-8, abstol=1e-8)
ys = vec(out[1, :]) ./ nsim["A"]
label = ("Sim. fit")
plot!(nsim["t"], ys; ls=ls, label=label)
ls = [:solid :dashdot :dashdotdot][2]
tmp_prob = remake(probbg; u0 = vcat(nsim["u0"], 0.0), p =vcat(nsim["Δp"], nsim["Rᵥ"], nsim["Iᵥ"]))
out = solve(tmp_prob, Tsit5(); saveat=nsim["trng"], reltol=1e-8, abstol=1e-8)
yb= vec(out[1, :]) ./ nsim["A"]
label = ("Sim. BG")
plot!(nsim["t"], yb; ls=ls, label=label, legend=:bottomright)
plot!(xlim=[0.1, 0.8], ylim=[3.0, 6.0])


metrics(nsim["y"], ys)
metrics(nsim["y"], yb)

# -----------------------------------------------------------------------------
# Save θ
filename = "./data/fit_msd_sim_doe5.jld2"
saveθ(filename, xlist, θlist, losslist, keylist, cvt)