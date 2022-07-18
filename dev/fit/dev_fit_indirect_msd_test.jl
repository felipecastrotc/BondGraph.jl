include("lib_dev_fit.jl")

# For file history go to the fit-hist/dev_fit_indirect_md_test.jl. 
# This is a cleaned file version of it, with some solved previous issues.

sim = Dict()

# -----------------------------------------------------------------------------
# Get data

# sim[:name] = "sim_doe3";
sim[:name] = "sim_doe5";
# sim[:npoitns], sim[:minRe], sim[:maxRe] = 30, 0, Inf
sim[:n], sim[:minRe], sim[:maxRe] = 30, 0, 2200

file = jldopen("../numerical/data/" * sim[:name] * ".jld2", "r")

sims_all = gendata(file, sim[:n], sim[:minRe], sim[:maxRe], remove_slow = 0.9);
sims = sims_all;

# -----------------------------------------------------------------------------
# Load model

filemdl = "./data/model_doe5.zdat"
m = loadzlib(filemdl)

x̄, x̃ = m[:x_mean], m[:x_std]
Ȳ = m[:y_mean]
Ỹ = m[:y_std]

m[:nnm]= cpu(m[:nnm])
m[:nnd]= cpu(m[:nnd])
m[:nns]= cpu(m[:nns])

# -----------------------------------------------------------------------------
# Approx functions

function convx(x)
    y =  copy(x)
    y[2, :] = log10.(y[2, :])
    y[3, :] = log10.(y[3, :])
    y[5, :] = log10.(y[5, :])
    y[7, :] = log10.(y[7, :])
    return y
end

# Mass model
function m̂!(u, p)
    return (m[:nnm]((convx(p) .- x̄) ./ x̃) .* Ỹ[:m] .+ Ȳ[:m])[1]
end;

# Damper model
function d̂!(u, p)
    return (m[:nnd]((convx(p) .- x̄) ./ x̃) .* Ỹ[:d] .+ Ȳ[:d])[1]*u[1]
end;

# Spring model
function k̂!(u, p)
    return (m[:nns]((convx(p) .- x̄) ./ x̃) .* Ỹ[:s] .+ Ȳ[:s])[1]*u[2]
end;

# -----------------------------------------------------------------------------
# Differential equation setup

function diffeqbg!(du, u, p, t)
    return du[1] = (p[1] - u[1] * p[2] )/ p[3]
end

function diffeq!(du, u, p, t)
    du[1] = (p[1] - k̂!(u, p[2]) - d̂!(u, p[2])) / m̂!(u, p[2])
    du[2] = u[1]
end

prob = ODEProblem(diffeq!, [0.0, 0.0], (0.0, 20.0), [0.0, zeros(7)])
probbg = ODEProblem(diffeqbg!, [0.0], (0.0, 20.0), [0.0, 0.0, 0.0])

# -----------------------------------------------------------------------------
# Set predict function

function predict(s, nsim, prob, t=nothing)
    saveat = isnothing(t) ? nsim["trng"] : t
    if s == :bg
        p = vcat(nsim["Δp"], nsim["Rᵥ"], nsim["Iᵥ"])
        tmp_prob = remake(prob; u0=nsim["u0"], p=p)
    else
        p = [nsim["Δp"], nsim["all"]]
        println(p)
        tmp_prob = remake(prob; u0 = vcat(nsim["u0"], 0.0), p = p)
    end
    out = solve(tmp_prob, Tsit5(); saveat=saveat, reltol=1e-8, abstol=1e-8)
    return vec(out[1, :]) ./ nsim["A"]
end

ksim = [k["key"] for k in sims]
klist = [findfirst(k .== ksim) for k in kl]

# Plot some simulation
plotsims(sims[[247]], [], prob)
plotsims(sims[[1212]], Dict(:bg => :bg, :s => []), Dict(:bg => probbg, :s => prob))

cors = []
for s in sims
    y = predict([], s, probt)
    if length(y) == length(s["y"])
        push!(cors, cor(y, s["y"]) - sqrt(mean((y - s["y"]).^2)))
    else
        push!(cors, 0.0)
    end
end

histogram(cors)

tmp = collect((1:length(cors)))[cors .> 0.9]
tmp = shuffle(tmp)[1:5]
plotsims(sims[tmp], [], probt)

tmp = [182, 550, 58, 297, 328, 858, 265]
# tmp = shuffle(klist)[1:4]
plotsims(sims[tmp], [], probt)
# [r["Re"] for r in sims[tmp]]



