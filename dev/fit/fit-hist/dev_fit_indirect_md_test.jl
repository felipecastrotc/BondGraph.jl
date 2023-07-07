include("lib_dev_fit.jl")
using Serialization
using Statistics

# -----------------------------------------------------------------------------
# Get data

# file = jldopen("../numerical/data/sim_doe1.jld2")
file = jldopen("../numerical/data/sim_doe4.jld2", "r")

sims_all = gendata(file, 30, 0, Inf, remove_slow=0.9);
Re = [i for (i, s) in enumerate(sims_all) if s["Re"] > 3000]
# Re = [i for (i, s) in enumerate(sims_all) if s["Re"] < 2200]
sims = sims_all[Re]  # turb

sims = [s for s in sims if s["key"] in kl]

# -----------------------------------------------------------------------------
# Approx functions

# cfg = deserialize("./data/model_sim_vars_doe4.dat")
# nd = deserialize("./data/model_sim_res_doe4.dat")
# nm = deserialize("./data/model_sim_ine_doe4.dat")
nd = deserialize("./data/model_sim_rest_doe4.dat")
nm = deserialize("./data/model_sim_inet_doe4.dat")
cfg = deserialize("./data/model_sim_varst_doe4.dat")

x̄, x̃ = cfg[:x_mean], cfg[:x_std]
ȳd, ỹd = cfg[:yd_mean], cfg[:yd_std]
ȳm, ỹm = cfg[:ym_mean], cfg[:ym_std]

# Some variables are log based 
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
    return (nm((convx(u) .- x̄) ./ x̃) .* ỹm .+ ȳm)[1]
end;

# Damper model
function d̂!(u, p)
    return (nd((convx(u) .- x̄) ./ x̃) .* ỹd .+ ȳd)[1]
end;

# -----------------------------------------------------------------------------
# Differential equation setup

function diffeqbg!(du, u, p, t)
    return du[1] = (p[1] - u[1] * p[2] )/ p[3]
end

function diffeql!(du, u, p, t)
    return du[1] = (p[1] - u[1] * d̂!(p[2], 1.0)) / m̂!(p[2], 1.0)
end

function diffeqt!(du, u, p, t)
    return du[1] = (p[1] - (u[1]^2) * d̂!(p[2], 1.0)) / m̂!(p[2], 1.0)
end

nsim = sims[1]
Δp = nsim["Δp"]
u0 = nsim["u0"]

probbg = ODEProblem(diffeqbg!, u0, (0.0, 20.0), [Δp, nsim["Rᵥ"], nsim["Iᵥ"]])
probl = ODEProblem(diffeql!, u0, (0.0, 20.0), [Δp, nsim["all"]])
probt = ODEProblem(diffeqt!, u0, (0.0, 20.0), [Δp, nsim["all"]])

# -----------------------------------------------------------------------------
# Solution and loss

function predict(s, nsim, prob)
    if s == :bg
        p = [nsim["Δp"], nsim["Rᵥ"], nsim["Iᵥ"]]
    else
        p = [nsim["Δp"], nsim["all"]]
    end
    tmp_prob = remake(prob; u0=nsim["u0"], p=p)
    return vec(solve(tmp_prob, Tsit5(); saveat=nsim["trng"])) ./ nsim["A"]
end


ksim = [k["key"] for k in sims]
klist = [findfirst(k .== ksim) for k in kl]


# plotsims(sims[[297, 1212, 906, 599]], [], probl)
plotsims(sims[[247]], [], probt)

plotsims(sims[[1212]], Dict(:bg => :bg, :s => []), Dict(:bg => probbg, :s => probl))


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


