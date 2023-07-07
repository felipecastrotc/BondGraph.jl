include("lib_dev_fit.jl")
using BSON: load

# -----------------------------------------------------------------------------
# Get data

# file = jldopen("../numerical/data/sim_doe1.jld2")
file = jldopen("../numerical/data/sim_doe4.jld2", "r")

sims_all = gendata(file, 10, 0, Inf, remove_slow=0.9);
sims_tmp = [s for s in sims_tmp if s["Rᵥ"] < 1e7];
# sims = sims_tmp
sims = sims_all;


# -----------------------------------------------------------------------------
# Approx functions

@load "./data/model_sim_doe4.bson" nn
@load "./data/model_sim_doe4.bson" nn

# Mass model

function m̂!(u, p)
    return (nn((u .- mean(x; dims=2)) ./ std(x; dims=2)) .* std(y) .+ mean(y))[1]
end;

function d̂!(u, p)
    return (nn((u .- mean(x; dims=2)) ./ std(x; dims=2)) .* std(y) .+ mean(y))[1]
end;


# -----------------------------------------------------------------------------
# Differential equation setup

function diffeq!(du, u, p, t)
    return du[1] = (p[1] - u[1] * d̂!(p[2], p[3])) / m̂!(p[4], p[5])
end

nsim = sims[1]
Δp = nsim["Δp"]
u0 = nsim["u0"]

diffeqt!(du, u, p, t) = diffeq!(du, u, [p[1], p[2], 1.0, p[3], 1.0], t)
probt = ODEProblem(diffeqm!, u0, (0.0, 20.0), [Δp, xd_in, xm_in])

# -----------------------------------------------------------------------------
# Solution and loss

function predict(p, nsim, prob)
    p = vcat(nsim["Δp"], )
    tmp_prob = remake(prob; u0=nsim["u0"], p=p)
    return vec(solve(tmp_prob, Tsit5(); saveat=nsim["trng"])) ./ nsim["A"]
end

function loss(p, prob)
    return sqrt(
        mean([mean((s["y"] .- predict(p, s, prob)) .^ 2) for s in sims])
    )
end

plotsims(sims, [], probt)
