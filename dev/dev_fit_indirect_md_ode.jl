include("lib_dev_fit.jl")
using BSON: @load

# -----------------------------------------------------------------------------
# Get data

# file = jldopen("../numerical/data/sim_doe1.jld2")
file = jldopen("../numerical/data/sim_doe4.jld2", "r")

sims_tmp = [s for s in sims_tmp if s["Rᵥ"] < 1e7];
sims_all = gendata(file, 10, 0, Inf, remove_slow=0.9);
# sims = sims_tmp
sims = sims_all;
Re = [i for (i, s) in enumerate(sims_all) if s["Re"] < 2200]
sims = sims_all[[6, 13, 22]]

# -----------------------------------------------------------------------------
# Approx functions

d̂!(u, p) = p[1];
m̂!(u, p) = p[1];

mθ = [1e6];
dθ = [1e6];

# -----------------------------------------------------------------------------
# Differential equation setup

function diffeq!(du, u, p, t)
    return du[1] = (p[1] - u[1] * d̂!(p[2], p[3])) / m̂!(p[4], p[5])
end

nsim = sims[1]
Δp = nsim["Δp"]
u0 = nsim["u0"]

diffeqd!(du, u, p, t) = diffeq!(du, u, [Δp, 1.0, p, 1.0, mθ], t)
diffeqm!(du, u, p, t) = diffeq!(du, u, [Δp, 1.0, dθ, 1.0, p], t)

probm = ODEProblem(diffeqm!, u0, (0.0, 20.0), mθ)
probd = ODEProblem(diffeqd!, u0, (0.0, 20.0), dθ)

# -----------------------------------------------------------------------------
# Train setup

function predict(p, nsim, prob)
    global Δp
    Δp = nsim["Δp"]
    tmp_prob = remake(prob; u0=nsim["u0"], p=p)
    return vec(solve(tmp_prob, Tsit5(); saveat=nsim["trng"])) ./ nsim["A"]
end

# Define loss functions
function loss(p, prob)
    return sqrt(
        mean([mean((s["y"] .- predict(p, s, prob)) .^ 2) for s in sims])
    )
end

# Generate gradient functions
∇d = gradient((p) -> loss(p, probd), dθ)
∇m = gradient((p) -> loss(p, probm), mθ)

# -----------------------------------------------------------------------------
# Train
# cfg = Dict()

tik = 0.0
c = 0.0
mlist, dlist, losslist, timelist = [], [], [], []
sims_train = sims_all[1:10]
sims_train = sims_all
for s in 1:length(sims_train)
    global sims, mθ, dθ

    mθ = [1e6]
    dθ = [1e6]

    sims = sims_train[[s]]

    nsim = sims[1]
    info = Dict(
        # :Iᵥ => Dict(:pred => () -> m̂!(1.0, mθ), :truth => nsim["Iᵥ"]),
        # :Rᵥ => Dict(:pred => () -> d̂!(1.0, dθ), :truth => nsim["Rᵥ"]),
        #:c => 0.0,
        :s => s,
        :time => tik,
    )

    λ = 1e-2
    α = 50
    c = 0.0
    it = 40000
    hist = zeros(it)
    tik = @elapsed for i in 1:it
        # i = 2
        ∇d = gradient((p) -> loss(p, probd), dθ)
        ∇m = gradient((p) -> loss(p, probm), mθ)
        
        λ = 5e-3  # Better than below according to a test with 20 samples
        # λ = i > 1 ? 10^(floor(log10(hist[i - 1])) - 1) : 5e-3
        # λ = i > 1 ? 10^(floor(log10(hist[i - 1])) - 1) : 1e-2
        # λ = i > 1 ? 10^(floor(log10(hist[i - 1])) - 2) : 5e-3
        # λ = 9e-3 
        # λ = i > 1 ? 10^(floor(log10(hist[i - 1])) - 1) : 9e-3

        dθ += -(λ * norm(dθ) / norm(∇d[1])) * ∇d[1]
        mθ += -(λ * norm(mθ) / norm(∇m[1])) * ∇m[1]

        hist[i] = loss(dθ, probd)

        # Stop criteria
        if stopcrit!(i, hist, α) break end
        #info[:c] = c

        # Print info it
        println(printit(i, hist[i], info))
    end

    push!(timelist, tik)
    push!(mlist, mθ[1])
    push!(dlist, dθ[1])

    i = findfirst(hist .== 0.0)
    push!(losslist, hist[i - 1])
end

extrema(losslist)
histogram(log10.(losslist))
mean(losslist)

# Save parameter
filename = "./data/fit_md_sim_doe4.jld2"
file = jldopen(filename, "w")

file["x"] = hcat([s["all"] for s in sims_train]...)
file["m"] = Float64.(mlist)
file["d"] = Float64.(dlist)
file["loss"] = Float64.(losslist)

close(file)


#FAST

# -----------------------------------------------------------------------------
# Approx functions

function d̂!(u, p)
    return (nd((u .- x̄) ./ x̃) .* ỹd .+ ȳd)[1]
end;

function m̂!(u, p)
    return (nm((u .- x̄) ./ x̃) .* ỹm .+ ȳm)[1]
end;


diffeqm!(du, u, p, t) = diffeq!(du, u, [Δp, p[1], 1.0, p[2], 1.0], t)


function predict(p, nsim, prob)
    global Δp
    Δp = nsim["Δp"]
    x = nsim["all"][[1, 2, 3, 4, 7]]
    x[2] = log10(x[2])
    x[3] = log10(x[3])
    x[5] = log10(x[5])
    p = [x, x]
    tmp_prob = remake(prob; u0=nsim["u0"], p=p)
    return vec(solve(tmp_prob, Tsit5(); saveat=nsim["trng"])) ./ nsim["A"]
end



predict(123, sims[123], probm)

d̂!(x[:, 1], 1)

x[:, 1]
sims[1]["all"]

yd[1]

idxs

x[:, 123]

k = id_train[j]
j = 392
k = id[j]
nsim = sims[k]
u = nsim["all"][[1, 2, 3, 4, 7]]
u[2] = log10(u[2])
u[3] = log10(u[3])
u[5] = log10(u[5])
x[:, j]- u

m̂!(u, 1.0)
nsim["Iᵥ"]
yd[k]
d̂!(u, 1.0)
nsim["Rᵥ"]
nsim["Re"]

plotsims(sims[[k]], [k], probm)

