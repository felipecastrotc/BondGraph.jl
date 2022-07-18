include("lib_dev_fit.jl")
using BSON: @save

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
        :c => 0.0,
        :s => s,
        :time => tik,
    )

    λ = 1e-2
    α = 50
    c = 0
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
        if stopcrit!(c, i, hist, α) break end
        info[:c] = c

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