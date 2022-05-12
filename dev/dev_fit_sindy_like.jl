using BSON: @load
using LsqFit
using SymPy

include("lib_dev_fit.jl")


# -----------------------------------------------------------------------------
# Get data

# file = jldopen("../numerical/data/sim_doe1.jld2")
file = jldopen("../numerical/data/sim_doe3.jld2", "r")

sims_all = gendata(file, 10, 0, Inf; remove_slow=0.9);
sims_all = gendata(file, 10, 1800, 4000, remove_slow=0.9);
sims_all = gendata(file, 30, 500, 1500; remove_slow=0.9);
sims = sims_all[[1, 2]]
sims = sims_all

# x = vcat([s["all"] for s in sims_all]'...)';
# x̄, x̃ = mean(x; dims = 2), std(x; dims = 2)
# Xs = (x .- x̄) ./ x̃

# for (i, s) in enumerate(sims_all)
#     s["scl"] = Xs[:, i]
# end

# sims = sims_all

# -----------------------------------------------------------------------------
# Like SINDy

# [1, 2, 3, 4, 5, 6, 7]
# [ρ, μ, d, L, ϵ, a, A]

m̂!(x, p) = p[1]*(x[1]*x[4]/(x[7])) + p[2]*(x[2]*x[4]/(x[3]^4)) + p[3]*(x[1]*x[2]*x[3]/(x[6]*x[7])) + p[4]*(x[3]*x[6]/x[1]);
d̂!(x, p) = p[1]*(x[1]*x[4]/(x[7])) + p[2]*(x[2]*x[4]/(x[3]^4)) + p[3]*(x[1]*x[2]*x[3]/(x[6]*x[7])) + p[4]*(x[3]*x[6]/x[1]);
128/π

dθ = ones(4)
mθ = ones(4)

# -----------------------------------------------------------------------------
# Differential equations
# Train setup

dθ = [0, 0, 0, 1.0]
mθ = [0, 0, 0, 1.0]

nsim = sims[1]
u = nsim["u0"]
p = [nsim["Δp"], dp, dθ, mp, mθ]

d̂!(vcat(p[2], u), p[3])
# min(

function diffeq!(du, u, p, t)
    # return du[1] = (p[1] - u[1] * d̂!(vcat(p[2], u), p[3])) / max(m̂!(vcat(p[4], u), p[5]), 2e3)
    return du[1] = (p[1] - u[1] * d̂!(vcat(p[2], u), p[3])) / m̂!(vcat(p[4], u), p[5])
end

nsim = sims[1]
mp = nsim["all"]
dp = nsim["all"]
Δp = nsim["Δp"]
u0 = nsim["u0"]

diffeqd!(du, u, p, t) = diffeq!(du, u, [Δp, dp, p, mp, mθ], t)
diffeqm!(du, u, p, t) = diffeq!(du, u, [Δp, dp, dθ, mp, p], t)

probm = ODEProblem(diffeqm!, u0, (0.0, 20.0), mθ)
probd = ODEProblem(diffeqd!, u0, (0.0, 20.0), dθ)

# Test
function predict(p, data, prob)
    global mp, dp, Δp
    tmp_prob = remake(prob; u0=data["u0"], p=p)
    mp = data["all"]
    dp = data["all"]
    # mp = data["scl"]
    # dp = data["scl"]
    # mp = data["mp"]
    # dp = data["dp"]
    # dp = vcat(data["all"], [1.0], 1.0 ./ data["all"])
    Δp = data["Δp"]
    return vec(solve(tmp_prob, Tsit5(); saveat=data["trng"])) ./ data["A"]
end

# -----------------------------------------------------------------------------
# Define loss functions

# function loss(p, prob)
#     return sqrt(
#         mean([mean((s["y"] .- predict(p, s, prob)) .^ 2) for s in sims])
#     ) + sum(log10.(abs.(dθ.^2))) + sum(log10.(abs.(mθ.^2)))
# end
# function loss(p, prob)
#     return sqrt(
#         mean([mean((s["y"] .- predict(p, s, prob)) .^ 2) for s in sims])
#     )
# end

# function loss(p, prob)
#     return sqrt(
#         mean([mean((s["y"] .- predict(p, s, prob)) .^ 2) for s in sims])
#     ) + 3*norm(dθ./norm(dθ), 1) + 3*norm(mθ./norm(mθ), 1)
# end
# function loss(p, prob)
#     return 0.9*mean([-cor(s["y"], predict(p, s, prob)) for s in sims]) + 0.01*sqrt(
#             mean([mean((s["y"] .- predict(p, s, prob)) .^ 2) for s in sims])
#         )
# end
function lossm(p, prob)
    return mean([-cor(s["y"], predict(p, s, prob)) for s in sims])
end

function lossd(p, prob)
    return sqrt(
            mean([mean(((s["y"] .- predict(p, s, prob))/s["y"][end]) .^ 2) for s in sims])
        )
end


# function loss(p, prob)
    # return mean([-cor(s["y"], predict(p, s, prob)) for s in sims])
# end

# -----------------------------------------------------------------------------
# Generate gradient functions
# ∇d = gradient((p) -> loss(p, probd), dθ)
# ∇m = gradient((p) -> loss(p, probm), mθ)
∇d = gradient((p) -> lossd(p, probd), dθ)
∇m = gradient((p) -> lossm(p, probd), mθ)

# -----------------------------------------------------------------------------
# Training
info = Dict(
    :lm => 0.0,
)


i = 1
λ = 5e-3
α = 50

sims = sims_all[[1, 2,3,4,5,6]]
# sims = sims_all[[3, 4]]
sims = sims_all[[2, 3]]
[s["Re"] for s in sims]

dθ = [0.0, 40, 0.0, 0]
mθ = [0.98, 0, 0, 0.0]

dθ = [0.0, 128/π, 0.0, 0]
mθ = [1.0, 0, 0, 0.0]

loss(dθ, probd)

dθ = Float64.(rand(4) .> 0.35)
mθ = Float64.(rand(4) .> 0.35)
dθ = rand(4)
mθ = rand(4)
# mθ = ones(4)
optm = ADAM(0.1);
optd = ADAM(0.1);
# plotsims(sims[[1]], dθ, probd)
it = 100000
hist = zeros(it)
for i in 1:it
    ∇d = gradient((p) -> lossd(p, probd), dθ)
    ∇m = gradient((p) -> lossm(p, probm), mθ)

    # Hd = Zygote.hessian((p) -> lossd(p, probd), dθ)
    # Hm = Zygote.hessian((p) -> lossm(p, probm), mθ)

    # dθ += -Hd*∇d[1]*1e-1
    # if norm(∇m[1]) > 0.0
    #     mθ += -Hm*∇m[1]*1e-1
    # else
    #     mθ += -Hd*∇d[1]*1e-1
    # end
    # hist[i] = loss(dθ, probd)
    # if norm(∇m[1]) > 0.0
    #     update!(optm, mθ, ∇m[1])
    # else
    #     update!(optm, dθ, ∇d[1])
    # end
    update!(optm, mθ, ∇m[1])
    update!(optd, dθ, ∇d[1])
    
    # # λ = i > 1 ? 10^(floor(log10(abs(hist[i - 1]))) - 1) : 5e-3
    # λ = i > 1 ? 10^(floor(log10(abs(hist[i - 1]))) - 2) : 5e-3
    # # λ = i > 1 ? 10^(floor(log10(abs(hist[i - 1]))) - 3) : 5e-3
    # # λ = i > 1 ? 10^(floor(log10(abs(hist[i - 1]))) - 4) : 5e-3
    # # λ = 1e-4
    # # λ = 1e-6
    # # (norm(dθ) / norm(∇d[1])) * ∇d[1]
    # # dθ += -(λ * norm(dθ) / norm(∇d[1])) * ∇d[1]
    # if norm(∇m[1]) > 0.0
    #     mθ += -(λ * norm(mθ) / norm(∇m[1])) * ∇m[1]
    # else
    #     mθ += -(λ * norm(dθ) / norm(∇d[1])) * ∇d[1]
    # end

    hist[i] = lossd(dθ, probd)
    info[:lm] = lossm(mθ, probm)
    
    # Print info it
    println(printit(i, hist[i], info))
    # Plot 
    plotit(i, sims[1], () -> predict(dθ, sims[1], probd); every=50)
    # plotit(i, nsim, () -> predict(dθ, nsim, probd); every=10)
end

dθ
mθ

plotsims(sims, dθ, probd)
plotsims(sims[[1]], dθ, probd)
plotit(1, sims[1], () -> predict(dθ, sims[1], probd); every=10)

nsim = sims[1]
nsim["Re"]
λ = nsim["y"][end]/predict(dθ, nsim, probd)[end]
cor(predict(dθ, nsim, probd), nsim["y"])
plot(nsim["y"])
# plot!(predict(dθ, nsim, probd).*λ)
plot!(predict(dθ, nsim, probd).*λ)


# Trying LsqFit

Xs

function fitcp(X, y)
    order, rnk = l
    # Set LstFit model
    model(X̂, p̂) = CP(X̂, p̂, rnk, order)

    return mdl
end

m̂!(x, mθ)

using LsqFit
sdy(x, p) = @. p[1]*(x[1, :]*x[4, :]/(x[7, :])) + p[2]*(x[2, :]*x[4, :]/(x[3, :]^4)) + p[3]*(x[1, :]*x[2, :]*x[3, :]/(x[6, :]*x[7, :])) + p[4]*(x[3, :]*x[6, :]/x[1, :]);

x


W = hcat((x[1, :].*x[4, :]./(x[7, :])), (x[2, :].*x[4, :]./(x[3, :].^4)), (x[1, :].*x[2, :].*x[3, :]./(x[6, :].*x[7, :])), (x[3, :].*x[6, :]./x[1, :]))


mθ = W\ym
dθ = W\yd

dθ = [0, 128/π, 0, 0]
scatter(W*mθ, ym)
scatter(W*dθ, yd)

# Initial guess
p₀ = ones(4)
# Fit
mdl = curve_fit(sdy, x, yd, p₀)
mdl.param

(sdy(x, mdl.param) .- ym)./ym
(sdy(x, mdl.param) .- yd)./yd
