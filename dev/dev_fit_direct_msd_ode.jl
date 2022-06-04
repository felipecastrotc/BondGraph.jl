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
sim[:npoitns], sim[:minRe], sim[:maxRe] = 10, 0, Inf
sim[:n], sim[:minRe], sim[:maxRe] = 10, 3000, 10000
sim[:n], sim[:minRe], sim[:maxRe] = 10, 0, 2200

file = jldopen("../numerical/data/" * sim[:name] * ".jld2", "r")

sims_all = gendata(file, sim[:n], sim[:minRe], sim[:maxRe], remove_slow = 0.9);
sims = sims_all;

for s in sims_all
    # Reynolds without velocity
    s["sRe"] = s["all"][1]*s["all"][3]/s["all"][2]
    s["θ"] = [0.0, s["Iᵥ"], 0.0]
end

# Flow functions
function colebrook(ϵ, d, Re)
    # Return the friction factor as 0 for reynolds number equals to 0
    if Re < 2300
        return min(64 / Re, 256)
    end
    # Niazkar approximation
    # https://link.springer.com/article/10.1007/s12205-019-2217-1
    A = -2 * log10((ϵ / d) / 3.7 + 4.5547 / (Re^0.8784))
    B = -2 * log10((ϵ / d) / 3.7 + 2.51 * A / (Re))
    C = -2 * log10((ϵ / d) / 3.7 + 2.51 * B / (Re))
    # Estimate f
    rhs = A - ((B - A)^2) / (C - 2 * B + A)
    f = 1 / (rhs^2)

    return f
end

# ρ, μ, d, L, ϵ

fric = [colebrook(s["all"][5], s["all"][3], s["Re"]) for s in sims_all]
histogram(fric, bins=50)

# xx = deepcopy(hcat([s["all"][1:5] for s in sims_all]...))
xx = deepcopy(hcat([vcat(s["all"][[3,5]], s["Re"]) for s in sims_all]...))
x = log10.(xx)

# x̄, x̃ = mean(x, dims=2), std(x, dims=2);
# ȳ, ỹ = mean(fric, dims=1)[1], std(fric, dims=1)[1];
X = (x .- x̄)./x̃;
Y = (fric .- ȳ)./ỹ

histogram(Y)
data = collect(zip(eachcol(X), Y))

# m = Chain(Dense(size(x, 1) => 5, tanh), Dense(5, 1))
m = Chain(Dense(size(x, 1) => 10, tanh), Dense(10, 1))
ps = Flux.params(m)

loss(x, y) = mean((y .- vec(m(x))).^2)

opt = ADAM(0.01)
opt = ADAM(0.001)
opt = ADAM(0.0001)

Flux.@epochs 400 Flux.train!(loss, ps, data, opt, cb = Flux.throttle(() -> println("Loss: ", sqrt(loss(X, Y))), 100))

loss(X, Y)

metrics(Y, vec(m(X)))

# -----------------------------------------------------------------------------
# Approx functions
# m = Chain(Dense(size(x, 1) => 15, tanh), Dense(15, 1))
γ, ml = Flux.destructure(m)

nn(u, p) = u[3] > 0 ? ml(p)((log10.(u) .- x̄)./x̃)[1]*ỹ + ȳ : 0.0

ρ = x -> x[1]; d = x -> x[2]; L = x -> x[3]; ϵ = x -> x[4]; re = x -> x[5]; γp = x -> x[6:end];

k̂!(u, p) = p[1] * u[2];
# d̂!(u, p) = nn(vcat(d(p), ϵ(p), re(p)), γp(p))*(ρ(p)*L(p)/d(p)) * (u[1])^2;

ap(u, p) = u[1] > 0 ? p[1]*exp(sum(p[2:end].*log.(u))) : 0.0
d̂!(u, p) =  ap(vcat( re(p)*abs(u[1]), d(p), ϵ(p), ρ(p), L(p)), γp(p))*(u[1]^2);
# function d̂!(u, p)
#     x = ap(vcat( re(p)*abs(u[1]), d(p), ϵ(p), ρ(p), L(p)), γp(p))*(u[1]^2);
#     # println(x)
#     return x
# end
m̂!(u, p) = p[1];

plot(nsim["t"], nsim["y"])
plot!(nsim["t"], predict(γ, sims[1], prob))

cvt = Dict(:m => x -> x[3], :s => x -> x[1], :d => x -> x[2])

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
i = 1:15
sims = deepcopy(sims_all[idxs[[i...]]])
sims[1]["Re"]
[s["Re"] for s in sims]

γ = ones(6)
opt = ADAM(0.001)

it = 400000
α = 1000
hist = zeros(it)
tik = @elapsed for i = 1:it
    
    ∇γ = gradient((p) -> loss(p, prob), γ)
    H = hessian((p) -> loss(p, prob), γ)

    λ = 1e-2;
    γ += -((H + λ*I)^-1)*∇γ[1]
    # update!(opt, γ, ∇γ[1])

    ρc, rmse  = lossout(γ, prob)
    hist[i] = ρc + rmse
    info[:ρ], info[:rmse] = ρc, rmse

    # Stop criteria
    if stopcrit!(i, hist, α)
        break
    elseif (ρc < -(1-1e-4)) & (rmse < 1e-2)
        break
    elseif (rmse < 1e-3)
        break
    end
    # Print info it
    println(printit(i, hist[i], info))
end

i = 300
nsim = sims_all[i]
# nsim = sims[i]
nsim["Re"]
plot(nsim["t"], nsim["y"])
# plot!(nsim["t"], predict(θ, nsim, prob), legend=false)
plot!(nsim["t"], predict(γ, nsim, prob), legend=false)



# -----------------------------------------------------------------------
# MSD

plot(rand(1000))
scatter(rand(1000), rand(1000))
histogram(randn(1000))

function diffeq!(du, u, p, t)

    m₁, m₂ = p[1], p[2];
    k₁, k₂ = p[3], p[4];
    c₁, c₂ = p[5], p[6];
    F₁, F₂ = p[7], p[8];
    ẋ₁, x₁ = u[1], u[2];
    ẋ₂, x₂ = u[3], u[4];

    du[1] = (F₁ - (k₁ + k₂)*x₁ + k₂*x₂ - (c₁ + c₂)*ẋ₁ + c₂*ẋ₂)/m₁
    du[2] = ẋ₁
    du[3] = (F₂ + abs(k₂)*x₁ - abs(k₂)*x₂ + abs(c₂)*ẋ₁ -abs(c₂)*ẋ₂)/abs(m₂)
    du[4] = ẋ₂
end

m = [4.0, 2.0];
k = [1.0, 0.5];
c = [0.05, 0.9];
F = [0.0, 0,0];
# u0 = [0.0, 0.2, 1.0, 0.0];
# u0 = [0.0, 0.0, 0.0, 0.0];
u0 = [1.0, 1.0, 0.0, 2.0];

ζ = c[1]/(2*sqrt(m[1]*k[1]))
ζ = γ[5]/(2*sqrt(γ[1]*γ[3]))

sqrt(γ[3]/γ[1])
sqrt(k[1]/m[1])

ts = 100.0
prob = ODEProblem(diffeq!, u0, (0.0, ts), vcat(m, k, c, F))
sol = solve(prob)
plot(sol)

trng = 0:0.01:ts
plot(trng2, sol(trng2)[var[1], :])
plot!(trng2, predict(γ, prob, trng2)[1, :], legend=false)


# -----------------------------------------------------------------------------
# Train setup

function predict(p, prob, t=nothing, full=false)
    saveat = isnothing(t) ? trng : t
    # tmp_prob = remake(prob; p = vcat(m[1], p[1], k[1], p[2], c[1], p[3], F))
    tmp_prob = remake(prob; p = vcat(p, F))
    out = solve(tmp_prob, Tsit5(); saveat=saveat, reltol=1e-8, abstol=1e-8)
    if full
        return out
    else
        return out[var, :]
        # return out[[1,2], :]
    end
    # return out[[1, 2], :]
end

function loss(p, prob, Y, t=nothing)
    Ŷ = predict(p, prob, t)
    # ρ = -cor(Ŷ, Y)
    ρ = -mean(cor.(eachrow(Y), eachrow(Ŷ)))
    rmse = sqrt(mean(((Y - Ŷ)./maximum(Y)).^2))
    return ρ + rmse + sum(1.0./exp.(3*p))
end

function lossout(p, prob, Y)
    Ŷ = predict(p, prob)
    # ρ = -cor(Ŷ, Y)
    ρ = -mean(cor.(eachrow(Y), eachrow(Ŷ)))
    rmse = sqrt(mean((Y - Ŷ).^2))
    return ρ, rmse, sum(1.0./exp.(3*p))
end

# γ = ones(size(vcat(m, k, c, F)))
# Generate gradient functions
trng = 0:1:100
var = [2]
Y = sol(trng)[var, :]
γ = ones(6)
∇γ = gradient((p) -> loss(p, prob, Y, trng), γ)

# -----------------------------------------------------------------------------
# Train

γ
tik = 0.0

info = Dict(
    :s => 0,
    :ρ => 0.0,
    :rmse => 0.0,
    :tik => tik,
    :eta => 0.0,
)

γ₀ = [0.1, 0.1, 0.1]
γ₀ = ones(6)*0.1
it = 10000
mit = 1
α = 10000
ρc, rmse, κ = Inf, Inf, Inf;
opt = ADAM(0.1)
hist = zeros(floor(Int, it/mit))
tik = @elapsed for i = 1:it
    ∇γ = gradient((p) -> loss(p, prob, Y), γ)[1]
    # # h = [1,2]
    # # h = [1,4]
    # h = [2,5]
    # # h = [3,6]
    # ζ = γ[h]
    # ∇ζ = ∇γ[h]
    # update!(opt, ζ, ∇ζ)
    # γ[h] = ζ
    update!(opt, γ, ∇γ)

    if i % mit == 0.0
        z = floor(Int, i / mit)
        ρc, rmse, κ  = lossout(γ, prob, Y)
        hist[z] = ρc + rmse
        info[:ρ], info[:rmse] = ρc, rmse

        # Stop criteria
        if stopcrit!(i, hist, α)
            break
        # elseif (ρc < -(1-1e-4)) & (rmse < 1e-2)
        elseif (ρc < -(1-1e-5)) & (rmse < 1e-3)
            println(printit(i, hist[z], info))
            break
        end

        # Print info it
        println(printit(i, hist[z], info))
    end
end

γ
trng2 = 0:0.01:100
# plot(trng2, sol(trng2)[1, :])
plot(trng2, sol(trng2)[var[1], :])
plot!(trng2, predict(γ, prob, trng2)[1, :], legend=false)
# plot!(trng2, predict(γ, prob, trng2)[2, :], legend=false)

plot(sol(trng2))
plot!(predict(γ, prob, trng2, true), legend=false)



