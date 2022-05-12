include("lib_dev_fit.jl")
using BSON: @load

function ArmijoLineSearch(f, d, α0; ρ=0.5, c1=1e-4)

    ϕ0 = f(0)
    dϕ0 = -1.0.*(d.^2)

    ϕ_a0 = f(α0)
    k = 1
    while (ϕ_a0 > (ϕ0 + c1*α0*dϕ0)) || (k < 10)
        α0 = α0 * ρ
        ϕ_a0 = f(α0)
        k += 1
    end
    
    return α0
end


function idxloss(p, loss, i, θ, prob)
    if i == 1
        return loss(vcat(p, θ[2:end]), prob)
    else
        return loss(vcat(θ[1:i-1], p, θ[i+1:end]), prob)
    end
end


# -----------------------------------------------------------------------------
# Get data

# file = jldopen("../numerical/data/sim_doe1.jld2")
file = jldopen("../numerical/data/sim_doe4.jld2", "r")

# sims_tmp = [s for s in sims_tmp if s["Rᵥ"] < 1e7];
sims_all = gendata(file, 30, 8000, 10000, remove_slow = 0.9);
sims_all = gendata(file, 30, 0, Inf, remove_slow = 0.9);
# sims = sims_tmp
sims = sims_all;
Re = [i for (i, s) in enumerate(sims_all) if s["Re"] > 3500]
Re = [i for (i, s) in enumerate(sims_all) if s["Re"] < 2200]
sims = sims_all[[1, 3, 8]]
sims = sims_all[[9, 12, 14, 15]]  # lam
sims = sims_all[[3, 47, 55]]  # turb
sims = sims_all[Re]  # turb

findfirst("1" .== [k["key"] for k in sims])

# -----------------------------------------------------------------------------
# Approx functions

# d̂!(u, p) = p[1];
# m̂!(u, p) = p[1];

# mθ = [1e6];
# dθ = [1e6];

# d̂!(u, p) = p[1]*u[1]^2;
d̂!(u, p) = p[1]*u[1];
d̂!(u, p) = Float64(nn((Float32.(xx) .- x̄)./x̃)[1]*ỹd + ȳd)*u[1];

# d̂!(u, p) = p[1]*u[1] + p[2]*u[1]^2;
# d̂!(u, p) = p[1] + p[2]*u[1] + p[3]*u[1]^2;
# m̂!(u, p) = p[1] + p[2]*u[1] + p[3]*u[1]^2;
# d̂!(u, p) = p[3]*u[1]^2;
m̂!(u, p) = p[1];

# mθ = [1e6, 0.0, 0.0];
# dθ = [0.0, 0.0, 1e6];
mθ = [1e6];
dθ = [1e6];
# dθ = [1e6, 0.0];
# mθ = [1e6, 0.0, 0.0];
# dθ = [0.0, 0.0, 1e6];

# -----------------------------------------------------------------------------
# Differential equation setup

# function diffeq!(du, u, p, t)
#     return du[1] = (p[1] - u[1] * d̂!(p[2], p[3])) / m̂!(p[4], p[5])
# end

function diffeq!(du, u, p, t)
    return du[1] = (p[1] - d̂!(u[1], p[3])) / m̂!(u[1], p[5])
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
    global Δp, xx
    Δp = nsim["Δp"]
    # nsim = sims[376]
    xx = copy(nsim["all"])
    xx[2] = log10.(xx[2])
    xx[3] = log10.(xx[3])
    xx[5] = log10.(xx[5])
    xx[7] = log10.(xx[7])
    tmp_prob = remake(prob; u0 = nsim["u0"], p = p)
    return vec(solve(tmp_prob, Tsit5(); saveat = nsim["trng"])) ./ nsim["A"]
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
    θ = sum(log10.(abs.(dθ))) + sum(log10.(abs.(mθ)))
    return ρ, rmse, θ
end


# loss(dθ, probd)
# Generate gradient functions
∇d = gradient((p) -> loss(p, probd), dθ)
∇m = gradient((p) -> loss(p, probm), mθ)

# -----------------------------------------------------------------------------
# Train
# cfg = Dict()
using Random

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

mlist, dlist, losslist, timelist, klist, xlist = [], [], [], [], [], []
sims_train = sims_all
sims_train = sims_all[[9, 12, 14, 15]]
sims_train = sims_all[[3, 47, 55, 63]]
sims_train = sims_all[Re]
# sims_train = sims_all[1:10]

mθ = [1e6]
dθ = [1e6]

# for s in 1:length(sims_train)
# idxs = shuffle(1:length(sims_train))[1:100]
# idxs = shuffle(1:length(sims_train))
# idxs = 1:2
idxs = 1:length(sims_train)
for (k, s) in enumerate(idxs)
    # k = 1
    global sims, mθ, dθ

    # optm = ADAM(1000);
    # optd = ADAM(1000);
    sims = sims_train[[s]]
    # nsim = sims[1]
    info[:s] = k

    it = 400000
    α = 1000
    hist = zeros(it)
    tik = @elapsed for i = 1:it
        ∇m = gradient((p) -> loss(p, probm), mθ)
        ∇d = gradient((p) -> loss(p, probd), dθ)
        # mθ[i] = p

        # ∇m[1][[2, 3]] .= 0.0
        # ∇d[1][[1, 2]] .= 0.0

        # update!(optm, mθ, ∇m[1])
        # update!(optd, dθ, ∇d[1])

        λ = 5e-3
        # λ = 1e-4  # Better than below according to a test with 20 samples
        # λ = 1e-3  # Better than below according to a test with 20 time samples
        # λ = 1e-5  # Better than below according to a test with 20 samples
        # λ = 1e-6  # Better than below according to a test with 20 samples
        # λ = i > 1 ? 10^(floor(log10(hist[i - 1])) - 1) : 5e-3
        # λ = i > 1 ? 10^(floor(log10(hist[i - 1])) - 1) : 1e-2
        # λ = i > 1 ? 10^(floor(log10(hist[i - 1])) - 2) : 5e-3
        # λ = 9e-3 
        # λ = i > 1 ? 10^(floor(log10(hist[i - 1])) - 1) : 9e-3

        mθ += -(λ * norm(mθ) / norm(∇m[1])) * ∇m[1]
        dθ += -(λ * norm(dθ) / norm(∇d[1])) * ∇d[1]

        # dθ[2] = 0.0;
        # λ = norm(∇m[1])^0.2
        # λ = norm(∇d[1])^0.9

        ρ, rmse, θ  = lossout(dθ, probd)
        hist[i] = ρ + rmse
        info[:ρ], info[:rmse], info[:θ] = ρ, rmse, θ

        # Stop criteria
        if stopcrit!(i, hist, α)
            break
        elseif (ρ < -(1-1e-4)) & (rmse < 3e-2)
            break
        end

        # Print info it
        println(printit(i, hist[i], info))
    end

    push!(mlist, mθ[1])
    push!(dlist, [0.0, dθ[1]])
    # push!(dlist, [dθ[1], 0.0])
    push!(xlist, sims[1]["all"])
    push!(klist, sims[1]["key"])

    i = findfirst(hist .== 0.0)
    push!(losslist, hist[i-1])
    push!(timelist, tik)

    info[:eta] = (length(idxs) - k) * (sum(timelist)/length(timelist))
    info[:tik] = tik
end


mθ = [ym[1]]
nsim = sims[376]
plot(nsim["t"], nsim["y"])
plot!(nsim["t"], predict(mθ, nsim, probm), legend=false)
# plot!(nsim["t"], predict([nsim["Iᵥ"], 0.0, 0.0, 0.0], nsim, probm))
# # plot!(nsim["t"], predict([nsim["Iᵥ"]], nsim, probm))

# metrics(nsim["y"], predict(mθ, nsim, probm))

# Turb u + u^2
# rmse: 3.5593e-02 mse: 1.2669e-03 ρ: 9.98355e-01 a: 1.00254e+00 b: 5.71178e-01
# 0.0
# 1.4557456294974292e7
# 1.0167066969264485e6

# Turb u^2
# rmse: 3.9396e-03 mse: 1.5520e-05 ρ: 9.99983e-01 a: 1.00126e+00 b: 5.71178e-01
# 0.0
# 0.0
# 5.9435415142057526e10

# Lam u + u^2
# rmse: 4.2286e-03 mse: 1.7881e-05 ρ: 9.99860e-01
# 0.0
# 4.694348794776163e7
# 32840.04644376066
# Lam u
# rmse: 4.2307e-03 mse: 1.7899e-05 ρ: 9.99860e-01 a: 9.98314e-01 b: 1.27581e+00
#  0.0
#  4.694290276424844e7
#  0.0

# plot(predict(dθ, nsim, probd))
# plot!(predict([nsim["Rᵥ"]], nsim, probd))

extrema(losslist)
histogram(log10.(abs.(losslist)))
histogram(losslist)

sum(losslist .> -0.98)

# histogram(losslist)
# mean(losslist)

# Save parameter
filename2 = "./data/fit_md_sim_doe4_all3.jld2"
file2 = jldopen(filename2, "r")

file2["key"]
file2["x"]
file2["m"]
file2["d"]
file2["loss"]
close(file2)

# Save parameter
filename = "./data/fit_md_sim_doe4_all.jld2"
file = jldopen(filename, "w")
# file = jldopen(filename, "a+")

file["key"] = vcat(klist, file2["key"])
file["x"] = hcat(hcat(xlist...), file2["x"])
file["m"] = vcat(Float64.(mlist), file2["m"])
file["d"] = hcat(hcat(dlist...), file2["d"])
file["loss"] = vcat(Float64.(losslist), file2["loss"])
close(file)

# Save parameter
filename = "./data/fit_md_sim_doe4_all2.jld2"
file = jldopen(filename, "w")
# file = jldopen(filename, "a+")

file["key"] = klist
file["x"] = hcat(xlist...)
file["m"] = Float64.(mlist)
file["d"] = hcat(dlist...)
file["loss"] = Float64.(losslist)
close(file)


M = Float64.(mlist)
D = Float64.(dlist)

iv = [s["Iᵥ"] for s in sims_all[idxs]]

plot(abs.(M - iv)./iv, yaxis=:log)
plot(abs.(M - iv)./iv)

n = argmax(abs.(M - iv)./iv)
sims_all[n]["Iᵥ"]
M[n]

mθ = M[[n]]
dθ = D[[n]]
losslist[n]

sims = sims_all[[idxs[n]]]
klist[n]
sims[1]["k"]

plot(sims[1]["t"], predict(mθ, sims[1], probm))
plot!(sims[1]["t"], sims[1]["y"])


# Testes
# Trying CGD
i = 1
∇mx = gradient((p) -> idxloss(p, loss, i, mθ, probm), mθ[i])
f̂m(α) = idxloss(mθ[i] - α*∇mx[1], loss, i, mθ, probm)
α = ArmijoLineSearch(f̂m, ∇mx[1], norm(∇mx[1])^-0.1)
mθ[i] = mθ[i] - α*∇mx[1]

for k in 1:100
    for i in [1,2,3]
        ∇mx = gradient((p) -> idxloss(p, loss, i, mθ, probm), mθ[i])
        f̂m(α) = idxloss(mθ[i] - α*∇mx[1], loss, i, mθ, probm)
        α = ArmijoLineSearch(f̂m, ∇mx[1], norm(∇mx[1])^-0.2)
        mθ[i] = mθ[i] - α*∇mx[1]
    end
    for i in [1,2,3]
        ∇dx = gradient((p) -> idxloss(p, loss, i, dθ, probd), dθ[i])
        f̂d(α) = idxloss(dθ[i] - α*∇dx[1], loss, i, dθ, probd)
        α = ArmijoLineSearch(f̂d, ∇dx[1], norm(∇dx[1])^-0.2)
        dθ[i] = dθ[i] - α*∇dx[1]
    end
    println(k)
end

#FAST

# -----------------------------------------------------------------------------
# Approx functions

function d̂!(u, p)
    return (nd((u .- x̄) ./ x̃).*ỹd.+ȳd)[1]
end;

function m̂!(u, p)
    return (nm((u .- x̄) ./ x̃).*ỹm.+ȳm)[1]
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
    tmp_prob = remake(prob; u0 = nsim["u0"], p = p)
    return vec(solve(tmp_prob, Tsit5(); saveat = nsim["trng"])) ./ nsim["A"]
end



predict(123, sims[123], probm)nsim["trng"]

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
x[:, j] - u

m̂!(u, 1.0)
nsim["Iᵥ"]
yd[k]
d̂!(u, 1.0)
nsim["Rᵥ"]
nsim["Re"]

plotsims(sims[[k]], [k], probm)
