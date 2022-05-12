
include("lib_dev_fit.jl")
using Serialization
using Random
using CUDA

# Get Reynolds keys

# -----------------------------------------------------------------------------
# Get data
# filefit = jldopen("./data/fit_md_sim_doe3.jld2", "r")
file = jldopen("../numerical/data/sim_doe4.jld2", "r")
filefit = jldopen("./data/fit_md_sim_doe4_lam_1000_2000.jld2", "r")
filefit = jldopen("./data/fit_md_sim_doe4_all.jld2", "r")
# close(filefit)

Re = getre(file, 0, 2200)
Re = getre(file, 3000, 10000)
lam = []
for k in keys(Re)
    result = findfirst(filefit["key"] .== k)
    if !isnothing(result)
        push!(lam, result)
    end
end

x = Float32.(filefit["x"])[:, lam]
# x = Float32.(filefit["x"])[[1, 2, 3, 4, 7], :]
# x = Float64.(filefit["x"])
yd = Float32.(filefit["d"])[2, lam]
# yd = Float32.(filefit["d"])[1, lam]
# yd = Float32.(filefit["d"])[1, :]
# yd = Float64.(filefit["d"])
ym = Float32.(filefit["m"])[lam]
# ym = Float64.(filefit["m"])
j = filefit["loss"][lam]
kl = filefit["key"][lam]
# j = filefit["loss"]
id = collect(1:size(x, 2))

nσ = 1
yv = yd[:][yd[:] .!= 0]
yv = yd
rm_idx = collect(1:length(yv))[yv.>(mean(yv)+nσ*std(yv))]
rm_idx = vcat(rm_idx, id[ym.>(mean(ym)+nσ*std(ym))])
rm_idx = vcat(rm_idx, id[j .> -0.98])
rm_idx = Set(rm_idx)

kp_idx = filter(x -> !(x in rm_idx), 1:length(yv))

yd = yd[kp_idx]
# yd = yd[:, kp_idx]
ym = ym[kp_idx]
x = x[:, kp_idx]
id = collect(1:length(id[kp_idx]))
kl = kl[kp_idx]
kl

x[2, :] = log10.(x[2, :])
x[3, :] = log10.(x[3, :])
x[5, :] = log10.(x[5, :])
x[7, :] = log10.(x[7, :])

# -----------------------------------------------------------------------------
# Data
x̄, x̃ = mean(x; dims = 2), std(x; dims = 2)
Xs = (x .- x̄) ./ x̃

ȳd, ỹd = mean(yd), std(yd)
ȳm, ỹm = mean(ym), std(ym)

Ysd = (yd .- ȳd) ./ ỹd
Ysm = (ym .- ȳm) ./ ỹm

# datad = [(i, j) for (i, j) in zip(eachcol(xs), ysd)]
# datam = [(i, j) for (i, j) in zip(eachcol(xs), ysm)]

γ = 0.8
idxs = shuffle(id)
ntrn = floor(Int, γ * length(id))
trnidx = idxs[1:ntrn]
tstidx = idxs[ntrn+1:end]

tstidx = id[[all(c) for c in eachcol(-1.5.< Xs .< 1.5)]]
tstidx = id[[all(c) for c in eachcol(-1.2 .< Xs .< 1.2)]]
if length(tstidx) > (1 - γ)*length(id)
    ntst = floor(Int, (1 - γ) * length(id))
    tstidx = shuffle(tstidx)[1:ntst]
else
    trnidx = collect(setdiff(Set(id), Set(tstidx)))
end

xs = gpu(Xs[:, trnidx])
# ysd = gpu(Ysd[:, trnidx])
ysd = gpu(Ysd[trnidx])
ysm = gpu(Ysm[trnidx])
kl[id[trnidx]]
xs[:, 1]

xt = gpu(Xs[:, tstidx])
ytd = gpu(Ysd[tstidx])
# ytd = gpu(Ysd[:, tstidx])
ytm = gpu(Ysm[tstidx])

# -----------------------------------------------------------------------------
# Damper
# Train
opt = RMSProp(0.0001)
opt = ADAM(0.01)

scatter(xs[1, :], xs[2, :])
scatter!(xt[1, :], xt[2, :])

scatter(xs[1, :], xs[3, :])
scatter!(xt[1, :], xt[3, :])

scatter(xs[1, :], xs[4, :])
scatter!(xt[1, :], xt[4, :])

scatter(xs[1, :], xs[5, :])
scatter!(xt[1, :], xt[5, :])

scatter(xs[1, :], xs[6, :])
scatter!(xt[1, :], xt[6, :])

scatter(xs[1, :], xs[7, :])
scatter!(xt[1, :], xt[7, :])


# nn = gpu(Chain(Dense(size(xs, 1) => 256, relu), Dense(256 => 1)))
# nn = gpu(Chain(Dense(size(xs, 1) => 32, tanh), Dense(32 => 32, tanh), Dense(32 => 1)))
# UHUUUUU - Laminar resistance
nn = gpu(Chain(Dense(size(xs, 1) => 14, tanh), Dense(14 => 1)))
nn = gpu(Chain(Dense(size(xs, 1) => 8, tanh), Dense(8 => 1)))
# UHUUUUU - Laminar mass
# nn = gpu(Chain(Dense(size(xs, 1) => 8, tanh), Dense(8 => 1)))
# nn = gpu(Chain(Dense(size(xs, 1) => 16, relu), Dense(16 => 1)))
# nn = gpu(Chain(Dense(size(xs, 1) => 256, relu), Dense(256 => 256, relu),Dense(256 => 1)))
#nn = gpu(Chain(Dense(size(xs, 1) => 10, relu), Dense(10 => 10, relu), Dense(10 => 1)))
# nn = gpu(
#     Chain(
#         Dense(size(xs, 1) => 16, relu),
#         Dense(16 => 16, relu),
#         Dense(16 => 16, relu),
#         Dense(16 => 1),
#     ),
# )

ps = Flux.params(nn)

# loss(x, y) = Flux.Losses.mse(nn(x), y)
loss(x, y) = sqrt(Flux.Losses.mse(nn(x), y))
loss(xs, ysd')

# loss(xs, ysd)
info = gpu(Dict(:c => 0.0))
it = 900000
hist = gpu(zeros(it, 2))
for i = 1:it
    # batch = gpu(shuffle(1:size(xs, 2))[1:30])
    batch = 1:size(xs, 2)
    # ∇p = gradient(() -> loss(xs[:, batch], ysd[:, batch]), ps)
    ∇p = gradient(() -> loss(xs[:, batch], ysd[batch]'), ps)
    # ∇p = gradient(() -> loss(xs[:, batch], ysm[batch]'), ps)

    update!(opt, ps, ∇p)

    # λ = i > 1 ? 10^(floor(log10(hist[i - 1])) - 2) : λ

    # for (p, g) in zip(ps, ∇p)
    #     p .+= -(λ*norm(p)/norm(g))*g
    # end

    # hist[i, 1] = loss(xs, ysm')
    # hist[i, 2] = loss(xt, ytm')
    hist[i, 1] = loss(xs, ysd')
    hist[i, 2] = loss(xt, ytd')
    # hist[i, 1] = loss(xs, ysd)
    # hist[i, 2] = loss(xt, ytd)
    info[:c] = hist[i, 2]

    # hist[i] = loss(xs, ysm')
    # info[:c] = loss(xt, ytm')

    # Stop criteria
    # if stopcrit!(c, i, hist, α) break end
    # oi = stopcrit!(c, i, hist, α)
    # info[:c] = c

    # Print info it
    println(printit(i, hist[i], info))
end

ŷ = vec(nn(xt));
# metrics(ŷ, vec(ytd))
metrics(ŷ, ytm)

ŷ = vec(nn(xs));
metrics(ŷ, vec(ysd))
metrics(ŷ, ysm)

# Plot history
idx = findfirst(hist[:, 1] .== 0.0)
i = isnothing(idx) ? size(hist, 1) + 1 : idx
plot(hist[1:(i-1), 1]; yaxis = :log10)
plot!(hist[1:(i-2), 2]; yaxis = :log10)


serialize("./data/model_sim_rest_doe4.dat", nn)
serialize("./data/model_sim_inet_doe4.dat", nn)
serialize("./data/model_sim_varst_doe4.dat", Dict(:x_mean => x̄, :x_std => x̃, :yd_mean => ȳd, :yd_std => ỹd, :ym_mean => ȳm, :ym_std => ỹm))


# -----------------------------------------------------------------------------
# Mass
nn = Chain(Dense(size(xs, 1) => 28, relu), Dense(28 => 1))
ps = Flux.params(nn)

loss(x, y) = Flux.Losses.mse(nn(x), y)
loss(x, y) = sqrt(Flux.Losses.mse(nn(x), y))
loss(xs, ysd')
loss(xs, ysm')

# Train
# opt = ADAM(0.01)
# opt = RMSProp(0.00001)

info = Dict(:c => 0.0)

λ = 5e-2vcat(Float64.(losslist), file2["loss"])
hist = zeros(it)
for i = 1:it
    ∇p = gradient(() -> loss(xs, ysd'), ps)

    λ = i > 1 ? 10^(floor(log10(hist[i-1])) - 2) : λ

    for (p, g) in zip(ps, ∇p)
        p .+= -(λ * norm(p) / norm(g)) * g
    end

    hist[i] = loss(xs, ysd')

    # Stop criteria
    # if stopcrit!(c, i, hist, α) break end
    oi = stopcrit!(c, i, hist, α)
    info[:c] = c

    # Print info it
    println(printit(i, hist[i], info))
end

ŷ = vec(nn(xt));
metrics(ŷ, ytd)

i = findfirst(hist .== 0.0)
plot(log10.(hist[1:(i-1)]))

ŷ2 = ŷ .* std(y) .+ mean(y)
metrics(ŷ2, y)

@save "./data/model_sim_doe3.bson" nn
