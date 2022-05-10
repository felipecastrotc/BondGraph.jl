
include("lib_dev_fit.jl")
using BSON: @load
using Random
using CUDA

# Get Reynolds keys


# -----------------------------------------------------------------------------
# Get data
# filefit = jldopen("./data/fit_md_sim_doe3.jld2", "r")
filefit = jldopen("./data/fit_md_sim_doe4_lam_1000_2000.jld2", "r")
filefit = jldopen("./data/fit_md_sim_doe4_all.jld2", "r")
Re = getre(file, 500, 2400)
lam = []
for k in keys(Re)
    result = findfirst(filefit["key"] .== k)
    if !isnothing(result)
        push!(lam, result)
    end
end

x = Float32.(filefit["x"])[:, lam]
#x = Float32.(filefit["x"])[[1, 2, 3, 4, 7], :]
# x = Float64.(filefit["x"])
yd = Float32.(filefit["d"])[lam]
# yd = Float64.(filefit["d"])
ym = Float32.(filefit["m"])[lam]
# ym = Float64.(filefit["m"])
id = 1:length(yd)

nσ = 1
rm_idx = collect(1:length(yd))[yd.>(mean(yd)+nσ*std(yd))]
rm_idx = vcat(rm_idx, collect(1:length(ym))[ym.>(mean(ym)+nσ*std(ym))])
kp_idx = filter(x -> !(x in rm_idx), 1:length(yd))

yd = yd[kp_idx]
ym = ym[kp_idx]
x = x[:, kp_idx]
id = id[kp_idx]

x[2, :] = log10.(x[2, :])
x[3, :] = log10.(x[3, :])
x[5, :] = log10.(x[5, :])

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

idxs = shuffle(1:size(Xs, 2))
ntrn = floor(Int, 0.9 * length(ym))
id_train = id[idxs[1:ntrn]]
xs = gpu(Xs[:, idxs[1:ntrn]])
ysd = gpu(Ysd[idxs[1:ntrn]])
ysm = gpu(Ysm[idxs[1:ntrn]])

id_test = id[idxs[ntrn:end]]
xt = gpu(Xs[:, idxs[ntrn:end]])
ytd = gpu(Ysd[idxs[ntrn:end]])
ytm = gpu(Ysm[idxs[ntrn:end]])

# -----------------------------------------------------------------------------
# Damper
# Train
opt = RMSProp(0.0001)
opt = ADAM(0.01)

nn = gpu(Chain(Dense(size(xs, 1) => 256, relu), Dense(256 => 1)))
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

info = gpu(Dict(:c => 0.0))

λ = 1e-1
λ = 5e-2
α = 10000
c = 0
it = 600000
hist = gpu(zeros(it, 2))
for i = 1:it
    # batch = gpu(shuffle(1:size(xs, 2))[1:30])
    batch = 1:size(xs, 2)
    # ∇p = gradient(() -> loss(xs[:, batch], ysd[batch]'), ps)
    ∇p = gradient(() -> loss(xs[:, batch], ysm[batch]'), ps)

    update!(opt, ps, ∇p)

    # λ = i > 1 ? 10^(floor(log10(hist[i - 1])) - 2) : λ

    # for (p, g) in zip(ps, ∇p)
    #     p .+= -(λ*norm(p)/norm(g))*g
    # end

    hist[i, 1] = loss(xs, ysm')
    hist[i, 2] = loss(xt, ytm')
    # hist[i, 1] = loss(xs, ysd')
    # hist[i, 2] = loss(xt, ytd')
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
metrics(ŷ, ytd)
metrics(ŷ, ytm)

ŷ = vec(nn(xs));
metrics(ŷ, ysd)
metrics(ŷ, ysm)

# Plot history
idx = findfirst(hist[:, 1] .== 0.0)
i = isnothing(idx) ? size(hist, 1) + 1 : idx
plot(hist[1:(i-1), 1]; yaxis = :log10)
plot!(hist[1:(i-1), 2]; yaxis = :log10)

ŷ2 = ŷ .* std(y) .+ mean(y)
metrics(ŷ2, y)

nd = deepcopy(nn)
@save "./data/model_sim_res_doe4.bson" nd
nm = deepcopy(nn)
@save "./data/model_sim_ine_doe4.bson" nm

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

λ = 5e-2
α = 50
c = 0
it = 100000
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
