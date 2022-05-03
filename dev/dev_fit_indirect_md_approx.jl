include("lib_dev_fit.jl")
using BSON: @load

# -----------------------------------------------------------------------------
# Get data
filefit = jldopen("./data/fit_md_sim_doe3.jld2", "r")

x = filefit["x"]
ym = filefit["d"]
yd = filefit["m"]

# -----------------------------------------------------------------------------
# fit
# Data
xs = (x .- mean(x, dims=2))./std(x, dims=2)
ysd = (yd .- mean(yd))./std(yd)
ysm = (ym .- mean(ym))./std(ym)

datad = [(i, j) for (i, j) in zip(eachcol(xs), ysm)]
datam = [(i, j) for (i, j) in zip(eachcol(xs), ysm)]

# NN
nn = Chain(Dense(size(xs, 1) => 28, relu), Dense(28 => 1))
ps = Flux.params(nn)

loss(x, y) = Flux.Losses.mse(nn(x), y)
loss(x, y) = sqrt(Flux.Losses.mse(nn(x), y))
loss(xs, ysd')
loss(xs, ysm')

# Train
opt = ADAM(0.01)
opt = RMSProp(0.00001)

it = 10000;
λ = 5e-2;
α = 50;
it = 1000000;
hist = zeros(it);
for i in 1:ity(() -> loss(xs, ys'), ps)

    λ = i > 1 ? 10^(floor(log10(hist[i - 1])) - 2) : λ
    for (p, g) in zip(ps, ∇p)
        p .+= -(λ*norm(p)/norm(g))*g
    end

    hist[i] = loss(xs, ys')
    @printf("It: %d - loss: %.9e\n", i, hist[i])
end

ŷ = vec(nn(xs))
metrics(ŷ, ys)

i = findfirst(hist .== 0.0)
plot(log10.(hist[1:i-1]))

ŷ2 = ŷ .* std(y) .+ mean(y)
metrics(ŷ2, y)

@save "./data/model_sim_doe3.bson" nn

