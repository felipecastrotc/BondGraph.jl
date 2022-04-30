using Plots, JLD2, LinearAlgebra, Statistics
using Flux
using Flux.Optimise: ADAM, update!
using Printf
using BSON: @save

# -----------------------------------------------------------------------------
# Support

# Set loss and model functions
lossrmse(x, y) = sqrt(mean((x - y).^2))

# Metrics
function metrics(y, ŷ)
    rmse = lossrmse(y, ŷ)
    a, b = y\hcat(ŷ, ones(length(y)))
    ρ = cor(y, ŷ)
    
    err = (1 .- abs.((y - ŷ)./y)).*100

    @printf("rmse: %.4e mse: %.4e ρ: %.5e a: %.5e b: %.5e ē: %.5e em: %.5e\n", rmse, rmse^2, ρ, a, b, mean(err), median(err))

    scatter(y, ŷ)
end

# -----------------------------------------------------------------------------
# Get data
filef = jldopen("./data/fit_md_sim_doe3.jld2", "r")

x = filef["x"]
y = filef["d"]
y = filef["m"]

# -----------------------------------------------------------------------------
# fit
# Data
xs = (x .- mean(x, dims=2))./std(x, dims=2)
ys = (y .- mean(y))./std(y)
data = [(i, j) for (i, j) in zip(eachcol(xs), ys)]

# NN
nn = Chain(Dense(size(xs, 1) => 28, relu), Dense(28 => 1))
ps = Flux.params(nn)

loss(x, y) = Flux.Losses.mse(nn(x), y)
loss(x, y) = sqrt(Flux.Losses.mse(nn(x), y))
loss(xs, ys')

# Train
opt = ADAM(0.01)
opt = RMSProp(0.00001)
opt = ADAM(0.000001)

loss(xs, ys')

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

