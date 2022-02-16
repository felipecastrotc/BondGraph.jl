using Flux, Statistics, Plots

# x = rand(-20:0.00001:20, 1000);
x = rand(-40:0.00001:40, 100);
y = 0.1 * x .^ 3 + 3 * x .^ 2 + x;
y = 0.1 * x .^ 3;

nn = 10
σ = tanh
σ = relu
m = Chain(Dense(1, nn, σ), Dense(nn, 1))
m = Chain(Dense(1, nn, σ), Dense(nn, nn, σ), Dense(nn, nn, σ), Dense(nn, 1))
m = Chain(Dense(1, nn, σ), Dense(nn, nn, σ), Dense(nn, 1))
loss(x, y) = Flux.Losses.mse(m(x), y)

ps = Flux.params(m)
data = zip(eachrow(x), eachrow(y))

# Flux.@epochs 1000 Flux.train!(loss, ps, data, ADAM(0.05))
Flux.@epochs 600 Flux.train!(loss, ps, data, ADAM(0.05))

x̄ = rand(-60:0.00001:60, 100);
ȳ = 0.1 * x̄ .^ 3 + 3 * x̄ .^ 2 + x̄;

scatter(x̄, m(x̄')')
scatter!(x̄, ȳ)
