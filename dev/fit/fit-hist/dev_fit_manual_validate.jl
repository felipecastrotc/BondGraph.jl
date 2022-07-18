include("lib_dev_fit.jl")

# -----------------------------------------------------------------------------
# Get data

# file = jldopen("../numerical/data/sim_doe1.jld2")
file = jldopen("../numerical/data/sim_doe3.jld2", "r")

sims_tmp = gendata(file, 10, 1000, 2000);
sims_tmp = gendata(file, 10, 1800, 4000);
sims = sims_tmp[[1, 2]]

plotsims(sims_tmp)
close(file)

# ------------------------------------------------------------------------------
# Manual

# ------------------------------------------------------------------------------
# Simple manual equation
function diffequ!(du, u, p, t)
    du[1] = (p[1] - u[1]*p[2])/(p[3])
end

# Set manual problem
u0 = [0.0]
probu = ODEProblem(diffequ!, u0, (0.0, 2.0), zeros(3))

ns = 2
Iᵥ = sims[ns]["Iᵥ"]
Rᵥ = sims[ns]["Rᵥ"]
Δp = sims[ns]["Δp"]
A = sims[ns]["A"]
ŷ = vec(solve(probu, Tsit5(), p=[Δp, Rᵥ, Iᵥ], saveat=sims[ns]["trng"]))./A

# Plot comparison
plot(sims[ns]["t"], sims[ns]["y"], label=".")
plot!(sims[ns]["t"], ŷ, label="~")
# Error
sqrt(mean((ŷ - sims[ns]["y"]).^2))

# -----------------------------------------------------------------------------
# Simple manual function test

# -----------------------------------------------------------------------------
# Dummy BG approximation
m̂!(u, p) = @. p[1]*u[1] * u[2] / (p[2]*u[3])
d̂!(u, p) = @. p[1]*u[1]*u[2]/(p[2]*u[3]^2)

# -----------------------------------------------------------------------------
# Simple manual m and d functions

function diffeq!(du, u, p, t)
    du[1] = (p[1] - u[1]*d̂!(p[2], p[3]))/m̂!(p[4], p[5])
end

# Set manual problem
du = [0.0]
u0 = [0.0]
probf = ODEProblem(diffeq!, u0, (0.0, 20.0))

# Set the parameters
ns = 2
mθ = [1.0, 1.0]
dθ = [128.0, π]
dp = sims[ns]["dp"]
mp = sims[ns]["mp"]
Δp = sims[ns]["Δp"]
A = sims[ns]["A"]
p = [Δp, dp, dθ, mp, mθ]

# Solve
ŷ = vec(solve(probf, p=p, saveat=sims[ns]["trng"]))./A

# Plot comparison
plot(sims[ns]["t"], sims[ns]["y"], label=".")
plot!(sims[ns]["t"], ŷ, label="~")
# Error
sqrt(mean((ŷ - sims[ns]["y"]).^2))

# ------------------------------------------------------------------------------
# Fitting

#-----------------------------------------------------------------------------
# Laurent 

# Laurent Layer
λ(X, W, n) = sum((X' * W[:, i]) .^ c for (i, c) in enumerate(n))
# λ(X, W) = (X' * W[:, 1]) .^ -2 + (X' * W[:, 2]) .^ -1 + (X .* 0)' * W[:, 3] + (X' * W[:, 4]) .^ 1 + (X' * W[:, 5]) .^ 2

# Layer settings
pol = -2:2
pol = 1:3
# pol = -4:4
# nvarm = 4
# nvard = 5
nvarm = 3
nvard = 3

mθ = rand(nvarm, length(pol))[:]
dθ = rand(nvard, length(pol))[:]

# Create a dummy function
m̂!(u, p) = λ(u, reshape(p, nvarm, length(pol)), pol)
d̂!(u, p) = λ(u, reshape(p, nvard, length(pol)), pol)
# m̂!(u, p) = λ(u, reshape(p, nvarm, 5))
# d̂!(u, p) = λ(u, reshape(p, nvard, 5))

#-----------------------------------------------------------------------------
# Find BG

# data["Rᵥ"] = 128 * μ * L / (π * d^4)
# data["Iᵥ"] = ρ * L / A

p = sims[1]["dp"]
p[0]*p[1]

p'*rand(3)

sims[1]["A"]
using Symbolics

@variables a b c


mθ = ones(3, length(pol))[:]
#     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
mθ = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# mθ = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
# mθ = [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]
m̂!([a, b, c], mθ)

expand(m̂!([a, b, 1/c], mθ))

([a, b, c]' * [0, 0, 0 ]).^0

# Laurent Layer
λ(X, W, n) = sum((X' * W[:, i]) .^ c for (i, c) in enumerate(n))
# λ(X, W) = (X' * W[:, 1]) .^ -2 + (X' * W[:, 2]) .^ -1 + (X .* 0)' * W[:, 3] + (X' * W[:, 4]) .^ 1 + (X' * W[:, 5]) .^ 2

# Layer settings
# pol = -2:2
# pol = [-2, -1, 1, 2]
pol = 1:4
# nvarm = 4
# nvard = 5
nvarm = 3
nvard = 3

mθ = rand(nvarm, length(pol))[:]
dθ = rand(nvard, length(pol))[:]

# Create a dummy function
m̂!(u, p) = λ(u, reshape(p, nvarm, length(pol)), pol)
d̂!(u, p) = λ(u, reshape(p, nvard, length(pol)), pol)
# m̂!(u, p) = λ(u, reshape(p, nvarm, 5))
# d̂!(u, p) = λ(u, reshape(p, nvard, 5))

# -----------------------------------------------------------------------------
# m̂ = FastChain(FastDense(4, 5, relu), FastDense(5, 1))
# d̂ = FastChain(FastDense(5, 5, relu), FastDense(5, 1))
m̂ = FastChain(FastDense(3, 10, relu), FastDense(10, 1))
d̂ = FastChain(FastDense(3, 10, relu), FastDense(10, 1))
# # m̂!(u, p) = dot(m̂(u, p), [1])
# # d̂!(u, p) = dot(d̂(u, p), [1])
m̂!(u, p) = abs(dot(m̂(u, p), [1]))*1000
d̂!(u, p) = abs(dot(d̂(u, p), [1]))*1000

mθ = initial_params(m̂)
dθ = initial_params(d̂)


# -----------------------------------------------------------------------------
# Train setup
function diffeq!(du, u, p, t)
    du[1] = (p[1] - u[1]*d̂!(p[2], p[3]))/m̂!(p[4], p[5])
end

data = sims[1]
mp = data["mp"]
dp = data["dp"]
Δp = data["Δp"]
u0 = sims[1]["u0"]

diffeqd!(du, u, p, t) = diffeq!(du, u, [Δp, dp, p, mp, mθ], t)
diffeqm!(du, u, p, t) = diffeq!(du, u, [Δp, dp, dθ, mp, p], t)

probm = ODEProblem(diffeqm!, u0, (0.0, 20.0), mθ)
probd = ODEProblem(diffeqd!, u0, (0.0, 20.0), dθ)

function predict(p, data, prob)
    global mp, dp, Δp
    tmp_prob = remake(prob, u0=data["u0"], p=p)
    mp = data["mp"]
    dp = data["dp"]
    Δp = data["Δp"]
    # diffeqd!(du, u, p, t) = diffeq!(du, u, [data["dp"], p, data["mp"], mθ, data["Δp"]], t)
    # diffeqm!(du, u, p, t) = diffeq!(du, u, [data["dp"], dθ, data["mp"], p, data["Δp"]], t)
    vec(solve(tmp_prob, Tsit5(), saveat=data["trng"]))./data["A"]
end

# -----------------------------------------------------------------------------
# Define loss functions

function loss(p, prob)
    # return sum((s["y"] .- predict_neuralode(p, s)) .^ 2)/n
    return sqrt(mean([mean((s["y"] .- predict(p, s, prob)) .^ 2) for s in sims]))
end

# -----------------------------------------------------------------------------
# Generate gradient functions
∇d = gradient((p) -> loss(p, probd), dθ)
# ∇m = gradient((p) -> loss(p, probm), mθ)

# -----------------------------------------------------------------------------
# Test the loss and predict functions
# Change the m and d functions to the BG approximation  and evaluate the results
# Generate random weights
# mθ = [1.0, 1.0]
# dθ = [128.0, π]

# plot()
# for (i, s) in enumerate(sims)
#     c = palette(:default)[i%16]
#     plot!(s["t"], s["y"], ls=:dash, linecolor=c)
#     plot!(s["t"], predict(dθ, s, probd), linecolor=c)
# end
# plot!()

# ns = 1
# mp = sims[ns]["mp"]
# dp = sims[ns]["dp"]

# Iᵥ = sims[ns]["Iᵥ"]
# Rᵥ = sims[ns]["Rᵥ"]
# m̂!(mp, mθ)
# d̂!(dp, dθ)

# -----------------------------------------------------------------------------
# Generate random weights
# mθ = rand(nvarm, length(pol))[:]
# dθ = rand(nvard, length(pol))[:]

# mθ = initial_params(m̂)
# dθ = initial_params(d̂)

# -----------------------------------------------------------------------------
# Plot initial
ns = 2
data = sims[ns]
plot(data["t"], data["y"])
plot!(data["t"], predict(dθ, sims[ns], probd))
# plot(data["t"], predict(dθ, sims[1], probd))

loss(dθ, probd)

# -----------------------------------------------------------------------------
# Training
optm = ADAM(1);
optd = ADAM(0.1);

# optd = AdaBelief(1);
# optm = AdaBelief(1);

# sims2 = gendata(file, 10, 1000, 6000)
sims = sims2
# sims = sims2[[2, 3]]
plotsims(sims)

mθ = [1.0, 1.0]
dθ = [128.0, π]

it = 24000;
hist = zeros(it);
for i = 1:it
    # i = 2
    ∇d = gradient((p) -> loss(p, probd), dθ)
    ∇m = gradient((p) -> loss(p, probm), mθ)

    update!(optm, mθ, ∇m[1])
    update!(optd, dθ, ∇d[1])
    # dθ += -(5e-2*norm(dθ)/norm(∇d[1]))*∇d[1]
    # dθ += -(5e-5*norm(dθ)/norm(∇d[1]))*∇d[1]
    # dθ += -(1e-6*norm(dθ)/norm(∇d[1]))*∇d[1]
    # dθ += -(1e-10*norm(dθ)/norm(∇d[1]))*∇d[1]
    # mθ += -(1e-1*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-4*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-6*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-7*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-8*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-2*norm(mθ)/norm(∇m[1]))*∇m[1]

    hist[i] = loss(dθ, probd)
    ns = 2
    mp = sims[ns]["mp"]
    dp = sims[ns]["dp"]
    Iᵥ = sims[ns]["Iᵥ"]
    Rᵥ = sims[ns]["Rᵥ"]
    @printf("It: %d - loss: %.9e m: %.5e/%.5e d: %.5e/%.5e\n", i, hist[i], m̂!(mp, mθ), Iᵥ,d̂!(dp, dθ), Rᵥ)
    print

    if i == 1
        p = plot(sims[ns]["t"], sims[ns]["y"])
    elseif ((i % 10) == 0)
        p = plot!(sims[ns]["t"], predict(dθ, sims[ns], probd))
        # p = plot()
        # for (i, s) in enumerate(sims)
        #     c = palette(:default)[i%16]
        #     plot!(s["t"], s["y"], ls=:dash, linecolor=c)sum
        #     p = plot!(s["t"], predict(dθ, s, probd), linecolor=c)
        # end
        display(p)
        
    end

    if sum(diff(hist[max(i - 6, 1):i]) .> 0) > 100
        break
    end

end

plotsims(sims, dθ)

dθ
mθ

i = findfirst(hist .== 0.0)
# i = length(hist)
plot(hist[1:i-1], yaxis=:log10)
# plot(hist, yaxis=:log10)


# -----------------------------------------------------------------------------
# Startup NN with BG

sims = gendata(file, 10, 0, 1e10)

# Data
x = Float32.(vcat([s["dp"] for s in sims]'...))'
y = Float32.(vcat([s["Rᵥ"] for s in sims]'...))

xs = (x .- mean(x, dims=2))./std(x, dims=2)
ys = (y .- mean(y))./std(y)
data = [(i, j) for (i, j) in zip(eachcol(xs), ys)]

# NN
nd = Chain(Dense(3 => 9, relu), Dense(9 => 1))
ps = Flux.params(nd)

lossd(x, y) = Flux.Losses.mse(nd(x), y)
lossd(xs, ys')

# Train
optd = ADAM(0.000001)

it = 10000;
hist = zeros(it);
for i in 1:it
    Flux.train!(lossd, ps, data, optd)
    hist[i] = lossd(xs, ys')
    @printf("It: %d - loss: %.9e\n", i, hist[i])
end

scatter(ys, nd(xs)')

# Tranform to the pattern being used
dθ, d̂ = Flux.destructure(nd);

x̄ = mean(x, dims=2)
x̃ = std(x, dims=2)
ȳ = mean(y)
ỹ = std(y)

d̂!(x, p) = (d̂(p)((x .- x̄)./x̃)[1])*ỹ + ȳ

# -----------------------------------------------------------------------------
# Startup NN with BG

sims = gendata(file, 10, 0, 1e10)
sims = gendata(file, 10, 0, 2400)

# Data
x = Float32.(vcat([s["dp"] for s in sims]'...))'
y = Float32.(vcat([s["Rᵥ"] for s in sims]'...))

xs = (x .- mean(x, dims=2))./std(x, dims=2)
ys = (y .- mean(y))./std(y)
data = [(i, j) for (i, j) in zip(eachcol(xs), ys)]

# NN
hs = 18
nd = Chain(Dense(3 => hs, relu), Dense(hs => hs, relu), Dense(hs => hs, relu), Dense(hs => 1))
ps = Flux.params(nd)

lossd(x, y) = sqrt(Flux.Losses.mse(nd(x), y))
lossd(xs, ys')

# Train
optd = ADAM(0.0001)

it = 3000;
hist = zeros(it);
for i in 1:it
    Flux.train!(lossd, ps, data, optd)
    hist[i] = lossd(xs, ys')
    @printf("It: %d - loss: %.9e\n", i, hist[i])
end

scatter(ys, nd(xs)')

# Tranform to the pattern being used
dθ, d̂ = Flux.destructure(nd);

x̄ = mean(x, dims=2)
x̃ = std(x, dims=2)
ȳ = mean(y)
ỹ = std(y)

d̂!(x, p) = (d̂(p)((x .- x̄)./x̃)[1])*ỹ + ȳ

# -----------------------------------------------------------------------------
# Startup NN with BG

sims = gendata(file, 10, 0, 2400)

# Data
x = Float32.(vcat([s["dp"] for s in sims]'...))'
x[3, :] = 1/x[3, :]
y = Float32.(vcat([s["Rᵥ"] for s in sims]'...))
xs = (x .- mean(x, dims=2))./std(x, dims=2)
ys = (y .- mean(y))./std(y)

lossl(x, y) = sqrt(mean((x - y).^2))

lossl(y, m̂!(x, mθ))

# Train
mθ = randn(nvarm, length(pol))[:]
optl = ADAM(1)

it = 100000;
hist = zeros(it);
for i in 1:it
    ∇l = gradient((p) -> lossl(m̂!(xs, p), ys), mθ)
    
    update!(optl, mθ, ∇l[1])
    # mθ += -(1e-5*norm(mθ)/norm(∇l[1]))*∇l[1]

    hist[i] = lossl(m̂!(x, mθ), y)
    @printf("It: %d - loss: %.9e\n", i, hist[i])
end

scatter(y, m̂!(x, mθ))

# Tranform to the pattern being used
dθ, d̂ = Flux.destructure(nd);

x̄ = mean(x, dims=2)
x̃ = std(x, dims=2)
ȳ = mean(y)
ỹ = std(y)

d̂!(x, p) = (d̂(p)((x .- x̄)./x̃)[1])*ỹ + ȳ

