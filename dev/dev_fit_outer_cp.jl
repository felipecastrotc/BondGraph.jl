using BSON: @load
using LsqFit
using SymPy

include("lib_dev_fit.jl")

# -----------------------------------------------------------------------------
# Symbolic functions

function clean_model(mdl)
    tmp = Dict()
    for i in mdl.atoms()
        if !(typeof(N(i)) <: Sym)
            tmp[i] = round(i; digits=6)
        end
    end
    return mdl.xreplace(tmp)
end

function get_sym(mdl, l, Xₛ)
    # Set LstFit model
    model(X̂, p̂) = outer(X̂, p̂, l)
    # Expression
    expr = expand(model(Xₛ, mdl.param)[1])
    # return clean_model(expr);
    return expr
end

function get_sym_cp(mdl, l, Xₛ)
    order, rnk = l
    # Set LstFit model
    model(X̂, p̂) = CP(X̂, p̂, rnk, order)
    # Expression
    expr = expand(model(Xₛ, mdl.param)[1])
    # return clean_model(expr);
    return expr
end

# -----------------------------------------------------------------------------
# Get data

# file = jldopen("../numerical/data/sim_doe1.jld2")
file = jldopen("../numerical/data/sim_doe3.jld2", "r")

sims_all = gendata(file, 30, 0, Inf);
sims_all = gendata(file, 10, 1800, 4000);
sims_all = gendata(file, 20, 1000, 2000; remove_slow=0.9);
sims = sims_all[[1, 2]]

plotsims(sims)

# Flat data
sims_tmp = gendata(file, 10, 0, Inf);
sims_tmp = [s for s in sims_tmp if s["Rᵥ"] < 1e7]
sims = sims_tmp
# x = vcat([s["dp"] for s in sims]'...)';
x = vcat([s["all"] for s in sims]'...)';
y = vcat([s["Rᵥ"] for s in sims]'...)

sims2 = sims
close(file)

x = vcat([s["all"] for s in sims_all]'...)';
x̄, x̃ = mean(x; dims = 2), std(x; dims = 2)
Xs = (x .- x̄) ./ x̃

for (i, s) in enumerate(sims_all)
    s["scl"] = Xs[:, i]
end

sims = sims_all

# -----------------------------------------------------------------------------
# Outer layer - Initial tests to fit the resistance

# Model to be used on LsqFit
function outer(X̂, p̂, l)
    P = reshape(view(p̂, 1:(l[1] * l[2])), (l[1], l[2]))
    n = length(P)
    Q = reshape(view(p̂, (n + 1):(n + l[1] * l[2])), (l[1], l[2]))
    n += length(Q)
    W = reshape(view(p̂, (n + 1):(n + l[2] * l[3])), (l[2], l[3]))
    n += length(W)
    # I replaced the b with p̂[n+1] directly due to Flux.jl issues
    # b = view(p̂, n + 1);
    return (W' * (((X̂' * P) .* (X̂' * Q))') .+ p̂[n + 1])[:]
end

# To be used with LsqFit library
function fit(X, y, l)
    # Set LstFit model
    model(X̂, p̂) = outer(X̂, p̂, l)
    # Initial guess
    p₀ = ones(2 * (l[1] * l[2]) + l[2] * l[3] + 1) ./ 100
    # Fit
    mdl = curve_fit(model, X, y, p₀)
    return mdl
end

# Add more variables to X
X = vcat(x, ones(1, size(x, 2)));
# X = vcat(x, zeros(1, size(x, 2)));
X = vcat(X, 1.0 ./ x);
# Model dimensions
l = [size(X, 1), 2, 1];
# Initial guess
p₀ = ones(2 * (l[1] * l[2]) + l[2] * l[3] + 1) ./ 100;

lossrmse(y, outer(X, p₀, l))
cor(y, outer(X, p₀, l))

# Manual fit
mdl = fit(X, y, l);

ŷ = outer(X, mdl.param, l)
metrics(y, ŷ)

# -----------------------------------------------------------------------------
# Symbolic analysis

# Initialise symbolic variables
ρₛ, μₛ, dₛ, Lₛ, ϵₛ, aₛ, Aₛ = symbols("ρ, μ, d, L, ϵ, a, A"; real=true);

# Generate X
X = vcat(x, ones(1, size(x, 2)));
X = vcat(X, (1.0 ./ x[3, :] .^ 2)');
# X = vcat(X, x.^2);

# Adjust model
l = [size(X, 1), 2, 1];     # model dimensions
mdl = fit(X, y, l);
metrics(y, outer(X, mdl.param, l))

# Symbolic expression
syms = [ρₛ, μₛ, dₛ, Lₛ, ϵₛ, aₛ, Aₛ, 1, 1 / dₛ^2]
get_sym(mdl, l, syms)
# Expression directly from the function
out = outer(syms, mdl.param, l)[1]
expand(out * out)

# -----------------------------------------------------------------------------
# CP layer

# CP Tensor decomposition for flat based algorithms
function CP(X, p₀, rnk, order)
    nvar = size(X, 1)
    # Reshape the parameters to a tensor ∈ ℝⁿˣᵐˣᵖ where n is the number of 
    # variables, m the matrix rank and p the polynomio order.
    P = reshape(view(p₀, 1:(length(p₀) - rnk)), (nvar, rnk, order))
    r = view(p₀, (length(p₀) - rnk + 1):length(p₀))
    # First order
    Oₙ = X' * P[:, :, 1]
    # Second to n order
    for i in 2:(order - 1)
        Oₙ .*= X' * P[:, :, i]
    end
    return Oₙ * r
end

function fitcp(X, y, l)
    order, rnk = l
    # Set LstFit model
    model(X̂, p̂) = CP(X̂, p̂, rnk, order)
    # Initial guess
    p₀ = ones(rnk * size(X, 1) * order + rnk)
    # Fit
    mdl = curve_fit(model, X, y, p₀)
    return mdl
end

# Tests with symbolic
order = 4;
rnk = 1;
p₀ = ones(rnk * 3 * order + rnk);
# Get output
out = CP([μₛ, Lₛ, 1 / dₛ^2], p₀, rnk, order)
expand(out)

# Fit
# X = vcat(x, ones(1, size(x, 2)));
# X = vcat(X, (1.0./x[3, :].^2)');
X = copy(x);
X = vcat(x, ones(1, size(x, 2)));
X = vcat(X, 1.0 ./ x);

# Model parameters
order = 6;
rnk = 1;
l = [order, rnk]
# Set initial values
p₀ = ones(rnk * size(X, 1) * order + rnk) ./ 10
# Fit model
mdl = fitcp(X, y, l);
# Get approximation
ŷ = CP(X, mdl.param, rnk, order);
metrics(y, ŷ)

# -----------------------------------------------------------------------------
# Pipe differential equation mass and damper functions
nsim = sims[1]

# Damping approximation
# dθ = [1e6];
# d̂!(u, p) = @. p[1];

order = 6;
rnk = 1;
d̂!(X̂, p̂) = CP(X̂, p̂, rnk, order);
# dθ = ones(rnk*size(X, 1)*order + rnk)./1000;
dθ = ones(rnk * (length(nsim["all"]) + 1) * order + rnk)

# Mass approximation
# mθ = [1e6];
# m̂!(u, p) = @. p[1];

# mθ = [1.0, 1.0];
# m̂!(u, p) = @. p[1]*u[1] * u[2] / (p[2]*u[3]);

order = 6;
rnk = 1;
m̂!(X̂, p̂) = CP(X̂, p̂, rnk, order);
mθ = ones(rnk * (length(nsim["all"]) + 1) * order + rnk)

nnm = Chain(Dense((length(nsim["all"])*2 + 1), 5, relu), Dense(5, 1), sum, abs)
nnd = Chain(Dense((length(nsim["all"])*2 + 1), 5, relu), Dense(5, 1), sum, abs)

mθ, mre = Flux.destructure(nnm);
dθ, dre = Flux.destructure(nnd);

m̂!(X̂, p̂) = mre(p̂)(X̂);
d̂!(X̂, p̂) = dre(p̂)(X̂);

m = m̂!(vcat(nsim["all"], 1.0./nsim["all"], nsim["u0"]), mθ)
d = d̂!(vcat(nsim["all"], 1.0./nsim["all"], nsim["u0"]), dθ)
# m̂!(vcat(nsim["all"], nsim["u0"]), dθ)
# d̂!(vcat(nsim["all"], nsim["u0"]), dθ)

# -----------------------------------------------------------------------------
# Like SINDy

[1, 2, 3, 4, 5, 6, 7]
[ρ, μ, d, L, ϵ, a, A]

m̂!(x, p) = p[1]*(x[1]*x[4]/(x[7])) + p[2]*(x[2]*x[4]/(x[3]^4)) + p[3]*(x[1]*x[2]*x[3]/(x[6]*x[7])) + p[4]*(x[3]*x[6]/x[1]);
d̂!(x, p) = p[1]*(x[1]*x[4]/(x[7])) + p[2]*(x[2]*x[4]/(x[3]^4)) + p[3]*(x[1]*x[2]*x[3]/(x[6]*x[7])) + p[4]*(x[3]*x[6]/x[1]);
128/π
dθ = ones(4)
mθ = ones(4)


# -----------------------------------------------------------------------------
# Differential equations
# Train setup

function diffeq!(du, u, p, t)
    # return du[1] = min((p[1] - u[1] * d̂!(vcat(p[2], u), p[3])) / m̂!(vcat(p[4], u), p[5]), 20)
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
    mp = vcat(data["all"], 1.0./data["all"])
    dp = vcat(data["all"], 1.0./data["all"])
    # mp = data["all"]
    # dp = data["all"]
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
#     )
# end


function loss(p, prob)
    return sqrt(
            mean([mean(((s["y"] .- predict(p, s, prob))/s["y"][end]) .^ 2) for s in sims])
        )
end

# function loss(p, prob)
#     return mean([-cor(s["y"], predict(p, s, prob)) for s in sims]) + 0.001*sqrt(
#         mean([mean((s["y"] .- predict(p, s, prob)) .^ 2) for s in sims])
#     )
# end

# -----------------------------------------------------------------------------
# Generate gradient functions
∇d = gradient((p) -> loss(p, probd), dθ)
∇m = gradient((p) -> loss(p, probm), mθ)

# -----------------------------------------------------------------------------
# Training

optm = ADAM(0.0000005);
optd = ADAM(0.0000005);
optm = ADAM(0.0001);
optd = ADAM(0.0001);
optm = ADAM(100);
optd = ADAM(100);
optm = ADAM(0.001);
optd = ADAM(0.001);
# sims = sims2[[1, 6, 82, 8]]
# sims = sims2[[82]]
sims = sims_all
sims = sims_all[[1, 2, 3, 4, 5, 6, 7, 8, 9]]
sims = sims_all[[1, 3, 5, 8]]
sims = sims_all[[1, 2]]
plotsims(sims)

# mθ = [1e8]
# mθ = [1.0, 1.0]
# mθ = ones(rnk*size(X, 1)*order + rnk)./100;

# dθ = copy(mdl.param);
# dθ = ones(rnk*size(X, 1)*order + rnk)./10;
# dθ = [1e8]
# dθ = ones(rnk * (length(nsim["all"]) + 1) * order + rnk) ./ 10;
# mθ = ones(rnk * (length(nsim["all"]) + 1) * order + rnk) ./ 10;

nsim = sims[1]
# dp = vcat(nsim["all"], [1.0], 1.0 ./ nsim["all"])
info = Dict(
    # :Iᵥ => Dict(:pred => () -> m̂!(nsim["all"], mθ), :truth => nsim["Iᵥ"]),
    # :Rᵥ => Dict(:pred => () -> d̂!(dp, dθ), :truth => nsim["Rᵥ"]),
    :c => 0.0,
)

i = 1
λ = 5e-3
α = 50
it = 10000
hist = zeros(it)

# dθ = [0.0, 40.74, 0.0, 0]./10
# mθ = [1.0, 0, 0, 0.0]./10

# dθ .= [0.0, 0.0, 0.0, 1.0]
# mθ .= [0.0, 0.0, 0.0, 1.0]

# dθ .= [1.0, 0.0, 0.0, 1.0]
# mθ .= [0.0, 0.0, 0.0, 1.0]
# dθ .= ones(4)

# ∇d = gradient((p) -> loss(p, probd), dθ)

# dθ = Float64.(rand(4) .> 0.5)
# mθ = Float64.(rand(4) .> 0.5)
# dθ = ones(4)
# mθ = ones(4)
optm = ADAM(0.1);
optd = ADAM(0.1);
# plotsims(sims[[1]], dθ, probd)
# plot()
for i in 1:it
    ∇d = gradient((p) -> loss(p, probd), dθ)
    ∇m = gradient((p) -> loss(p, probm), mθ)

    # Zygote.hessian((p) -> loss(p, probd), dθ)

    update!(optm, mθ, ∇m[1])
    update!(optd, dθ, ∇d[1])
    
    # λ = i > 1 ? 10^(floor(log10(abs(hist[i - 1]))) - 2) : 5e-3
    # λ = i > 1 ? 10^(floor(log10(abs(hist[i - 1]))) - 4) : 5e-3
    # λ = 1e-4
    # λ = 1e-6
    # (norm(dθ) / norm(∇d[1])) * ∇d[1]
    # dθ += -(λ * norm(dθ) / norm(∇d[1])) * ∇d[1]
    # mθ += -(λ * norm(mθ) / norm(∇m[1])) * ∇m[1]

    hist[i] = loss(dθ, probd)

    # Stop criteria
    # if i > (α + 1)
    #     c = mean(diff(hist[(i - α):(i - 1)]))
    #     info[:c] = c
    #     if abs(c) < 1e-5
    #         # print("oi")
    #         break
    #     end
    # end

    # Print info it
    println(printit(i, hist[i]))
    # Plot 
    plotit(i, sims[1], () -> predict(dθ, sims[1], probd); every=5)
    # plotit(i, nsim, () -> predict(dθ, nsim, probd); every=10)
end

dθ
mθ

plotsims(sims, dθ, probd)
plotsims(sims[[4]], dθ, probd)


nsim = sims[1]
λ = nsim["y"][end]/predict(dθ, nsim, probd)[end]
cor(predict(dθ, nsim, probd), nsim["y"])
plot(nsim["y"])
plot!(predict(dθ, nsim, probd).*λ)

dθ
mθ

# Plot history
idx = findfirst(hist .== 0.0)
i = isnothing(idx) ? length(hist)+1 : idx
# plot(hist[1:(i - 1)]; yaxis=:log10)
plot(1:i-1 , hist[1:(i - 1)])

# Plot simulation
plotsims(sims, dθ, probd)
