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

sims_all = gendata(file, 10, 0, Inf);
sims_tmp = gendata(file, 10, 1000, 2000);
sims_tmp = gendata(file, 10, 1800, 4000);
sims = sims_tmp[[1, 2]]

plotsims(sims_tmp)

# Flat data
sims_tmp = gendata(file, 10, 0, Inf);
sims_tmp = [s for s in sims_tmp if s["Rᵥ"] < 1e7]
sims = sims_tmp
# x = vcat([s["dp"] for s in sims]'...)';
x = vcat([s["all"] for s in sims]'...)';
y = vcat([s["Rᵥ"] for s in sims]'...)

sims2 = sims
close(file)

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

# Damping approximation
# dθ = [1e6];
# d̂!(u, p) = @. p[1];

order = 6;
rnk = 1;
d̂!(X̂, p̂) = CP(X̂, p̂, rnk, order);
# dθ = ones(rnk*size(X, 1)*order + rnk)./1000;
dθ = copy(mdl.param);

# Mass approximation
mθ = [1e6];
m̂!(u, p) = @. p[1];

# mθ = [1.0, 1.0];
# m̂!(u, p) = @. p[1]*u[1] * u[2] / (p[2]*u[3]);

# m̂!(X̂, p̂) = CP(X̂, p̂, rnk, order);

# -----------------------------------------------------------------------------
# Differential equations
# Train setup

function diffeq!(du, u, p, t)
    return du[1] = (p[1] - u[1] * d̂!(p[2], p[3])) / m̂!(p[4], p[5])
end

data = sims[1]
mp = data["mp"]
dp = data["dp"]
Δp = data["Δp"]
u0 = data["u0"]

diffeqd!(du, u, p, t) = diffeq!(du, u, [Δp, dp, p, mp, mθ], t)
diffeqm!(du, u, p, t) = diffeq!(du, u, [Δp, dp, dθ, mp, p], t)

probm = ODEProblem(diffeqm!, u0, (0.0, 20.0), mθ)
probd = ODEProblem(diffeqd!, u0, (0.0, 20.0), dθ)

# Test
function predict(p, data, prob)
    global mp, dp, Δp
    tmp_prob = remake(prob; u0=data["u0"], p=p)
    # mp = data["all"]
    # dp = data["all"]
    mp = data["mp"]
    # dp = data["dp"]
    dp = vcat(data["all"], [1.0], 1.0 ./ data["all"])
    Δp = data["Δp"]
    return vec(solve(tmp_prob, Tsit5(); saveat=data["trng"])) ./ data["A"]
end

# -----------------------------------------------------------------------------
# Define loss functions

function loss(p, prob)
    return sqrt(
        mean([mean((s["y"] .- predict(p, s, prob)) .^ 2) for s in sims])
    )
end

# -----------------------------------------------------------------------------
# Generate gradient functions
∇d = gradient((p) -> loss(p, probd), dθ)

# -----------------------------------------------------------------------------
# Training

optm = ADAM(0.0000005);
optd = ADAM(0.0000005);
optm = ADAM(0.0001);
optd = ADAM(0.0001);
optm = ADAM(100);
optd = ADAM(100);
# sims = sims2[[1, 6, 82, 8]]
sims = sims2[[82]]
plotsims(sims)

mθ = [1e8]
mθ = [1.0, 1.0]
# mθ = ones(rnk*size(X, 1)*order + rnk)./100;

dθ = copy(mdl.param);
# dθ = ones(rnk*size(X, 1)*order + rnk)./10;
# dθ = [1e8]

nsim = sims[1]
dp = vcat(nsim["all"], [1.0], 1.0 ./ nsim["all"])
info = Dict(
    :Iᵥ => Dict(:pred => () -> m̂!(nsim["all"], mθ), :truth => nsim["Iᵥ"]),
    :Rᵥ => Dict(:pred => () -> d̂!(dp, dθ), :truth => nsim["Rᵥ"]),
    :c => 0.0,
)

λ = 5e-3
α = 50
it = 10000
hist = zeros(it)
for i in 1:it
    # i = 2
    ∇d = gradient((p) -> loss(p, probd), dθ)
    ∇m = gradient((p) -> loss(p, probm), mθ)

    # update!(optm, mθ, ∇m[1])
    # update!(optd, dθ, ∇d[1])
    
    # λ = i > 1 ? 10^(floor(log10(hist[i - 1])) - 2) : 5e-3
    λ = 5e-3
    dθ += -(λ * norm(dθ) / norm(∇d[1])) * ∇d[1]
    mθ += -(λ * norm(mθ) / norm(∇m[1])) * ∇m[1]

    hist[i] = loss(dθ, probd)

    # Stop criteria
    if i > (α + 1)
        c = mean(diff(hist[(i - α):(i - 1)]))
        info[:c] = c
        if abs(c) < 1e-5
            break
        end
    end

    # Print info it
    printit(i, hist[i], info)
    # Plot 
    plotit(i, nsim, ŷ; every=500)
end

sims[1]

plotsims(sims[[1]], dθ, probd)

# Plot history
idx = findfirst(hist .== 0.0)
i = isnothing(idx) ? length(hist)+1 : idx
plot(hist[1:(i - 1)]; yaxis=:log10)

# Plot simulation
plotsims(sims, dθ, probd)
