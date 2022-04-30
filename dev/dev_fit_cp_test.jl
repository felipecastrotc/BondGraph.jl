include("lib_dev_fit.jl")
using BSON: @load

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
# y = vcat([s["Rᵥ"] for s in sims]'...)
x = vcat([s["all"] for s in sims]'...)';

sims2 = sims
close(file)
# -----------------------------------------------------------------------------
# Outer layer

# Model to be used on LsqFit
function model_cp(X̂, p̂, l)
    P = reshape(view(p̂, 1:l[1]*l[2]), (l[1], l[2]));
    n = length(P);
    Q = reshape(view(p̂, (n + 1):(n + l[1]*l[2])), (l[1], l[2]));
    n += length(Q);
    W = reshape(view(p̂, (n + 1):(n + l[2]*l[3])), (l[2], l[3]));
    n += length(W);
    # I replaced the b with p̂[n+1] directly due to Flux.jl issues
    # b = view(p̂, n + 1);
    return (W' * (((X̂' * P) .* (X̂' * Q))') .+ p̂[n+1])[:];
end

function fit(X, y, l)
    # Set LstFit model
    model(X̂, p̂) = model_cp(X̂, p̂, l);
    # Initial guess
    p₀ = ones(2*(l[1]*l[2]) + l[2]*l[3] + 1)./100;
    # Fit
    mdl = curve_fit(model, X, y, p₀);
    return mdl
end

# Add more variables to X
# X = vcat(x, ones(1, size(x, 2)));
X = vcat(x, zeros(1, size(x, 2)));
X = vcat(X, 1.0./x);
# Model dimensions
l = [size(X, 1), 2, 1];
# Initial guess
p₀ = ones(2*(l[1]*l[2]) + l[2]*l[3] + 1)./100;

lossmse(y, model_cp(X, p₀, l))
cor(y, model_cp(X, p₀, l))

# Manual fit
mdl = fit(X, y, l);

p = mdl.param
ŷ = model_cp(X, p, l)
lossmse(y, ŷ)
cor(y, ŷ)
scatter(y, ŷ)

# -----------------------------------------------------------------------------
# Symbolic analysis
using SymPy

# mp = [ρ, L, A]
# dp = [μ, L, d]
# Initialise symbolic variables
(μₛ, Lₛ, dₛ, Aₛ) = symbols("μ, L, d, A", real=true);

function clean_model(mdl)
  tmp = Dict()
  for i in mdl.atoms()
    if !(typeof(N(i))<:Sym)
      tmp[i] = round(i, digits=6);
    end
  end
  return mdl.xreplace(tmp);
end

function get_sym(mdl, l, Xₛ)
    # Set LstFit model
    model(X̂, p̂) = model_cp(X̂, p̂, l);
    # Expression
    expr = expand(model(Xₛ, mdl.param)[1]);
    # return clean_model(expr);
    return expr;
end

# Define the system 1p₀
y = f(x[1, :], x[2, :], x[3, :]);
# Generate X
X = vcat(x, ones(1, size(x, 2)));
X = vcat(X, (1.0./x[3, :].^2)');
# X = vcat(X, x.^2);

# Adjust model
l = [size(X, 1), 2, 1];     # model dimensions
mdl = fit(X, y, l);
lossmse(y, model_cp(X, mdl.param, l))
get_sym(mdl, l, [μₛ, Lₛ, dₛ, 1, 1/dₛ^2,])

out = model_cp([μₛ, Lₛ, dₛ, 1, 1/dₛ^2,], mdl.param, l)[1]
# get_sym(mdl, l, [μₛ, Lₛ, dₛ, 1, 1/μₛ, 1/Lₛ, 1/dₛ, μₛ^2, Lₛ^2, dₛ^2])

out
expand(out*out)

# -----------------------------------------------------------------------------
# CP layer

# CP Tensor decomposition for flat based algorithms
function CP(X, p₀, rnk, order)
    nvar = size(X, 1)
    # Reshape the parameters to a tensor ∈ ℝⁿˣᵐˣᵖ where n is the number of 
    # variables, m the matrix rank and p the polynomio order.
    P = reshape(view(p₀, 1:(length(p₀) - rnk)), (nvar, rnk, order))
    r = view(p₀, (length(p₀)-rnk +1):length(p₀))
    # First order
    Oₙ = X' * P[:, :, 1];
    # Second to n order
    for i in 2:order-1
        Oₙ .*= X' * P[:, :, i];
    end
    return Oₙ * r
end

function fitcp(X, y, l)
    order, rnk = l
    # Set LstFit model
    model(X̂, p̂) = CP(X̂, p̂, rnk, order);
    # Initial guess
    p₀ = ones(rnk*size(X, 1)*order + rnk);
    # Fit
    mdl = curve_fit(model, X, y, p₀);
    return mdl
end

function get_symcp(mdl, l, Xₛ)
    order, rnk = l
    # Set LstFit model
    model(X̂, p̂) = CP(X̂, p̂, rnk, order);
    # Expression
    expr = expand(model(Xₛ, mdl.param)[1]);
    # return clean_model(expr);
    return expr;
end

# Tests with symbolic
order = 4;
rnk = 1;
p₀ = ones(rnk*size(X, 1)*order + rnk);
# Get output
out = CP([μₛ, Lₛ, 1/dₛ^2], p₀, rnk, order)
expand(out)

# Fit
# X = vcat(x, ones(1, size(x, 2)));
# X = vcat(X, (1.0./x[3, :].^2)');
X = copy(x);
X = vcat(x, ones(1, size(x, 2)));
X = vcat(X, 1.0./x);

# Model parameters
order = 6;
rnk = 1;
l = [order, rnk]
# Set initial values
p₀ = ones(rnk*size(X, 1)*order + rnk)./10
# Fit model
y = Float64.(mlist)
y = Float64.(dlist)
mdl = fitcp(X, y, l);
# Get approximation
ŷ = CP(X, mdl.param, rnk, order);
metrics(y, ŷ)

# -----------------------------------------------------------------------------
# CP setup

# order = 7;
# rnk = 3;
# Set LstFit model
d̂!(X̂, p̂) = CP(X̂, p̂, rnk, order);

# Initial guess
# p₀ = ones(rnk*size(X, 1)*order + rnk)./1000;
dθ = copy(mdl.param);

mθ = [1.0, 1.0];
m̂!(u, p) = @. p[1]*u[1] * u[2] / (p[2]*u[3]);
m̂!(X̂, p̂) = CP(X̂, p̂, rnk, order);

# -----------------------------------------------------------------------------
# CP setup

mθ = [1e6];
dθ = [1e6];
m̂!(u, p) = @. p[1];
d̂!(u, p) = @. p[1];

# -----------------------------------------------------------------------------
# NN import setup

@load "./data/model_sim_doe3.bson" nn

m̂!(u, p) = (nn((u .- mean(x, dims=2))./std(x, dims=2)).* std(y) .+ mean(y))[1];
d̂!(u, p) = u;
mθ = [1e6];
dθ = [1e6];

m̂!(sims[1]["all"], 1)

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
    # mp = vcat(data["mp"], [1], 1.0./data["mp"])
    # dp = vcat(data["dp"], [1], 1.0./data["dp"])
    Δp = data["Δp"]
    # diffeqd!(du, u, p, t) = diffeq!(du, u, [data["dp"], p, data["mp"], mθ, data["Δp"]], t)
    # diffeqm!(du, u, p, t) = diffeq!(du, u, [data["dp"], dθ, data["mp"], p, data["Δp"]], t)
    vec(solve(tmp_prob, Tsit5(), saveat=data["trng"]))./data["A"]
end

# Test
function predict(p, data, prob)
    global mp, dp, Δp
    tmp_prob = remake(prob, u0=data["u0"], p=p)
    mp = data["all"]
    # dp = data["dp"]
    Δp = data["Δp"]
    # vec(solve(tmp_prob, Tsit5(), saveat=data["trng"]))./data["A"]
    vec(solve(tmp_prob, Tsit5(), saveat=0:0.01:5))./data["A"]
end

sims = sims2
ns = 44
s = sims[ns]
plot(s["t"], s["y"], label=".")
# plot!(s["t"], predict(mdl.param, s, probd), label="~")
dp = filef["d"][ns]
plot!(0:0.01:5, predict(dθ, s, probd), label="~")
s["Re"]

abs(d̂!(vcat(s["dp"], [1], 1.0./s["dp"]), mdl.param) - s["Rᵥ"])/s["Rᵥ"]

[i for (i, s) in enumerate(sims2) if s["Re"] < 2400]

# -----------------------------------------------------------------------------
# Define loss functions

function loss(p, prob)
    # return sum((s["y"] .- predict_neuralode(p, s)) .^ 2)/n
    return sqrt(mean([mean((s["y"] .- predict(p, s, prob)) .^ 2) for s in sims]))
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
# sims2 = copy(sims)
# sims = sims2[[1, 6, 82, 8]]
sims = sims2[[82]]
plotsims(sims)

mθ = [1.0, 1.0]
dθ = copy(mdl.param);
dθ = ones(rnk*size(X, 1)*order + rnk)./10;
mθ = ones(rnk*size(X, 1)*order + rnk)./100;

mθ = [1e8]
dθ = [1e8]

it = 40000;
hist = zeros(it);
for i = 1:it
    # i = 2
    ∇d = gradient((p) -> loss(p, probd), dθ)
    ∇m = gradient((p) -> loss(p, probm), mθ)

    # update!(optm, mθ, ∇m[1])
    # update!(optd, dθ, ∇d[1])
    # dθ += -(5e-2*norm(dθ)/norm(∇d[1]))*∇d[1]
    dθ += -(5e-3*norm(dθ)/norm(∇d[1]))*∇d[1]
    # dθ += -(1e-3*norm(dθ)/norm(∇d[1]))*∇d[1]
    # dθ += -(1e-6*norm(dθ)/norm(∇d[1]))*∇d[1]
    # dθ += -(1e-10*norm(dθ)/norm(∇d[1]))*∇d[1]
    # mθ += -(1e-1*norm(mθ)/norm(∇m[1]))*∇m[1]
    mθ += -(5e-3*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-3*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-6*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-7*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-8*norm(mθ)/norm(∇m[1]))*∇m[1]
    # mθ += -(1e-2*norm(mθ)/norm(∇m[1]))*∇m[1]

    hist[i] = loss(dθ, probd)
    ns = 1
    mp = vcat(sims[ns]["mp"], [1], 1.0./sims[ns]["mp"])
    dp = vcat(sims[ns]["dp"], [1], 1.0./sims[ns]["dp"])
    Iᵥ = sims[ns]["Iᵥ"]
    Rᵥ = sims[ns]["Rᵥ"]
    @printf("It: %d - loss: %.9e m: %.5e/%.5e d: %.5e/%.5e\n", i, hist[i], m̂!(mp, mθ), Iᵥ,d̂!(dp, dθ), Rᵥ)
    # @printf("It: %d - loss: %.9e\n", i, hist[i])

    if i == 1
        p = plot(sims[ns]["t"], sims[ns]["y"])
    elseif ((i % 500) == 0)
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


plotsims(sims[[1]], dθ, probd)

i = findfirst(hist .== 0.0)
# i = length(hist)
plot(hist[1:i-1], yaxis=:log10)
# plot(hist, yaxis=:log10)
plotsims(sims, dθ, probd)


# ----------------------------------

mlist = []
dlist = []
losslist = []
for s in 1:length(sims2)
    global sims, mθ, dθ

    mθ = [1e6]
    dθ = [1e6]

    sims = sims2[[s]]

    λ = 5e-3;
    α = 50;
    it = 40000;
    hist = zeros(it);
    for i = 1:it
        # i = 2
        ∇d = gradient((p) -> loss(p, probd), dθ)
        ∇m = gradient((p) -> loss(p, probm), mθ)

        
        λ = i > 1 ? 10^(floor(log10(hist[i - 1])) - 2) : 5e-3

        dθ += -(5e-3*norm(dθ)/norm(∇d[1]))*∇d[1]
        mθ += -(5e-3*norm(mθ)/norm(∇m[1]))*∇m[1]

        hist[i] = loss(dθ, probd)

        if i > (α+1)
            c = mean(diff(hist[i-α:i-1]))
            @printf("s: %d - It: %d - loss: %.9e - c: %.9e\n", s, i, hist[i], c)
            # if (c < 1e-2) & (c > -1e-2)
            if abs(c) < 1e-5
                break
            end
        else
            @printf("It: %d - loss: %.9e \n", i, hist[i])
        end

    end
    i = findfirst(hist .== 0.0)
    push!(mlist, mθ[1])
    push!(dlist, dθ[1])
    push!(losslist, hist[i-1])
end

plot(mlist)
plot(dlist)
# plotsims(sims, dθ, probd)
i = findfirst(hist .== 0.0)
# plot(hist[1:i-1])
# plot(diff(hist[1:i-1]))
# mean(diff(hist[i-α:i-1]))

extrema(losslist)
histogram(log10.(losslist))

10^((mean(log10.(losslist)) + median(log10.(losslist)))/2)
(10^((mean(log10.(losslist)) + median(log10.(losslist)))/2))^2


filename = "./data/fit_md_sim_doe3.jld2"
file = jldopen(filename, "w")

file["x"] = x
file["m"] = Float64.(mlist)
file["d"] = Float64.(dlist)
file["loss"] = Float64.(losslist)

close(file)