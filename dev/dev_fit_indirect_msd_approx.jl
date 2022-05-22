using Random
using CUDA

include("lib_dev_fit.jl")

# For file history go to the fit-hist/dev_fit_indirect_md_approx.jl. 
# This is a cleaned file version of it, with some solved previous issues.

# -----------------------------------------------------------------------------
# Load data

# Parameters file
filefit = "./data/fit_msd_sim_doe5.jld2";

# Set reynolds filter
filesim = "../numerical/data/sim_doe5.jld2"
re_rng = [0, 2200]

# Load 
x, Y, losslist, keylist = loadθ(filefit; filesim=filesim, re_rng=re_rng);

# -----------------------------------------------------------------------------
# Prepare data

x = Float32.(x);
Y = Dict(s => Float32.(Y[s]) for s in keys(Y))
id = collect(1:size(x, 2));

# High error simulations

msk = id[losslist .< -0.98];

# Apply mask
x = x[:, msk];
Y = Dict(s => Y[s][repeat([:],length(size(Y[s])) - 1)..., msk] for s in keys(Y))
losslist = losslist[msk];
keylist = keylist[msk];
id = collect(1:size(x, 2));

# -----------------------------------------------------------------------------
# Data transformations

# Log transformation on log spaced variables
x[2, :] = log10.(x[2, :]);
x[3, :] = log10.(x[3, :]);
x[5, :] = log10.(x[5, :]);
x[7, :] = log10.(x[7, :]);

# Scale data to 0 mean and 1 std
x̄, x̃ = mean(x; dims = 2), std(x; dims = 2);
Xs = (x .- x̄) ./ x̃;

Ȳ = Dict(s => mean(Y[s], dims=length(size(Y[s]))) for s in keys(Y));
Ỹ = Dict(s => std(Y[s], dims=length(size(Y[s]))) for s in keys(Y));
Ys = Dict(s => (Y[s] .- Ȳ[s])./Ỹ[s] for s in keys(Y));


# -----------------------------------------------------------------------------
# Data split

# train percentage = β -> test percentage = (1 - β)
β = 0.7;

# Random split
idxs = shuffle(id);
ntrn = floor(Int, β * length(id));
trn_idx = idxs[1:ntrn];
tst_idx = idxs[ntrn+1:end];

# Split by the scaled value to guarantee interpolations
lvl = 1.5;
lvl = 1.3;
msk = hcat([-lvl .< Ys[p] .< lvl for s in keys(Ys)]...)
msk = hcat([all(c) for c in eachcol(-lvl .< Xs .< lvl)], msk)
tst_idx = id[vec(sum(msk, dims=2) .== 4)];

if length(tst_idx) > (1 - β)*length(id)
    n_tst = floor(Int, (1 - β) * length(id))
    tst_idx = shuffle(tst_idx)[1:n_tst]
else
    trn_idx = collect(setdiff(Set(id), Set(tst_idx)))
end

# Check x split visually
scatter(Xs[1, trn_idx], Xs[2, trn_idx])
scatter!(Xs[1, tst_idx], Xs[2, tst_idx])

scatter(Xs[1, trn_idx], Xs[3, trn_idx])
scatter!(Xs[1, tst_idx], Xs[3, tst_idx])

scatter(Xs[1, trn_idx], Xs[4, trn_idx])
scatter!(Xs[1, tst_idx], Xs[4, tst_idx])

scatter(Xs[1, trn_idx], Xs[5, trn_idx])
scatter!(Xs[1, tst_idx], Xs[5, tst_idx])

scatter(Xs[1, trn_idx], Xs[6, trn_idx])
scatter!(Xs[1, tst_idx], Xs[6, tst_idx])

scatter(Xs[1, trn_idx], Xs[7, trn_idx])
scatter!(Xs[1, tst_idx], Xs[7, tst_idx])

# -----------------------------------------------------------------------------
# Send data to GPU

p = :s;
# Train dataset
xs = gpu(Xs[:, trn_idx]);
# ys = gpu(Ys[p][:, trn_idx]);
ys = gpu(Ys[p][trn_idx]);

# Test dataset
xt = gpu(Xs[:, tst_idx]);
# yt = gpu(Ys[p][:, tst_idx]);
yt = gpu(Ys[p][tst_idx]);

# -----------------------------------------------------------------------------
# Train

# Define model
nn = gpu(Chain(Dense(size(xs, 1) => 14, tanh), Dense(14 => 1)))
nn = gpu(Chain(Dense(size(xs, 1) => 8, tanh), Dense(8 => 1)))

ps = Flux.params(nn);

# Optimizer
opt = ADAM(0.01);

# Loss function
loss(x, y) = sqrt(Flux.Losses.mse(nn(x), y));
loss(xs, ys')

info = gpu(Dict(:c => 0.0))
it = 900000
hist = gpu(zeros(it, 2))
for i = 1:it
    # batch = gpu(shuffle(1:size(xs, 2))[1:30])
    batch = 1:size(xs, 2)
    ∇p = gradient(() -> loss(xs[:, batch], ys[batch]'), ps)

    update!(opt, ps, ∇p)

    hist[i, 1] = loss(xs, ys')
    hist[i, 2] = loss(xt, yt')
    info[:c] = hist[i, 2]

    # Stop criteria
    # if stopcrit!(c, i, hist, α) break end
    # oi = stopcrit!(c, i, hist, α)
    # info[:c] = c

    # Print info it
    println(printit(i, hist[i], info))
end

# -----------------------------------------------------------------------------
# Analyse

ŷ = vec(nn(xt));
metrics(ŷ, vec(yt))

ŷ = vec(nn(xs));
metrics(ŷ, vec(ys))

# Plot history
idx = findfirst(hist[:, 1] .== 0.0)
i = isnothing(idx) ? size(hist, 1) + 1 : idx
plot(hist[1:(i-1), 1]; yaxis = :log10)
plot!(hist[1:(i-2), 2]; yaxis = :log10)

# -----------------------------------------------------------------------------
# Store info

out = Dict(:y_mean => Ȳ, :y_std => Ỹ, :x_mean => x̄, :x_std => x̃)
out[:nnm] = deepcopy(cpu(nn))
out[:nnd] = deepcopy(cpu(nn))
out[:nns] = deepcopy(cpu(nn))

savezlib("./data/model_doe5.zdat", out)