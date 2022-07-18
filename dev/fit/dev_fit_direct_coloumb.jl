using Random
using CUDA
include("lib_dev_fit.jl")

# -----------------------------------------------------------------------------
# Get data

function diffeq!(du, u, p, t)

    Fc, Fs, νs, σ0, σ1, σ2, = p[:Fc], p[:Fs], p[:νs], p[:σ0], p[:σ1], p[:σ2]
    m, c, k, F = p[:m], p[:c], p[:k], p[:F]

    ẋ = u[1]
    x = u[2]
    z = u[3]
    
    j = 2;
    g = Fc + (Fs - Fc)*exp(-(abs(ẋ / νs)^j))
    ż = ẋ - σ0*(ẋ/g)*z
    Fl = σ0*z + σ1*ż + σ2*ẋ

    du[1] = (F(t) - c*ẋ - k*x - Fl)/m 
    du[2] = u[1]
    du[3] = ż
    du[4] = Fl
    # println(Fl)
end

pr = Dict()
pr[:σ0] = 1e5;
pr[:σ1] = sqrt(1e5);
pr[:σ2] = 0.4
pr[:Fc] = 1.0;
pr[:Fs] = 1.5;
pr[:νs] = 0.001;
pr[:m] = 1.0
pr[:c] = 0.0
pr[:k] = 2.0
pr[:F] = x -> min(5, x*η)
η = 0.5

sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8);
trng = 0:0.01:13
Y = sol(trng)[[1,2], :]
plot(trng, Y')
plot!(trng, predict(γ, probn, trng)', legend=false)

ts = 20.0
prob = ODEProblem(diffeq!, [0.0, 0.0, 0.0, 0.0], (0, ts), pr)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8);
# trng = 0:0.01:13
# solt = sol(trng)
# plot(solt[4, :], solt[1, :])
# plot!(trng, solt[1, :])
# plot(trng, solt[2, :])
# plot(trng, solt[4, :])
# sol = solve(prob, Tsit5(); saveat=0:0.01:10, reltol=1e-8, abstol=1e-8);
# sol = solve(prob)

trng = 0:1:13
trng = 0:0.01:13
Y = sol(trng)[[1,2], :]

plot(trng, Y[[1,2],:]')
plot!(trng, predict(γ, probn, trng)', legend=false)

Y4 = sol(trng)[4, :]
plot(trng, Y4)
plot!(trng, predictf(γ, probn, trng), legend=false)

# -----------------------------------------------------------------------------
# Approx functions

hardtanh(x) = max(-1, min(1, x))
# m = Chain(Dense(2, 10, tanh), Dense(10, 2))
# m = Chain(Dense(1, 5, x -> exp.(-(x.^2))), Dense(5, 1))
m = Chain(Dense(2, 5, tanh), Dense(5, 1))
# m = Chain(Dense(1, 5, tanh), Dense(5, 2))
γ, ml = Flux.destructure(m)

# F̂!(u, p) = ml(p)([u])[1];
F̂!(u, p) = ml(p)(u)[1];
# F̂!(u, p) = ml(p)(u);

# -----------------------------------------------------------------------------
# Differential equation setup

function diffeqn!(du, u, p, t)
    # x = F̂!(u[[1, 3]], p)
    # du[1] = (pr[:F](t) - pr[:c]*u[1] - pr[:k]*u[2] - x[1])/pr[:m]
    # du[2] = u[1]
    # du[3] = x[2]
    x = F̂!(vcat(u[1], pr[:F](t)), p)
    # x = F̂!(u[[1]], p)
    du[1] = (pr[:F](t) - pr[:c]*u[1] - pr[:k]*u[2] - x)/pr[:m]
    du[2] = u[1]
    # x = F̂!(u[1], p)
    # du[1] = (pr[:F](t) - pr[:c]*u[1] - pr[:k]*u[2] - x)/pr[:m]
    # du[2] = u[1]
    # du[3] = x
end

# probn = ODEProblem(diffeqn!, [0.0, 0.0, 0.0], (0.0, 50.0), γ)
probn = ODEProblem(diffeqn!, [0.0, 0.0], (0.0, 50.0), γ)

# -----------------------------------------------------------------------------
# Train setup

function predict(p, prob, t=nothing)
    saveat = isnothing(t) ? trng : t
    tmp_prob = remake(prob; p = p)
    out = solve(tmp_prob, Tsit5(); saveat=saveat, reltol=1e-8, abstol=1e-8)
    # return vec(out[[1, 2], :])
    return out[[1, 2], :]
end

function predictf(p, prob, t=nothing)
    saveat = isnothing(t) ? trng : t
    tmp_prob = remake(prob; p = p)
    out = solve(tmp_prob, Tsit5(); saveat=saveat, reltol=1e-8, abstol=1e-8)
    # return vec(out[[1, 2], :])
    return vec(out[3, :])
end

# Define loss functions

function loss(p, prob, Y)
    Ŷ = predict(p, prob)
    # ρ = -cor(Ŷ[:,], Y)
    ρ = -mean(cor.(eachrow(Y), eachrow(Ŷ)))
    rmse = sqrt(mean((Y - Ŷ).^2))
    return ρ + rmse
end

function lossout(p, prob, Y)
    Ŷ = predict(p, prob)
    ρ = -mean(cor.(eachrow(Y), eachrow(Ŷ)))
    rmse = sqrt(mean((Y - Ŷ).^2))
    return ρ, rmse
end

# Generate gradient functions
∇γ = gradient((p) -> loss(p, probn, Y), γ)

# -----------------------------------------------------------------------------
# Train

tik = 0.0

info = Dict(
    :s => 0,
    :ρ => 0.0,
    :rmse => 0.0,
    :tik => tik,
    :eta => 0.0,
)

pr[:F] = x -> min(5, x*0.75)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8);
trng = 0:0.01:13
Y = sol(trng)[[1,2], :]

opt = ADAM(0.03)
it = 400000
α = 3000
hist = zeros(it)
tik = @elapsed for i = 1:it
    
    ∇γ = gradient((p) -> loss(p, probn, Y), γ)

    update!(opt, γ, ∇γ[1])

    ρc, rmse  = lossout(γ, probn, Y)
    hist[i] = ρc + rmse
    info[:ρ], info[:rmse] = ρc, rmse

    # Stop criteria
    if stopcrit!(i, hist, α)
        break
    elseif (ρc < -(1-1e-4)) & (rmse < 1e-2)
        break
    end

    # Print info it
    println(printit(i, hist[i], info))
end

trng2 = sol.t
trng2 = trng
plot(trng2, sol(trng2)[2, :])
plot!(trng2, predict(γ, probn, trng2)[2, :], legend=false)

plot(trng2, sol(trng2)[1, :])
plot!(trng2, predict(γ, probn, trng2)[1, :], legend=false)


plot(sol(trng2)[2, :], sol(trng2)[1, :])
plot!(predict(γ, probn, trng2)[2, :], predict(γ, probn, trng2)[1, :], legend=false)




