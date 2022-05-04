using DiffEqFlux, DifferentialEquations, Plots
using JLD2, LinearAlgebra, Statistics

using Flux
using DiffEqSensitivity
using Flux.Optimise: ADAM, update!
using Printf

# Generate data
g = 9.81        # m/s^2 - Gravity

# -----------------------------------------------------------------------------
# Flow functions

function calcRe(V, ρ, μ, d)
    return ρ .* V .* d ./ μ
end

# -----------------------------------------------------------------------------
# Data functions

function getre(file, minre=1000, maxre=10000)
    out = Dict()
    for i in keys(file)
        ex = file[i]
        ρ, μ, d, L, ϵ, a = ex["ρ"], ex["μ"], ex["d"], ex["L"], ex["ϵ"], ex["a"]

        A = π * d^2 / 4
        # Rᵥ = 128 * μ * L / (π * d^4)
        # Iᵥ = ρ * L / A

        Hin = ex["Hr"] - ex["Hf"]

        Re = calcRe(ex["V"][end], ρ, μ, d)
        if (Re > minre) & (Re < maxre)
            # println(Re, "   ", i)
            out[i] = Re
        end
    end

    return out
end

function gendata(file, n=10, minre=1000, maxre=10000; remove_slow=nothing)
    K = keys(getre(file, minre, maxre))

    if !isnothing(remove_slow)
        K = [k for k in K if cor(file[k]["V"], file[k]["t"]) < remove_slow]
    end

    # Res = ["162", "163"]
    # Res = ["162"]
    sims = []
    for k in K
        # for k in Res
        data = Dict()

        ex = file[k]
        ρ, μ, d, L, ϵ, a = ex["ρ"], ex["μ"], ex["d"], ex["L"], ex["ϵ"], ex["a"]
        A = π * d^2 / 4

        # Neural networks inputs
        data["A"] = A
        # data["mp"] = [ρ, L, d, A]
        # data["dp"] = [μ, ϵ, L, d, A]
        data["mp"] = [ρ, L, A]
        data["dp"] = [μ, L, d]
        data["Δp"] = (ex["Hr"] - ex["Hf"]) * ρ * g
        data["Re"] = calcRe(ex["V"][end], ρ, μ, d)
        data["all"] = [ρ, μ, d, L, ϵ, a, A]

        data["Rᵥ"] = 128 * μ * L / (π * d^4)
        data["Iᵥ"] = ρ * L / A

        # Simulation settings
        tsim = ex["t"] .< 4
        stp = Int(round(length(ex["t"][tsim]) / n))
        rng = 1:stp:length(ex["t"][tsim])
        y = ex["V"][rng]
        t = ex["t"][rng]
        trng = 0:step(t):(step(t) * (length(rng) - 1))

        data["u0"] = [y[1] .* A]
        data["y"] = y
        data["t"] = t
        data["trng"] = trng

        push!(sims, data)
    end
    return sims
end

function scaledata(sims)
    # Scaling
    for p in ["mp", "dp"]
        for i in 1:length(sims[1][p])
            x = [s[p][i] for s in sims]
            x̄ = mean(x)
            ẋ = std(x)
            for s in sims
                s[p][i] = (s[p][i] - x̄) / ẋ
            end
        end
    end
end

# -----------------------------------------------------------------------------
# Plot functions

function plotsims(sims, p=missing, prob=missing)
    plt = plot()
    for (i, s) in enumerate(sims)
        # Color palette selection
        c = palette(:default)[i % 15 + 1]
        # Plot the reference
        plt = plot!(s["t"], s["y"]; ls=:dash, linecolor=c, label="Ref."*string(i))
        # Simulate and plot the simulation
        if !ismissing(p)
            plot!(s["t"], predict(p, s, prob); linecolor=c, label="Sim."*string(i))
        end
    end
    return plot!()
end

# -----------------------------------------------------------------------------
# Fit functions

# Set loss and model functions
lossrmse(x, y) = sqrt(mean((x - y) .^ 2))

# Metrics
function metrics(y, ŷ)
    rmse = lossrmse(y, ŷ)
    a, b = y \ hcat(ŷ, ones(length(y)))
    ρ = cor(y, ŷ)

    err = (1 .- abs.((y - ŷ) ./ y)) .* 100

    @printf(
        "rmse: %.4e mse: %.4e ρ: %.5e a: %.5e b: %.5e ē: %.5e em: %.5e\n",
        rmse,
        rmse^2,
        ρ,
        a,
        b,
        mean(err),
        median(err)
    )

    return scatter(y, ŷ, ylabel="Predicted", xlabel="True", legend=false)
end

function plotit(i, t, y, ŷ; every=500)
    if i == 1
        p = plot(t, y; label=string(i))
    elseif (i % every) == 0
        p = plot!(t, ŷ; label=string(i))
        display(p)
    end
end

function plotit(i::Int, data::Dict, ŷ::Vector; every=500)
    return plotit(i, data["t"], data["y"], ŷ; every=every)
end

function printit(i, loss, info=missing)
    # rep = 
    base = @sprintf("It: %d - Loss: %.9e ", i, loss)
    if isa(info, Dict)
        for k in keys(info)
            aux = ""
            if isa(info[k], Dict)
                aux *= @sprintf("%.3e/%.3e", info[k][:pred](), info[k][:truth])
            else
                if isa(info[k], Number)
                    aux *= @sprintf("%.3e", info[k])
                else
                    aux *= @sprintf("%.3e", info[k]())
                end
            end
            base *= @sprintf("- %s: %s ", string(k), aux)
        end
    else
        base *= @sprintf(" %.3e", info[k])
    end
    return base
end

function stopcrit!(i, hist, α=50, tol=1e-5)
    if i > (α + 1)
        c = mean(diff(hist[(i - α):(i - 1)]))
        return abs(c) < tol
    else
        return false
    end
end