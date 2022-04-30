using DiffEqFlux, DifferentialEquations, Plots
using JLD2, LinearAlgebra, Statistics

using Flux, LsqFit
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

function gendata(file, n=10, minre=1000, maxre=10000)
    Res = getre(file, minre, maxre)
    # Res = ["162", "163"]
    # Res = ["162"]
    sims = []
    for k in keys(Res)
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
        trng = 0:step(t):step(t)*(length(rng)-1)

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
                s[p][i] = (s[p][i] - x̄)/ẋ
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
        c = palette(:default)[i%15 + 1]
        # Plot the reference
        plt = plot!(s["t"], s["y"], ls=:dash, linecolor=c)
        # Simulate and plot the simulation
        if !ismissing(p)
            plot!(s["t"], predict(p, s, prob), linecolor=c)
        end
    end
    plot!()
end