using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D
using DifferentialEquations
using ModelingToolkit
using XLSX, DataFrames, Polynomials
using Plots, LsqFit, YAML

# Path used to store valve info
PATH_VALVE = "/Users/fctc/me/doutorado/dev/pump/valve-docs/"

# Valve data
xlsx = XLSX.readxlsx(PATH_VALVE * "info_data.xlsx")
vlv = DataFrame(xlsx["data"][:][:, 2:end]', xlsx["data"][:][:, 1])
vlv = float.(vlv)

col = "Fd"
col = "Cv"
pol = fit(vlv[!, "prct"], vlv[!, col], 5)

# Graph analysis of Cv polynomial fit
x = 0.1:0.001:1
plot(x, pol.(x))
plot!(vlv[!, "prct"], vlv[!, col])

# Store info
yml_file = PATH_VALVE * "fit_valve_pol_params.yml"
ps = YAML.load_file(yml_file)
ps[Symbol(col)] = pol.coeffs
f = open(PATH_VALVE * "fit_valve_pol_params.yml", "w")
YAML.write(f, ps)
close(f)

# Params to calculate the Reynolds of the valve
N₂ = 2.14 * 1e-3
N₄ = 7.6 * 1e-2
Fd = 0.32
Q = 63.75
d = 3 * 25.4
Fl = 0.82
C = 4.32
μ = 0.189
ρ = 850

Re = ((ρ * N₄ * Fd * Q) / (μ * sqrt(C * Fl))) * ((Fl^2 * C^2) / (N₂ * d^4) + 1)^(1 / 4)

latexify("1 + 0.33 * sqrt(Fl) / (n^0.25) * log10(ReV / 10000)")

latexify(" N_2 / (C / d^2)^2")
latexify("((Q/(Cᵥ*N₁*Fₚ))^2)*(ρ₁/ρ₀)")

# Friction factor
function cheng(Re, ϵ, d)
    a = 1 / (1 + (Re / 2720)^9)
    b = 1 / (1 + (Re / (160 * (d / ϵ)))^2)
    invf = ((Re / 64)^a) * ((1.8 * log10(Re / 6.8))^(2 * (1 - a) * b)) * ((2.0 * log10(3.7 * d / ϵ))^(2 * (1 - a) * (1 - b)))
    return 1 / invf
end

latexify("Δs * cheng(Re, ϵ, d) * (ρ / (2 * d)) * ((1 / A)^2)")

using SymPy, Latexify


μ, ρ, ϵ, d, v, L,q = symbols("μ, ρ, ϵ, d, v, L, q")

simplify(1/sqrt(-1.8*log10(6.9/(ρ*v*d/μ) + (ϵ/d/3.7)^1.11)))

A = π * (d / 2)^2
latexify((ρ * (q / A) * d) / μ)


μ, ρ, ϵ, d, v, L,q = symbols("μ, ρ, ϵ, d, v, L, q")
ω, r, β, A, Q, k, k2 = symbols("ω, r, β, A, Q, k, k_2")

k_2 = 1/(A*tan(β))
cux = ω*r - k2*Q

simplify(ρ * k * (cux)^2)