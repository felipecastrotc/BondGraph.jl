using XLSX, DataFrames, Polynomials
using Plots, LsqFit, YAML
using SymPy, Latexify

# Path used to store valve info
PATH_VALVE = "/Users/fctc/me/doutorado/dev/pump/valve-docs/"

# Valve data
xlsx = XLSX.readxlsx(PATH_VALVE * "info_data.xlsx")
vlv = DataFrame(xlsx["data"][:][:, 2:end]', xlsx["data"][:][:, 1])
vlv = float.(vlv)

# C_v fit
col = "Cv"
pol = fit(vlv[!, "prct"], vlv[!, col], 5)

# Graph analysis of Cv polynomial fit
x = 0.1:0.001:1
plot(x, pol.(x))
plot!(vlv[!, "prct"], vlv[!, col])

# Store info
yml_file = "./sims/cfg/valve_params.yml"
ps = YAML.load_file(yml_file)
ps = Dict()
ps[:C_v] = pol.coeffs
ps[:C_v_func] = "cv_func"
f = open(yml_file, "w")
YAML.write(f, ps)
close(f)

# Fd fit
col = "Fd"
pol = fit(vlv[!, "prct"], vlv[!, col], 5)

# Graph analysis of Cv polynomial fit
x = 0.1:0.001:1
plot(x, pol.(x))
plot!(vlv[!, "prct"], vlv[!, col])

# Store info
ps = YAML.load_file(yml_file)
ps[:F_d] = pol.coeffs
ps[:F_d_func] = "fd_func"
f = open(yml_file, "w")
YAML.write(f, ps)
close(f)

# Store other parameters from valve
N₁ = 8.65 * 1e-2
N₂ = 2.14 * 1e-3
N₄ = 7.6 * 1e-2
Fl = 0.82
d = 3 * 25.4

ps = YAML.load_file(yml_file)
ps[Symbol(:N_1)] = N₁
ps[Symbol(:N_2)] = N₂
ps[Symbol(:N_4)] = N₄
ps[Symbol(:F_l)] = Fl
ps[Symbol(:d)] = d
f = open(yml_file, "w")
YAML.write(f, ps)
close(f)

# Params to calculate the Reynolds of the valve
Fd = 0.32
Q = 63.75
d = 3 * 25.4
C = 4.32
μ = 0.189
ρ = 850

Re = ((ρ * N₄ * Fd * Q) / (μ * sqrt(C * Fl))) * ((Fl^2 * C^2) / (N₂ * d^4) + 1)^(1 / 4)

# Get some equations
μ, ρ, ϵ, d, v, L, q = symbols("μ, ρ, ϵ, d, v, L, q")
ΔP, ρ₁, ρ₀, Fₚ, N₁, C = symbols("ΔP, ρ₁, ρ₀, Fₚ, N₁, C")

# Get valve equation
# Eq. 1 from IS75.01.01
isa = q ~ C*N₁*Fₚ*sqrt(ΔP/(ρ₁/ρ₀))
p = (C * N₁ * Fₚ)
isa = isa.lhs / p ~ isa.rhs/p
isa = isa.lhs^2 ~ isa.rhs^2
p = (ρ₁ / ρ₀)
isa = isa.rhs*p ~ isa.lhs*p
isa.rhs == ΔPValve(q, ρ₁, C, N₁, Fₚ)

# Also checked with the tool below and I got the same result
# https://www.swagelok.com/en/toolbox/cv-calculator

# Other equations
latexify("1 + 0.33 * sqrt(Fl) / (n^0.25) * log10(ReV / 10000)")
latexify(" N_2 / (C / d^2)^2")
latexify("((Q/(Cᵥ*N₁*Fₚ))^2)*(ρ₁/ρ₀)")
latexify("Δs * cheng(Re, ϵ, d) * (ρ / (2 * d)) * ((1 / A)^2)")

simplify(1 / sqrt(-1.8 * log10(6.9 / (ρ * v * d / μ) + (ϵ / d / 3.7)^1.11)))

A = π * (d / 2)^2
latexify((ρ * (q / A) * d) / μ)


μ, ρ, ϵ, d, v, L, q = symbols("μ, ρ, ϵ, d, v, L, q")
ω, r, β, A, Q, k, k2 = symbols("ω, r, β, A, Q, k, k_2")

k_2 = 1 / (A * tan(β))
cux = ω * r - k2 * Q

simplify(ρ * k * (cux)^2)