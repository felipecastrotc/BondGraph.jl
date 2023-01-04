using LinearAlgebra, Statistics


# Impeller parameters range
# 20 grau ok entre 

# 200cp
# 1cp
# rodolfo

function γaf(r₁, r₂)
    return r₂^2 - r₁^2
end

function γbf(r₁, r₂, h₁, h₂, β₁, β₂)
    return r₂ / (Ah(2 * r₂, h₂) * tand(β₂)) - r₁ / (Ah(2 * r₁, h₁) * tand(β₁))
end

function Ah(d, h)
    return π * d * h
end

function A(d)
    return (π * (d / 2)^2)
end

function q2v(q, d)
    return q / A(d)
end

function darcyDp(f, ρ, q, d, l)
    return l * f * (ρ / 2) * ((q2v(q, d)^2) / d)
end

function ReQ(ρ, q, d, μ)
    return ρ*q*d/(μ*A(d))
end

function f_lam(ρ, q, d, μ)
    return 64/ReQ(ρ, q, d, μ)
end


rp = Dict()
rp["r₁"] = 60.5 / 2000
rp["r₂"] = 108 / 2000
rp["β₁"] = 16.83
rp["β₂"] = 31.23
rp["h₁"] = 24 / 1000
rp["h₂"] = 13 / 1000

# Unidades -> mm 
#                              P23     P47     P100
# Número de aletas              7       7       7
# Espessura da aleta            3       2       13
# Diâmetro externo              108     108     108
# Espessura do disco            2       3       2
# Diâmetro interno              40      54      60.5
# Altura da aleta na entrada    6       9       24
# Altura da aleta na saída      6       9       13
# Diâmetro do eixo              28.5    28.5    28.5
# Ângulo de saída               46.73   36.45   31.23

# Look file kamal_angles.txt I used librecad to obtain the angles
p23_ang_entr = 90 .- [62.25, 61.95, 65.67, 61.48, 61.99]
mean(p23_ang_entr)
p100_ang_entr = 90 .- [71.7926, 72.9019, 75.3047, 74.8563, 71.3025, 72.8062]
mean(p100_ang_entr)

# Generating a range for impeller variables
# The range I used here are from the Jorge Thesis rounded up and down
# I assumed a range for the external diameter of 150 mm
imp_range = Dict(
    "β₁" => [15, 29],         # degree
    "β₂" => [30, 47],         # degree
    "h₁" => [6, 24] / 1000,    # m
    "h₂" => [6, 13] ./ 1000,    # m
    "r₁" => [40, 61] ./ 2000,     # m
    # "r₁" => [29, 29] ./ 2000,     # m
    "r₂" => [108, 150] ./ 2000,      # m
)

# Geometrical parameters - getting the calculated range
# Cross sectional area
imp_range["A₁"] = 2 * π * imp_range["r₁"] .* imp_range["h₁"]
imp_range["A₂"] = 2 * π * imp_range["r₂"] .* imp_range["h₂"]
# Geometry
γ_a = collect(extrema(imp_range["r₂"]' .^ 2 .- imp_range["r₁"] .^ 2))
# Get the γ_b range
r₁, r₂ = imp_range["r₁"], imp_range["r₂"]
A₁, A₂ = imp_range["A₁"], imp_range["A₂"]
tβ₁, tβ₂ = tand.(imp_range["β₁"]), tand.(imp_range["β₂"])
γ_b = (r₂./(A₂*tβ₂'))[:] .- (r₁./(A₁.*tβ₁'))[:]'
γ_b = collect(extrema(γ_b))

# P100
rp["γ_a"] = γaf(rp["r₁"], rp["r₂"])
rp["γ_b"] = γbf(rp["r₁"], rp["r₂"], rp["h₁"], rp["h₂"], rp["β₁"], rp["β₂"])

# BEP flow rate -> ω*r*tan(β)*A
# Shock loss design flow rate without rpm range
Qₛ = collect(extrema(((r₁.*tβ₁)*A₁')[:]))
rp["Qₛ"] = (rp["r₁"] * tand(rp["β₁"])) * Ah(2*rp["r₁"], rp["h₁"])

# Estimation of the constant of the shock loss, using the results from Biazussi 
# 2014, page 129.
kₛ = [0.126, 0.14] * 4
kₛ = [0.1, 0.2] * 4       # Adjusted for a broader range
rp["kₛ"] = 0.14

# Parameters range
prop_range = Dict(
    # Flow parameters
    "ρ" => [800, 1010],    # kg/m^3
    "Ω" => [0, 1],         # -
)
rp["ρ"] = 1000
rp["ρ"] = 850

# Estimate the Liquid inertia range
ρ = prop_range["ρ"]
# Valores inicialmente aproximados
d = 3 * 25.4 / 1000           # m Diâmetro da linha de tubulação
l = 30                    # m Comprimento da linha
I_l = ρ * 30 / ((pi * d^2) / 4)   # Momento de inércia

rp["d"] = d
rp["I_l"] = rp["ρ"] * l / ((π * rp["d"]^2) / 4)

# Bulk modulus fluids
# Tabela A.3 do White 6ed pt-br pag. 829
b_range = [
    1.31 * 1e9,   # Oil SAE10W
    2.19 * 1e9,   # Water
]
# Compliance 
# It is calculated using the equations 4.36 and 4.39 from
# Karnop 2012 book.
# I added the 100mm diameter to increase slightly the range
d = [d, 0.1]
V₀ = π * ((d ./ 2) .^ 2) * l
k = collect(extrema(V₀ / b_range))
k = [2.7e-11, 8e-11]

# Water and 3" diameter
rp["b"] = b_range[2]
rp["k"] = (π * ((rp["d"] ./ 2) .^ 2) * l/4) / rp["b"]

# Friction range
# From the moody diagram i selected the friction 0.017~0.1
# The lower value was defined to be after the transition if it
# was considered the laminar friction. The system won`t face 
# higher Reynolds values so it should be ok
f = [0.017, 0.1]
Ap = π .* (d ./ 2) .^ 2
fₓ = collect(extrema(ρ' * l ./ (d*(Ap .^ 2)'*2)[:]))

rp["fₓ"] = darcyDp(1, rp["ρ"], 1, rp["d"], l)

# Inertia
# Equation A7 hackel -ρ∫rcot(β(r))dr
# Assuming β constant
Iᵢ = (ρ.*((r₂ .^ 2.0./2)[:]'.-(r₁ .^ 2.0./2)[:])[:]')[:]

Iᵢ = collect(extrema(Iᵢ' .* cotd.(vcat(imp_range["β₁"], imp_range["β₂"]))))
Iᵢ = [0.6, 10.0]

rp["Iᵢ"] = rp["ρ"] * ((rp["r₂"]^2.0 / 2) - (rp["r₁"]^2.0 / 2))

# Shaft
# The material properties was from medium steel alloy
# https://www.matweb.com/search/DataSheet.aspx?MatGUID=f7666326ceb3482f87a9f41ace1d1fb0&ckck=1
ρₐ = [6600, 7860]       # kg/m^3 
lₐ = [0.4, 1]           # m shaft length approximate
dₐ = 0.0285             # m from Biazussi 2014
Iₐ = collect(extrema(((π*(dₐ/2)^2.0.*lₐ)[:]' .* ρₐ) .* (dₐ / 2) .^ 2 / 2))
Iₐ = [1e-4, 6e-4]

rp["Iₐ"] = 5e-4

# Damping
cₐ = [0.01, 1]
rp["cₐ"] = 0.03
# Torque
# The range set was obtained using the following command
# It uses the steady state data mean 
# T = [np.mean(DataMask(i)["ESP_torque"]) for i in sd.paths["process"]]
τ = [0.5, 85]

param_range = Dict(
    "Iᵢ" => Iᵢ,
    "Iₐ" => Iₐ,
    "Iₚᵢ" => I_l,
    "Iₚₒ" => I_l,
    "f" => fₓ,
    "fₛ" => kₛ,
    "Qₛ" => Qₛ,
    "fₚᵢ" => fₓ,
    "fₚₒ" => fₓ,
    "fᵥ" => fₓ,
    "kᵥ" => k,
    "kᵤ" => k,
    "cₐ" => cₐ,
    "γ_a" => γ_a,
    "γ_b" => γ_b,
    "Pₜ" => [0, 0],
    "τ" => τ,
    "ρ" => prop_range["ρ"],
)

# Considering laminar only for now
# (16 * μ * π * d) / (ρ * Q)
# (16*1e-3*π*0.108)/(1000)
# (16*1e-3*π*3*0.0254)/(1000)
# equations(sys)


rp["μ"] = 1e-3
rpo = Dict(
    "Iᵢ" => rp["Iᵢ"],
    "Iₐ" => rp["Iₐ"],
    "Iₚᵢ" => rp["I_l"],
    "Iₚₒ" => rp["I_l"],
    "f" => rp["fₓ"] * f_lam(rp["ρ"], 1.0, 2*rp["r₂"], rp["μ"]),
    "fₛ" => rp["kₛ"],
    "Qₛ" => rp["Qₛ"],
    "fₚᵢ" => rp["fₓ"] * f_lam(rp["ρ"], 1.0, rp["d"], rp["μ"]),
    "fₚₒ" => rp["fₓ"] * f_lam(rp["ρ"], 1.0, rp["d"], rp["μ"]),
    "fᵥ" => rp["fₓ"] * f_lam(rp["ρ"], 1.0, rp["d"], rp["μ"]),
    "kᵥ" => rp["k"],
    "kᵤ" => rp["k"],
    "cₐ" => rp["cₐ"],
    "γ_a" => rp["γ_a"],
    "γ_b" => rp["γ_b"],
    "Pₜ" => 0.0,
    "τ" => 1.0,
    "ρ" => rp["ρ"],
)

rpo