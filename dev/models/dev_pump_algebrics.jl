using SymPy

ρ, μ, d, q, ω,ΔP = symbols("ρ, μ, d, q, ω, ΔP")

k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11 = symbols("k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11")


# Jorge thesis 2014
Ẇ_euler = 
string(Ẇ_euler)

# Eq 3.22
Ẇ_mecanicas = d^3*ΔP*ω*k7
string(Ẇ_mecanicas)

# Eq 3.23
Ẇ_disco = k8*μ*ω^2*d^3 + k9*ρ*ω^3*d^5
string(Ẇ_disco)

# Eq 3.25
Ẇ_arrasto = d^5*ρ*ω^3*(1/2 - 2*q*k1/(d^3*ω))^2*k10*(μ/(d^2*ρ*ω) + (1/2 - 2*q*k1/(d^3*ω))*k11)
Ẇ_arrasto = simplify(Ẇ_arrasto)
string(Ẇ_arrasto)


