using SymPy

ρ, μ, d, q, ω,ΔP = symbols("ρ, μ, d, q, ω, ΔP")

k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11 = symbols("k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11")

n = symbols("n")

# ==============================================================================
# Jorge thesis 2014
# ==============================================================================

# ------------------------------------------------------------------------------
# ΔP Equation

# -----------------------
# Friction losses

# Eq 6.12 or 3.11
ΔPatrito = k2*(k3*μ*d/(ρ*q) + k4*(μ*d/(ρ*q))^n)*(ρ*ω^2*d^2*(q/(ω*d^3))^2)
expand(simplify(ΔPatrito))
# Organizing
ΔPatrito_o = k2*k3*q*μ/d^3 + k2*k4*q^(2-n)*ρ^(1 - n)*d^(n - 4)*μ^n
ΔPatrito_o

# Testing if they are equivalents
t = []
for i in 1:10000
    v = rand(8)
    push!(t, subs(ΔPatrito, k2 => v[1], k3 =>v[2], k4=>v[3], ρ=> v[4], μ=> v[5], d => v[6], q => v[7], n => v[8]) - subs(ΔPatrito_o, k2 => v[1], k3 =>v[2], k4=>v[3], ρ=> v[4], μ=> v[5], d => v[6], q => v[7], n => v[8]))
end
sum(abs.(t) .> 1e-4)

# Generating model equations
# Eq 3.11
ΔPatrito311 = (q^2*ρ/d^4)*((d*μ*k2/(q*ρ)) + k3*((d*μ/(q*ρ))^n))
expand(simplify(ΔPatrito311))
string(expand(simplify(ΔPatrito311)))
ΔPₐ₁ = k2*q*μ/d^3
ΔPₐ₂ = k3*q^2*ρ*(d*μ/(q*ρ))^n/d^4
ΔPₐ₂ = k3*q^2*ρ*d^n*μ^n/(d^4*q^n*ρ^n)
ΔPₐ₂ = k3*q^(2-n)*ρ^(1-n)*d^(n-4)*μ^n

# -----------------------
# Vortex losses

# Eq 6.15 analysis
r1, r2, b1, b2, β1, β2, π1 = symbols("r_1 r_2 b_1 b_2 β_1 β_2 π")

q1 = 2*π1*b1*ω*r1^2*tan(β1)
q2 = 2*π1*b2*ω*r2^2*tan(β2)

ΔPchoques = ρ*ω^2*r1^2*(1 - q/q1)/2 + ρ*ω^2*r2^2*(1 - q/q2)/2
expand(ΔPchoques)
ΔPchoques = 

ρ*ω^2
(r1^2 - r1^2*q/q1) + (r2^2 - r2^2*q/q2)
ΔPchoques = ρ*ω^2*r1^2*(1 - q/q1) + ρ*ω^2*r2^2*(1 - q/q2)

d^2 + d^2*r2^2/r1^2

# Generating model equations
# Eq. 3.15
ΔPturb = ρ*ω^2*d^2*k4*(1 - k5*q/(ω*d^3))^2
expand(simplify(ΔPturb))
string(expand(simplify(ΔPturb)))
ΔPₜ₁ = d^2*k4*ρ*ω^2
ΔPₜ₂ = - 2*k4*k5*q*ρ*ω/d
ΔPₜ₃ = k4*k5^2*q^2*ρ/d^4

# -----------------------
# Localized loss 

# Generating model equations
# Eq. 3.16
ΔPloc = k6*ρ*ω^2*d^2*(q/(ω*d^3))^2
string(simplify(ΔPloc))


# -----------------------
# Model comparison

# Comparing 6.18 and 3.17
# Third term -> The same
eq618 = k2*(k3*(μ*d/(ρ*q)) + k4*(μ*d/(ρ*q))^n)*ρ*ω^2*d^2*(q/(ω*d^3))^2
eq317 = q^2*ρ*(d*μ*k2/(q*ρ) + k3*(μ*d/(ρ*q))^n)/d^4
# Fourth term -> The same


# ------------------------------------------------------------------------------
# Ẇ Equation

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


