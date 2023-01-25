# Support functions 

function AEq(d)
    return π * (d / 2)^2
end

# -----------------------------------------------------------------------------
# Flow functions

function ReQ(ρ, q, d, μ)
    return (4*ρ*q)/(μ*π*d)
    # A = π * (d / 2)^2
    # return (ρ * (q / A) * d) / μ
end

function darcyq_haal(q, Δs, d, ϵ, μ, ρ)
    return ((2*Δs*ρ)/(π*d^3))*haaland(q, d, ϵ, μ, ρ)
    # return (4*Δs*ρ)/(2*d*π*d^2)
end

function darcyq_lam(Δs, d, μ)
    # using SymPy
    return 128 * Δs * μ / (pi * d^4)
end

# -----------------------------------------------------------------------------
# Friction functions
function haaland(q, d, ϵ, μ, ρ)
    return 1/((-1.8*log10(6.9/ReQ(ρ, q, d, μ) + ((ϵ/d)/3.7)^1.11))^2)
end

# -----------------------------------------------------------------------------
# Valve

# Valve Cv function of opening
function cv_func(a, p)
    return p[1] + p[2] * a + p[3] * a^2 + p[4] * a^3 + p[5] * a^4 + p[6] * a^5
end

# -----------------------------------------------------------------------------
# ISA-75 equations

# Water specific mass @ 15C 
ρ₀ = 999.1026214670904

# For laminar flow where reynolds <= 10000 the ISA-75 standard, it has a 
# correction factor, Fr.

# Eq. 1 -> Get the flow-rate given a pressure difference
function QValve(ΔP, ρ₁, Cᵥ, N₁, Fₚ=1.0)
    return Cᵥ * N₁ * Fₚ * sqrt(ΔP / (ρ₁ / ρ₀))
end

# Eq. 1 -> Get the pressure difference for a given flow-rate
function ΔPValve(Q, ρ₁, Cᵥ, N₁, Fₚ=1.0)
    return ((Q/(Cᵥ*N₁*Fₚ))^2)*(ρ₁/ρ₀)
end

# Eq. 23 -> Get the valve Reynolds number
function ReValve(ρ, q, d, μ, C, N_2, N_4, Fd, Fl)
    return ((ρ * N_4 * Fd * q) / (μ * sqrt(C * Fl))) * ((Fl^2 * C^2) / (N_2 * d^4) + 1)^(1 / 4)
end

# Eq A.6 for ReValve  >= 10
function Frl10(n, ReV, Fl)
    return min((0.026 / Fl) * sqrt(n * ReV), 1.0)
end

# Eq A.7 for ReValve  >= 10
function Frg10(n, ReV, Fl)
    return min(min(1 + 0.33 * sqrt(Fl) / (n^0.25) * log10(ReV / 10000), Frl10(n, ReV, Fl)), 1.0)
end

# Eq. A.8a
function n_fulltrim(C, d, N_2)
    return N_2 / (C / d^2)^2
end

# Eq. A.8b
function n_redtrim(C, d, N_32)
    return 1.0 + N_32*(C/d^2)^(2/3)
end