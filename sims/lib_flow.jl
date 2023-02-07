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

function darcyq_haal(q, l, d, ϵ, μ, ρ)
    return haaland(q, d, ϵ, μ, ρ)*(8*l*ρ/(π^2*d^5))
    # return (4*l*ρ)/(2*d*π*d^2)
end

function darcyq_cheng(q, l, d, ϵ, μ, ρ)
    Re = max(ReQ(ρ, abs(q), d, μ), 0.1)
    return cheng(Re, ϵ, d)*(8*l*ρ/(π^2*d^5))
end

function darcyq_lam(l, d, μ)
    # using SymPy
    return 128 * l * μ / (pi * d^4)
end

# -----------------------------------------------------------------------------
# Friction functions
function haaland(q, d, ϵ, μ, ρ)
    Re = max(ReQ(ρ, abs(q), d, μ), 0.1)
    if Re < 2400
        return 64/Re
    else
        return 1/((-1.8*log10(6.9/Re + ((ϵ/d)/3.7)^1.11))^2)
    end
end

function cheng(Re, ϵ, d)
    a = 1 / (1 + (Re / 2720)^9)
    b = 1 / (1 + (Re / (160 * (d / ϵ)))^2)
    invf = ((Re / 64)^a) * ((1.8 * log10(Re / 6.8))^(2 * (1 - a) * b)) * ((2.0 * log10(3.7 * d / ϵ))^(2 * (1 - a) * (1 - b)))
    return 1 / invf
end

# -----------------------------------------------------------------------------
# Valve

# Valve Cv function of opening
function cv_func(a, p)
    return p[1] + p[2] * a + p[3] * a^2 + p[4] * a^3 + p[5] * a^4 + p[6] * a^5
end

# Valve Cv function of opening
function fd_func(a, p)
    return p[1] + p[2] * a + p[3] * a^2 + p[4] * a^3 + p[5] * a^4 + p[6] * a^5
end

# -----------------------------------------------------------------------------
# Progressing cavity pump
function pcavity(Q, μ, ω, k1, k2, P_in=0)
    return μ*(k1*ω - Q)/k2 + P_in
end

function ptwin(Q, μ, ω, k1, k2, P_in=0)
    # return k2*μ*(k1*ω - Q) + P_in
    return k1*μ*ω - k2*μ*Q + P_in
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