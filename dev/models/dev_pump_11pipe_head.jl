using DifferentialEquations
using ModelingToolkit
using Plots

# ==========================================================
# Utils

function ReEq(ρ, q, d, μ)
    A = π * (d / 2)^2
    return (ρ * (q / A) * d) / μ
end

function darcyq(q, Δs, d, ϵ, ρ, μ)
    Re = ReEq(ρ, q, d, μ)
    # Re = ReEq(ρ, q, d, μ) + 0.01
    # return Δs * cheng(Re, ϵ, d) * (ρ / (2 * d)) * ((1 / A)^2)
    if Re > 0.0
        return Δs * cheng(Re + 0.01, ϵ, d) * (ρ / (2 * d)) * ((1 / A)^2)
    else
        return Δs * cheng(0.01, ϵ, d) * (ρ / (2 * d)) * ((1 / A)^2)
    end
end

function colebrook(Re, ϵ, d)
    # Return the friction factor as 0 for reynolds number equals to 0
    if Re < 2300
        return min(64 / Re, 64)
    end
    # Niazkar approximation
    # https://link.springer.com/article/10.1007/s12205-019-2217-1
    A = -2 * log10((ϵ / d) / 3.7 + 4.5547 / (Re^0.8784))
    B = -2 * log10((ϵ / d) / 3.7 + 2.51 * A / (Re))
    C = -2 * log10((ϵ / d) / 3.7 + 2.51 * B / (Re))
    # Estimate f
    rhs = A - ((B - A)^2) / (C - 2 * B + A)
    f = 1 / (rhs^2)

    # # Set the Niazkar approximation as an initial guess to iterate
    # f₀ = f
    # for i in 1:20
    #     # Colebrook equation
    #     f = (1 / (-2 * log10(ϵ / (3.7 * d) + 2.51 / (Re * sqrt(f₀)))))^2
    #     # Absolute error with tolerance of 1e-7
    #     if abs(f - f₀) < 1e-7
    #         continue
    #     end
    #     f₀ = f
    # end
    return f
end

# Friction factor
function cheng(Re, ϵ, d)
    a = 1 / (1 + (Re / 2720)^9)
    b = 1 / (1 + (Re / (160 * (d / ϵ)))^2)
    invf = ((Re / 64)^a) * ((1.8 * log10(Re / 6.8))^(2 * (1 - a) * b)) * ((2.0 * log10(3.7 * d / ϵ))^(2 * (1 - a) * (1 - b)))
    return 1 / invf
end

# -------------------------------------------------------------------
# Model manual -  1 upstream 1 downstream - Source of effort
function pump_pipe!(du, u, p, t)
    Qi, ω, Q1, Q3, P2, P4 = u[:]
    Δsi, di, Δs, d, ϵ, ρ, μ, γa, γb, Ii, Ia, ca, Kf, Kc, L = p[:]
    # du = zeros(size(u))
    # Impeller
    du[1] = ((P2 - P4)*g*ρa + ρ * (γa * Qi - γb * ω) * ω - darcyq(Qi, Δs, d, ϵ, ρ, μ) * (Qi^2)) / Ii
    # Shaft
    du[2] = (τ(t) - ca * ω - ρ * (γa * Qi - γb * ω) * Qi) / Ia
    # Pipeline
    # Upstream
    # Q1
    du[3] = ((Pi - P2)*g*ρa - Kf * darcyq(Q1, Δs, d, ϵ, ρ, μ) * (Q1^2)) / L
    # Downstream
    # Q3
    du[4] = ((P4 - Po)*g*ρa - Kf * darcyq(Q3, Δs, d, ϵ, ρ, μ) * (Q3^2)) / L
    # Pressures
    du[5] = (Q1 - Qi) / (Kc * C*g*ρa)           # P2 Upstream inlet
    du[6] = (Qi - Q3) / (Kc * C*g*ρa)           # P4 Downstream coupling pump
end

# Properties
g = 9.81
ρa = 1000
# Fluid
# a = 1500                # m/s water
# ρ = 998                 # kg/mˆ3 water
# μ = 1e-3                # Pa*s water
a = 1290                # m/s petroleum
ρ = 850                 # kg/mˆ3 petroleum
μ = 800e-3              # Pa*s petroleum
# Pipeline
d = 75 / 1000             # m Pipe diameter
Δs = 7.5                # m Pipe length
ϵ = 0.1 * 1e-3            # m Steel rugosity
A = π * (d / 2)^2       # m^2 Pipe area
# Pipeline dynamics -> using Tanaka pipe
# https://sci-hub.se/10.1109/IECON.2000.972506
L = ρ * Δs / A
C = A * Δs / (ρ * a^2)
R = 32 * μ / (A * d^2) * Δs
# Impeller
Δsi = Δs       # Equivalent length
di = d      # Impeller's diameter

# State conditions
q = 1e-4
# Reynolds check
Re = (ρ * (q / A) * d) / μ
ReEq(ρ, q, d, μ)

# Friction check
# Check these values, including Reynolds with
# http://www.druckverlust.de/Online-Rechner/dp.php
cheng(Re, ϵ, d)
colebrook(Re, ϵ, d)
Δs * cheng(Re, ϵ, d) * (ρ / (2 * d)) * ((q / A)^2)
darcyq(q, Δs, d, ϵ, ρ, μ) * (q^2)

# P100
γa = 0.002
γb = -1.7317
# γa = 0.2
# γb = -10.7317
Ii = 0.8504
Ia = 5e-4
ca = 0.1
Kf = 2
Kc = 50

# Inputs
τ = (t) -> t > 0.5 ? 80 : 0.01    # Torque
# q = 0                      # flow-rate
# q = 1e-4                  # flow-rate
Po = 0
Pi = 0

# q = 1e-5
q = 0
q0 = q * 1.5
ω0 = 0.0
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
p = [Δsi*40, di, Δs, d, ϵ, ρ, μ, γa, γb, Ii, Ia, ca, Kf, Kc, L]

tspan = (0.0, 2.5)
prob = ODEProblem(pump_pipe!, u0, tspan, p)
sol = solve(prob, reltol=1e-8, abstol=1e-8)

sol
# Q, ω
plot(sol.t, sol[1, :])
plot(sol.t, sol[2, :])
# Q1, Q3
plot(sol.t, sol[3, :])
plot!(sol.t, sol[4, :])
# Q5, Q7
plot(sol.t, sol[5, :])
plot!(sol.t, sol[6, :])
# P2, P4
plot(sol.t, sol[7, :])
plot!(sol.t, sol[8, :])
# P6, P8
plot(sol.t, sol[9, :])
plot!(sol.t, sol[10, :])

t = 0.499:0.0001:1
nsol = sol(t)

Qs = [1, 3, 4]
tl_q = ["Flow rate impeller", "Flow rate pipe 1", "Flow rate pipe 2", "Flow rate pipe 3", "Flow rate pipe 4"]

Ps = [5, 6]
tl_p = ["Pressure drop pipe 1", "Pressure drop pipe 2", "Pressure drop pipe 3", "Pressure drop pipe 3", "Pressure drop pipe 4"]

cfg = Dict(:titlefontsize => 8, :size => (600, 180 * length(Qs)), :legend => false, :xlabel => "Time (s)", :xguidefontsize => 8, :ylabel => "Q (m^3/s)", :yguidefontsize => 8)
plot(t, nsol[Qs, :]', layout=(length(Qs), 1), titles=reshape(tl_q, 1, 5); cfg...)

cfg = Dict(:titlefontsize => 8, :size => (600, 140 * length(Qs)), :legend => false, :xlabel => "Time (s)", :xguidefontsize => 8, :ylabel => "Pressure (Pa)", :yguidefontsize => 8)
plot(t, nsol[Ps, :]', layout=(length(Ps), 1), titles=reshape(tl_p, 1, 5); cfg...)

cfg = Dict(:titlefontsize => 8, :legend => false, :xlabel => "Time (s)", :xguidefontsize => 8, :ylabel => "Rotation (rad/s)", :yguidefontsize => 8)
plot(t, nsol[2, :], titles="Shaft rotation velocity"; cfg...)

using MAT

file = matopen("matfile_head.mat", "w")

t = 0:0.0001:2.0
nsol = sol(t)
plot(nsol)

write(file, "sol", Matrix(nsol))
write(file, "t", collect(t))
close(file)
