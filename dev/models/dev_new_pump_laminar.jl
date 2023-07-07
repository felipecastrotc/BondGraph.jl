using DifferentialEquations
using Plots

# ==========================================================
# Utils

function AEq(d)
    return π * (d / 2)^2
end

function ReEq(ρ, q, d, μ)
    A = π * (d / 2)^2
    return (ρ * (q / A) * d) / μ
end

function darcyq(Δs, d, μ)
    # using SymPy
    # ρ, μ, q, d, l, f, A = symbols("ρ, μ, q, d, l, f, A")

    # Ad = A
    # Ad = π*(d/2)^2
    # Re = (ρ * (q/Ad) * d / μ)
    # fd = f
    # fd = 64/Re
    # l*fd*(ρ/(2*d))*(q/Ad)^2
    # string(l * fd * (ρ / (2 * d)) * (q / Ad)^2)
    return 128 * Δs * μ / (pi * d^4)
end

# -------------------------------------------------------------------
# Model manual -  1 upstream 1 downstream - Source of effort
function pump_pipe!(du, u, p, t)
    Pi = 0.0
    Po = 0.0
    Qi, ω, Q1, Q3, P2, P4 = u[:]
    Δsi, di, Δs, d, ϵ, ρ, μ, γa, γb, Ii, Ia, ca, Kf, Kc, L = p[:]
    # du = zeros(size(u))
    # Impeller
    du[1] = (P2 - P4 + ρ * (γa * Qi - γb * ω) * ω - darcyq(Δs, d, μ)*Qi) / Ii
    # Shaft
    du[2] = (τ(t) - ca * ω - ρ * (γa * Qi - γb * ω) * Qi) / Ia
    # Pipeline
    # Upstream
    # Q1
    du[3] = (Pi - P2 - Kf * darcyq(Δs, d, μ)*Q1) / L
    # Downstream
    # Q3
    du[4] = (P4 - Po - Kf * darcyq(Δs, d, μ)*Q3) / L
    # Pressures
    du[5] = (Q1 - Qi) / (Kc * C)           # P2 Upstream inlet
    du[6] = (Qi - Q3) / (Kc * C)           # P4 Downstream coupling pump
end

# Properties
# Fluid
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
# Impeller
Δsi = 0.2       # Equivalent length
di = 0.108      # Impeller's diameter

# P100
γa = 0.002
γb = -1.7317
Ii = 0.8504
Ia = 5e-4
ca = 0.1
Kf = 1
Kc = 1

0.7-0.25
# Inputs
# τ = (t) -> (tanh.(50 .*(t .- 0.31)) .+ 1)*(80/2)
τ = (t) -> (tanh.(100 .*(t .- 0.275)) .+ 1.0001)*(80/2)
τ = (t) -> (tanh.(100 .*(t .- 0.275)) .+ 1.0001)*(80/2)
τ = (t) -> t > 0.25 ? min(80, 10000*(t - 0.25)) : 1e-4    # Torque
τ = (t) -> t > 0.25 ? 80 : 1e-3    # Torque
τ = (t) -> (tanh.(200 .*(t .- 0.265)) .+ 1.0001)*(80/2)

# t = collect(0.23:0.0001:0.3)
# plot(t, τ.(t))
# plot!(t, 40 .*(tanh.(500 .*(t .- 0.2525)) .+ 1))
# t = collect(0.24:0.0001:0.27)
# plot(t, 40 .*(tanh.(250 .*(t .- 0.26)) .+ 1))

#      Qi,   ω,  Q1,  Q3,  P2, P4
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
p = [Δsi * 40, di, Δs, d, ϵ, ρ, μ, γa, γb, Ii, Ia, ca, Kf, Kc, L]

tspan = (0.0, 1.0)
prob = ODEProblem(pump_pipe!, u0, tspan, p)
sol = solve(prob, reltol=1e-8, abstol=1e-8)

plot(sol.t, sol[1, :])
plot!(twinx(), sol.t, τ.(sol.t), c=:red)

# Q, ω
plot(sol.t, sol[1, :])
plot(sol.t, sol[2, :])
plot!(twinx(), sol.t, τ.(sol.t), c=:red)
# Q1, Q3
plot(sol.t, sol[3, :])
plot!(twinx(), sol.t, τ.(sol.t), c=:red)
plot!(sol.t, sol[4, :])
# Q5, Q7
plot(sol.t, sol[5, :])
plot(sol.t, sol[6, :])
plot!(twinx(), sol.t, τ.(sol.t), c=:red)
# P2, P4
plot(sol.t, sol[7, :])
plot!(sol.t, sol[8, :])
# P6, P8
plot(sol.t, sol[9, :])
plot!(sol.t, sol[10, :])

t = 0:0.0001:1
nsol = sol(t)

idxs = Int(0.24 / 0.0001):Int(0.255 / 0.0001)
plot(nsol.t[idxs], nsol[2, idxs])
plot!(twinx(), nsol.t[idxs], τ.(nsol.t[idxs]), c=:red)

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

file = matopen("matfile.mat", "w")

t = 0:0.0001:2.0
nsol = sol(t)
plot(nsol)

write(file, "sol", Matrix(nsol))
write(file, "t", collect(t))
close(file)


64 / (ρ * q * d / μ) * q^2