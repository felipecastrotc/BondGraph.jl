using FFTW, Plots, DifferentialEquations #  BenchmarkTools

# How to couple PDEs
# https://www.youtube.com/watch?v=kiC7nLZsEwA

function BurgersModEq(dudt, u, p, t)
    κ, ν, ρ, dP = p     # unpack parameters
    # Spectral approximation of derivatives
    û = pfft * u
    dû = @. 1im * κ * û
    ddû = @. -(κ .^ 2) * û
    # Going back to the point domain
    du = real(pifft * dû)
    ddu = real(pifft * ddû)
    # Burger's equation with pressure
    @. dudt = -u * du + ν * ddu + (1 / ρ) * dP
end

function PressureEq(dPdt, P, p, t)
    ρ, a, du, u = p     # unpack parameters
    # Spectral approximation of derivatives
    P̂ = pfft * P
    dP̂ = @. 1im * κ * P̂
    # Going back to the point domain
    dP = real(pifft * dP̂)
    # Pressure equation
    @. dPdt = -u * dP + (ρ * a^2) * du
end

ν = 0.01           # Diffusion coefficient
ρ = 1000           # Specific mass
a = 1000           # Wave speed
L = 20             # Domain length
N = 1000           # Number of discretization points

Δx = L / N          # x discretization
x = -L/2:Δx:L/2     # define x domain

# Define discrete wavenumbers
κ = 2 * π * fftfreq(N + 1, 1 / Δx)
# Initial conditions
u0 = zeros(length(x))
P0 = zeros(length(x))

# Faster FFT
pfft = plan_fft(u0; flags = FFTW.ESTIMATE, timelimit = Inf)
pifft = plan_ifft(pfft * u0; flags = FFTW.ESTIMATE, timelimit = Inf)



tspan = 100
prob = ODEProblem(BurgersEq, u0, (0.0, tspan), [κ, ν])
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
# sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)

@gif for tᵢ in 0:1:tspan
    # plot(x, h; args...)
    plot(x, sol(tᵢ), title = tᵢ, ylim = (0, 1))
end every 1



