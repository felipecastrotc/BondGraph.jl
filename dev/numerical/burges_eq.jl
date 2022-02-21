using FFTW, Plots, DifferentialEquations
#  BenchmarkTools

# Examples
# https://github.com/SciML/SciMLBenchmarks.jl
# https://benchmarks.sciml.ai/html/MOLPDE/burgers_spectral_wpd.html

ν = 0.01           # Diffusion coefficient
L = 20              # Domain length
N = 1000            # Number of discretization points

Δx = L / N          # x discretization
x = -L/2:Δx:L/2     # define x domain

# Define discrete wavenumbers
κ = 2 * π * fftfreq(N + 1, 1 / Δx)
# Initial conditions
u0 = 1.0 ./ cosh.(x)
# u0 = cos.(cos.(x .- 0.1))

# Faster FFT
pfft = plan_fft(u0; flags = FFTW.ESTIMATE, timelimit = Inf)
pifft = plan_ifft(pfft * u0; flags = FFTW.ESTIMATE, timelimit = Inf)

function BurgersEq(dudt, u, p, t)
    κ, ν = p        # unpack parameters
    # Spectral approximation of derivatives
    û = pfft * u
    dû = @. 1im * κ * û
    ddû = @. -(κ .^ 2) * û
    # Going back to the point domain
    du = real(pifft * dû)
    ddu = real(pifft * ddû)
    # Burger's equation
    # @. dudt = -u * du
    @. dudt = -u * du + ν * ddu
end

tspan = 100
prob = ODEProblem(BurgersEq, u0, (0.0, tspan), [κ, ν])
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
# sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)

@gif for tᵢ in 0:1:tspan
    # plot(x, h; args...)
    plot(x, sol(tᵢ), title = tᵢ, ylim = (0, 1))
end every 1



