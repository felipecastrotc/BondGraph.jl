using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D
using DifferentialEquations, JLD2


function colebrook(ϵ, d, Re)
    # Return the friction factor as 0 for reynolds number equals to 0
    if Re < 2300
        return min(64/Re, 64)
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

function calcRe(V, ρ, μ, d)
    return ρ .* V .* d ./ μ
end

# -----------------------------------------------------------------------------
# Unit checking
@variables m, kg, s

uA = m^2
uV = m^3
ua = m / s^2
uN = kg * ua
uPa = uN / uA
uμ = uPa * s
uQ = m^3 / s
uρ = kg / uV

m * uρ * ua

uPa/(uρ*m/uA)


# -----------------------------------------------------------------------------
# Bond-graph linear 
g = 9.81        # m/s^2 - Gravity

l = 1.0         # m - Pipe segment length
d = 0.01        # m - Pipe segment diameter
A = π * d^2 / 4     # m^2 - Pipe section area
μ = 1e-3        # Pa*s - Viscosity
ρ = 998         # kg/m³ - Water
B = 2.2 * 1e9       # Pa -> Fluid bulk modulus

# H = P / (ρ * g)
Hin = 0.01        # m - Pressure head in
Hout = 1    # m - Pressure head OUT

@named Pin = Se(Hin * ρ * g)
@named Pout = Se(Hout * ρ * g)

Rᵥ = 128 * μ * l / (π * d^4)
Cᵥ = A * l / B
Iᵥ = ρ * l / A

@named C = Spring(k = Cᵥ)
@named R = Damper(c = Rᵥ)
@named I = Mass(m = Iᵥ)

1/Cᵥ
Δh = 0.0001
@named Pin = Se(Δh * ρ * g)
@named j1 = Junction1(Pin, -I, -R, sgn = -1)
@named j0 = Junction0(j1, -C, couple = false)
equations(j0)
equations(j1)
@named sys = reducedobs(structural_simplify(j0))
equations(sys)

ts = (0.0, 10.0)
prob = ODEProblem(sys, [], ts)
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)

plot(sol.t, sol[I.f]./ A)

tr = ts[1]:(diff(collect(ts))/100)[1]:ts[2]
s = sol(tr)

plot(s.t, s[1, :] ./ A)
plot!(s.t, s[j1.I.f] ./ A)

# plot(sol)
(Hin - Hout) * ρ * g

# (s[j1.I.f]./A)[end]
(sol[I.f]./A)[end]

1/A

v = 27
ϵ = 1e-5            # m Roughness
for k in 1:20
    Re = calcRe(v, ρ, μ, d)
    f = colebrook.(ϵ, d, Re)
    v = sqrt(Δh * 2 * g * d / (l * f))
end
v

# -----------------------------------------------------------------------------
# Laminar without fluid capacitance

g = 9.81        # m/s^2 - Gravity

l = 0.1         # m - Pipe segment length
d = 0.01        # m - Pipe segment diameter
A = π * d^2 / 4     # m^2 - Pipe section area
μ = 1e-3        # Pa*s - Viscosity
ρ = 998         # kg/m³ - Water
B = 2.2 * 1e9       # Pa -> Fluid bulk modulus

# Pressure drop 0.005m/m
Δh = 0.05*l     # m - Pressure head in
Δh = 0.005*l     # m - Pressure head in
Δh = 0.003*l     # m - Pressure head in
Δh = 0.5
@named Pin = Se(Δh * ρ * g)

Rᵥ = 128 * μ * l / (π * d^4)
Iᵥ = ρ * l / A

@named R = Damper(c = Rᵥ)
@named I = Mass(m = Iᵥ)

@named j1 = Junction1(Pin, -I, -R, sgn = -1, couple=false)
equations(j1)
@named sys = reducedobs(structural_simplify(j1))
equations(sys)

ts = (0.0, 0.3)
ts = (0.0, 20)
prob = ODEProblem(sys, [], ts)
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)

tr = ts[1]:(diff(collect(ts))/100)[1]:ts[2]
s = sol(tr)

plot(s.t, s[1, :] ./ A)

v

v = 0.1
ϵ = 1e-5            # m Roughness
for k in 1:20
    Re = calcRe(v, ρ, μ, d)
    f = colebrook.(ϵ, d, Re)
    v = sqrt(Δh * 2 * g * d / (l * f))
end
v - (s[1, :] ./ A)[end]

f = jldopen("../numerical/data/sim_laminar.jld2")
plot(s.t, s[1, :] ./ A)
plot!(f["t"], f["V"])

# -----------------------------------------------------------------------------
# Multiple

ir = Junction1(-I, -R, sgn = -1, name = :ir)

M = [Junction1(Pin, -I, -R, sgn = -1, name = :s11)]
push!(M, Junction0(-C, subsys = M[1], name = :s10))
push!(M, Junction1(ir, subsys = M[2], name = :s21))
push!(M, Junction0(-C, subsys = M[3], name = :s20))
push!(M, Junction1(ir, subsys = M[4], name = :s31))
push!(M, Junction0(-C, subsys = M[5], name = :s30))
push!(M, Junction1(ir, subsys = M[6], name = :s41))
push!(M, Junction0(-C, subsys = M[7], name = :s40, couple = false))

model = ODESystem(equations(M), t; name = :tst)

@named sys = reducedobs(structural_simplify(model))
equations(sys)

prob = ODEProblem(sys, [], (0.0, 4.0))
sol = solve(prob)
