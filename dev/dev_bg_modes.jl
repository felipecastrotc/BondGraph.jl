# using BondGraph
using Plots, Symbolics.Latexify
# import BondGraph: t, D
using DifferentialEquations, JLD2

# -----------------------------------------------------------------------------
# Laminar without fluid capacitance

@named Pin = Se(2)

@named C = Spring(k = 1.0)
@named R = Damper(c = 0.5)
@named I = Mass(m = 2)

@named j1 = Junction1(Pin, -I, -R, -C, sgn = -1, couple = false)
@named sys = reducedobs(structural_simplify(j1))
equations(sys)

ts = (0.0, 100)
prob = ODEProblem(sys, [], ts)
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)

tn = ts[1]:0.01:ts[2]
soln = sol(tn)
plot(soln.t, soln[I.f] ./ A)
plot(soln.t, soln[C.q])

# Start
@named Pin = Se(2)
@named up = Junction0(Pin)

# 1
@named C1 = Spring(k = 1.0)
@named R1 = Damper(c = 0.5)
@named I1 = Mass(m = 2)
# @named j11 = Junction1(-I, -R, -C, subsys=[up], sgn = -1)
@named j11 = Junction1(up, -I, -R, -C, sgn = -1)

# 2
@named C2 = Spring(k = 10.0)
@named R2 = Damper(c = 0.5)
@named I2 = Mass(m = 1)
# @named j12 = Junction1(-I, -R, -C, subsys=[up], sgn = -1)
@named j12 = Junction1(up, -I, -R, -C, sgn = -1)

# End
# @named Pout = Se(1)
# @named down = Junction0(-Pout, subsys=[j11, j12], couple=false)
@named down = Junction0(j11, j12, couple = false)
# @named down = Junction0(-Pout, j11, j12, couple=false)

# model2 = ODESystem(equations([down, up, j11, j12]), t; name = :tst)
# equations(model2)
# structural_simplify(model2)

@named sys = reducedobs(structural_simplify(down))
equations(sys)

ts = (0.0, 100)
prob = ODEProblem(sys, [], ts)
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)

tn = ts[1]:0.01:ts[2]
soln = sol(tn)
plot(soln)
plot(soln[j11.C.q])
plot!(soln[j12.C.q])
plot(soln[j11.I.f])
plot!(soln[j12.I.f])

diffeq

u1 = ẋ1
u2 = x1
u3 = ẋ2
u4 = x2

# m*ẍ + k*u[1] + c*u[2] = Pin  - Pout - 

@variables x1(t), x2(t), ẋ1(t), ẋ2(t), F(t), F1(t), F2(t)
@parameters k1 k2 c1 c2 m1 m2

F = 2
eqs = [
    D(ẋ1) ~ -k1 / m1 * x1 - c1 / m1 * ẋ1 + F / m1,
    D(x1) ~ ẋ1,
    D(ẋ2) ~ -k2 / m2 * x2 - c2 / m2 * ẋ2 + F / m2,
    D(x2) ~ ẋ2,
]

sys = ODESystem(eqs, t; name = :tst)
mdl = structural_simplify(sys)

ts = (0.0, 20.0)
# prob = ODEProblem(mdl, [0.0, 0.0, 0.0, 0.0], ts, Dict(k1 => 1.0, k2 => 40.0, c1 => 0.7, c2 => 0.001, m1 => 1.0, m2 => 0.001))
prob = ODEProblem(
    mdl,
    [0.0, 0.0, 0.0, 0.0],
    ts,
    Dict(k1 => 1.0, k2 => 20.0, c1 => 0.4, c2 => 0.5, m1 => 1.0, m2 => 1),
)
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)

plot(sol[t], sol[x1])
plot!(sol[t], sol[x2])

# n = 200
n = length(sol[t])
# plot(sol[t][1:n], sol[x2][1:n])
plot(sol[t][1:n], sol[x2][1:n] .+ sol[x1][1:n])

plot(sol[t], sol[x2])
plot(sol[t], sol[ẋ1])



function diffeq(du, u, p, t)
    # Mode1
    du[1] = p[1] / p[4] - (p[2] / p[4]) * u[2] - (p[3] / p[4]) * u[1]
    du[2] = u[1]
    # Mode2
    du[3] = p[1] / p[7] - (p[5] / (p[7])) * u[4] - ((p[6]) / (p[7])) * u[3]
    du[4] = u[3]
    # Mode3
    du[5] = p[1] / p[10] - (p[8] / (p[10])) * u[6] - ((p[9]) / (p[10])) * u[5]
    du[6] = u[5]
    # Mode4
    du[7] = p[1] / p[13] - (p[11] / (p[13])) * u[8] - ((p[12]) / (p[13])) * u[7]
    du[8] = u[7]
end


p = [2.0, 1.0, 0.7, 1.0, 40.0, 0.001, 0.001, 42.0, 0.01, 0.001, 100.0, 0.00001, 0.0001]
prob = ODEProblem(diffeq, zeros(8), ts, p)
# p = [2.0,   1.0, 0.7, 1.0,   40.0, 0.001, 0.001]
# prob = ODEProblem(diffeq, zeros(4), ts, p)
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)

plot(sol)

u = hcat(sol.u...)
n = length(sol)
n = 1000
plot(sol.t[1:n], sum([u[i, :] for i = 2:2:5], dims = 1)[1][1:n])
# plot(sum([u[i, :] for i in 1:2:8], dims=1))


# -----------------------------------------------------------------------------
# Modes fit
