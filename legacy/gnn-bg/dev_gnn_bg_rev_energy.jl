using CSV
using LinearAlgebra
using Plots
using LaTeXStrings
using ScikitLearn
using Statistics
using DataFrames
using Serialization
using Random
using Symbolics
using LsqFit
using DifferentialEquations


function msd!(du, u, p, t)
    A = [0 1; -k/m -b/m];
    tmp = A*u;
    du[1] = tmp[1];
    du[2] = tmp[2];
end

k = 20;
m = 0.5;
b = 0.5;

u = [1, 0];
tₛ = (0.0, 100.0);

du = [0.0, 0.0]
prob = ODEProblem(msd!, u, tₛ);
sol = solve(prob);

Δt = 0.01
# t = 0:Δt:20;
t = 0:Δt:4;
X = sol(t);
plot(t, X')

q = X[1, :];
f = X[2, :];
f =  vcat([0], diff(X[1, :])/Δt);

plot(f)

# e = k₃.*(q.^3);
e = k.*q;
p = m.*f;

Eₖ = 0.5.*q.*e;
Eₚ = 0.5.*f.*p;
plot(t, Eₖ)
plot!(t, Eₚ)
plot!(t, cumsum((b.*f).*f.*Δt))

plot(t, cumsum((b.*f).*f.*Δt))
plot!(t, (Eₖ + Eₚ)[1] .- (Eₖ + Eₚ))

plot(cumsum((b.*f).*f.*Δt) .+ (Eₖ + Eₚ))
plot(t, Eₖ + Eₚ)


# Analytical

k = 10;
m = 0.2;
b = 0.45;

A, β, ω = 1, -(b/(2*m)), sqrt(k/m);
x(t) = A*exp.(β.*t).*cos.(ω*t)
dx(t) = A*(β*exp.(β*t).*cos.(ω*t) .- ω*exp.(β*t).*sin.(ω*t))

plot(t, x(t))
plot!(t, dx(t))

q = x(t);
f = dx(t)
e = k.*q;
p = m.*f;

Eₖ = 0.5.*q.*e;
Eₚ = 0.5.*f.*p;
plot(t, Eₖ)
plot!(t, Eₚ)

plot(t, Eₖ + Eₚ )
plot(t, cumsum((b.*f).*f.*Δt))

plot(t, )


plot(b*f)
plot!(k*q)
plot(m*diff(f)/Δt)

plot((k*q)[2:end] .+ (b*f)[2:end] .+ m*diff(f)/Δt)


M()

plot(k*q.*f)
plot!((m*diff(f)/Δt).*f[1:end-1])
plot!((b*f).*f)




