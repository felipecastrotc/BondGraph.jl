using SymPy
using MAT

file = matopen("matfile.mat")
keys(file)

file["sol"]
file["sol"]

# plt.plot(sol[:, 0]*100)
# plt.plot(sol[:, 1]/10)
# plt.plot(sol[:, 2]*100)
# plt.plot(sol[:, 3]*100)
# plt.plot(sol[:, 4]/500000)
# plt.plot(sol[:, 5]/500000)
# plt.plot(sol[:, 5]/5e5)
# plt.plot(sol[:, 4])
# Qi, ω, Q1, Q3, P2, P4, dQi, dω, dQ1, dQ3, dP2, dP4 = u

ρ, μ, d, ϵ, Δs, γa, γb, Ii, K = symbols("ρ, μ, d, ϵ, Δs, γa, γb, Ii, K")
P2, P4, Qi, ω = symbols("P2, P4, Qi, ω")


500000 - 5e5
5000000000 - 5e5 * 1e4
# du[1] = 

e = (P2 * 5e5 - P4 * 5e5 + ρ * (γa * Qi / 100 - γb * ω * 10) * ω * 10 - K * ((Qi / 100)^2)) / Ii

string(simplify(e * 100))
simplify(e * 100)

darcyq(Qi, Δs, d, ϵ, ρ, μ)

m, s, kg, h = symbols("m, s, kg, h")


a = m/s^2
N = kg*a
Pa = N/(m^2)
Q = m^3/s
Qh = m^3/h
ρ = kg/m^3
g = a
cvt = 1/(g*ρ)
H = Pa*cvt

1/9000
I = Pa/(Q/s)
H*(g*ρ)/(I)

# I/

I/((g*ρ)*(3600*s/h))

9.81*1000*3600

Ih = H/(Qh/s)

Ih/cvt

