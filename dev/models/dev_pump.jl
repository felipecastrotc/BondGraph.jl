using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D
using DifferentialEquations, JLD2


# -----------------------------------------------------------------------------
# SIMPLER

@parameters ρ, T, γa, γb, Is, If
@variables ω(t), Q(t)

# Fluid
# dᵥ = 0.137    # m - impeller external diameter
ρᵥ = 1000;
d1ᵥ = 4*2*0.0254;       # in - impeller internal diameter
d2ᵥ = 7*2*0.0254;       # in - impeller external diameter
r1ᵥ = d1ᵥ/2;      # m - impeller internal radiu
r2ᵥ = d2ᵥ/2;      # m - impeller external radiu
β1 = 30;       # degree - impeller internal angle
β2 = 20;       # degree - impeller external angle
b1 = 1.75*0.0254;     # in  - internal constant
b2 = 1.75*0.0254;     # in  - external constant
γaᵥ = ρᵥ*(r2ᵥ^2 - r1ᵥ^2);
γbᵥ = (ρᵥ/(2*π))*(cotd(β2)/b2 - cotd(β1)/b1);

gyp = γa*ω - γb*Q
gyp = γa*ω - γb*Q

# Shaft
# T = 10
@named sT = Se(T)
@named sI = Mass(m=1.0)
@named sR1 = Damper(c=10.0)
@named sJ1 = Junction1(sT, -sI, -sR1, sgn=-1)

# Impeller
@parameters RL, RT

@named iI = Mass(m=1.0)
@named iRL = Damper(c = 10)
@named iRT = Damper(c = RT)
@named iJ1 = Junction1(-iI, -iRL, -iRT, sgn=-1)

# mGY Connection
@named gy = mGY(sJ1, iJ1, g = gyp)

equality = [sJ1.f ~ ω, iJ1.f ~ Q]
@named sys = ODESystem(vcat(equations(gy), equality)) 
@named mdl = reducedobs(structural_simplify(sys))

equations(mdl)

# human = Dict(sJ1.sI.f => ω, iJ1.iI.f => Q, sJ1.sI.I => Is, iJ1.iI.I => If, sJ1.sT.T => T, iJ1.iRT.RT => RT, iJ1.iRL.R => RL)
human = Dict(sJ1.sI.f => ω, iJ1.iI.f => Q, sJ1.sI.I => Is, iJ1.iI.I => If, sJ1.sT.T => T, iJ1.iRT.R => RT*Q, iJ1.iRL.R => RL)
mdl = renamevars(mdl, human)

equations(mdl)

# [12, 60]/8

A = (π*((0.05/2)^2))
Ifᵥ = ρᵥ*10/A

τᵥ = -2;

# vals = [T => τᵥ, ρ => ρᵥ, γa => γaᵥ, γb => γbᵥ, sJ1.sR1.R => 0.1, Is => 0.02, If => 0.0005, RT => 1000, RL => 0]
vals = [T => τᵥ, ρ => ρᵥ, γa => γaᵥ, γb => γbᵥ, sJ1.sR1.R => 1, Is => 50, If => Ifᵥ, RT => 1e5, RL => 0]
prob = ODEProblem(mdl, [0.0, 0.0], (0.0, 10000.0), vals)

equations(mdl)
sol = solve(prob, reltol=1e-8, abstol=1e-8)
plot(sol)

plot(sol.t, -sol[2, :]/A)
plot(sol.t, -sol[1, :])

ωᵥ, Qᵥ = -7.698944611469807, -0.003926997971959516
(-1e5*Qᵥ - 1e5*(Qᵥ^2) - γaᵥ*ωᵥ*ωᵥ + γbᵥ*Qᵥ*ωᵥ) / Ifᵥ

# -----------------------------------------------------------------------------
# SIMPLER + pipe

@parameters ρ, μ, T, L, γa, γb, Is, If
@variables ω(t), Q(t)

# Properties
ρᵥ = 1000;         # kg/m^3 - Water
μᵥ = 1e-3;         # Pa*s   - Water
Lᵥ = 0.1;          # m - Pipe length
dᵥ = 0.05;         # m - Pipe diameter
Aᵥ = π*(dᵥ^2)/4;   # m^2 - Pipe area
# d1ᵥ = 0.0605;    # m - P100 internal diameter
# d2ᵥ = 0.108;     # m - P100 external diameter
d1ᵥ = 0.07;        # m - PBG1 internal diameter
d2ᵥ = 0.137;       # m - PBG1 external diameter
r1ᵥ = d1ᵥ/2;       # m - impeller internal radius
r2ᵥ = d2ᵥ/2;       # m - impeller external radius
# β1 = 30;        # degree - P100 internal angle
# β2 = 20;        # degree - P100 external angle
β1 = 20;          # degree - PBG1 internal angle
β2 = 15;          # degree - PBG1 external angle
# b1 = 1.75*0.5*0.0254;     # m - internal constant
# b2 = 1.75*0.5*0.0254;     # m - external constant
b1 = 0.0093;     # m - PBG1 internal constant
b2 = 0.0093;     # m - PBG1 external constant
γaᵥ = ρᵥ*(r2ᵥ^2 - r1ᵥ^2)
γbᵥ = (ρᵥ/(2*π))*(cotd(β2)/b2 - cotd(β1)/b1)

gyp = γa*ω - γb*Q

# -------------------------
# Electrical motor
@named dUₐ = Se(T)
@named dL = Mass(m = 0.1)
@named dR = Damper(c = 1.0)

@named je = Junction1(dUₐ, -dR, -dL, sgn = -1)
@named dgy = mGY(je, g = 0.01)

# -------------------------
# Pipe
@parameters RP, Ip
@variables Qₚ₁(t), Qₚ₂(t)

Rpᵥ = 128 * μᵥ * Lᵥ / (π * dᵥ^4);
Ipᵥ = ρᵥ * Lᵥ / Aᵥ;

@named pI = Mass(m=Ipᵥ)
@named pR = Damper(c=Rpᵥ)
@named pJ1 = Junction1(-pR, -pI)

# # Coupling Pipe impeller
# @named p2i = Junction0(pJ0, sgn=-1)
# @named i2p = Junction0(-pJ0)

# -------------------------
# Shaft
# T = 10
# @named sT = Se(T)
@named sI = Mass(m=1.0)
@named sR1 = Damper(c=10.0)
# @named sJ1 = Junction1(sT, -sI, -sR1, sgn=-1)
@named sJ1 = Junction1(dgy, -sI, -sR1, sgn=-1)

# -------------------------
# Impeller
@parameters RL, RT

# Impeller system
@named iI = Mass(m=1.0)
@named iRL = Damper(c = 10)
@named iRT = Damper(c = RT)

# Impeller recirculation
@named iRR = Damper(c = 10)
@named iRes = Junction1(-iRR)
@named iIn = Junction0(-pJ1, subsys = [iRes])
@named iOut = Junction0(pJ1, subsys = [-iRes])

# Impeller junction
# @named iJ1 = Junction1(-iI, -iRL, -iRT)
# @named iJ1 = Junction1(-iI, -iRL, -iRT, p2i, i2p)
@named iJ1 = Junction1(-iI, -iRL, -iRT, iIn, -iOut, couple=false)
structural_simplify(iJ1)


# mGY Connection
@named gy = mGY(sJ1, iJ1, g = gyp)

equality = [sJ1.f ~ ω, iJ1.f ~ Q]
@named sys = ODESystem(vcat(equations(gy), equality)) 
@named mdl = reducedobs(structural_simplify(sys))

# human = Dict(sJ1.sI.f => ω, iJ1.iI.f => Q, sJ1.sI.I => Is, iJ1.iI.I => If, sJ1.sT.T => T, iJ1.iRT.RT => RT, iJ1.iRL.R => RL)
# human = Dict(sJ1.sI.f => ω, iJ1.iI.f => Q, sJ1.sI.I => Is, iJ1.iI.I => If, sJ1.sT.T => T, iJ1.iRT.R => RT*Q, iJ1.iRL.R => RL)
human = Dict(sJ1.sI.f => ω, iJ1.iI.f => Q, sJ1.sI.I => Is, iJ1.iI.I => If, iJ1.iRT.R => RT*Q, iJ1.iRL.R => RL, sJ1.dgy.je.dUₐ.T => T, sJ1.sR1.R =>sJ1.sR1.R*ω)
human = merge(human, Dict(iJ1.p2i.pJ0.pI.f => Qₚ₁, iJ1.p2i.pJ0.pI.I => Ip, iJ1.p2i.pJ0.pR.R => RP))
human = merge(human, Dict(iJ1.i2p.pJ0.pI.f => Qₚ₂, iJ1.i2p.pJ0.pI.I => Ip, iJ1.i2p.pJ0.pR.R => RP))
# human = merge(human, Dict(T => τf(t)))
mdl = renamevars(mdl, human)

equations(mdl)
parameters(mdl)
# [12, 60]/8

# Ifᵥ = ρᵥ*10/A
Ifᵥ = 0.02
τᵥ = 30.0;
τᵥ = 300;

# @register τf(t)
# τf(t) = max(30, 1e-1*t)

vals = [T => τᵥ, ρ => ρᵥ, γa => γaᵥ, γb => γbᵥ, sJ1.sR1.R => 0.0005, Is => 0.02, If => 0.0005, RT => 0e8, RL => 0.0, RP => Rpᵥ, Ip => Ipᵥ/100]
# vals = [ρ => ρᵥ, γa => γaᵥ, γb => γbᵥ, sJ1.sR1.R => 100, Is => 1000.0, If => Ifᵥ, RT => 100, RL => 100, RP => Rpᵥ, Ip => Ipᵥ]
# prob = ODEProblem(mdl, [0.0, 0.0, 0.0, 0.0], (0.0, 100.0), vals)
# prob = ODEProblem(mdl, [0.0, 0.0], (0.0, 100.0), vals)
prob = ODEProblem(mdl, [0.0, 0.0, 0.0, 0.0, 0.0], (0.0, 5.0), vals)

equations(mdl)
sol = solve(prob, reltol=1e-8, abstol=1e-8)
plot(sol)

0.01*sol[sJ1.dgy.je.dL.f]


plot(sol.t, sol[1, :])
plot(sol.t, 60*sol[2, :]/(2*π))
plot(sol.t, sol[3, :]/Aᵥ)
plot!(sol.t, sol[4, :]/Aᵥ)
plot!(sol.t, sol[5, :]/Aᵥ)
# plot!(sol.t, sol[3, :]/Aᵥ)
# plot!(sol.t, sol[4, :]/Aᵥ)

i = 1
(sol[i+2, end]*60 , sol[i+1, end]/(2*π)*60)
(sol[i+2, end]/Aᵥ , sol[i+1, end]/(2*π)*60)

ωᵥ, Qᵥ = sol[2:3, end]
(-1e5*Qᵥ - 1e5*(Qᵥ^2) - γaᵥ*ωᵥ*ωᵥ + γbᵥ*Qᵥ*ωᵥ) / Ifᵥ

u₁ = ωᵥ*r1ᵥ
u₂ = ωᵥ*r2ᵥ
Vₙ₁ = Qᵥ/(2*π*b1)
Vₙ₂ = Qᵥ/(2*π*b2)
ρᵥ*Qᵥ*(u₂*Vₙ₂*cotd(β2) - u₁*Vₙ₁*cotd(β1))

ρᵥ*Qᵥ*(u₂*Vₙ₂*cotd(β2) - u₁*Vₙ₁*cotd(β1))/ωᵥ


(τᵥ - 0.1*ωᵥ)*ωᵥ

0.1*ωᵥ

τᵥ*ωᵥ


ωᵥ, Qᵥ = sol[:, end]
# Δt = 0.00001
Δt = 1
RLᵥ = 0.0
RTᵥ = 0.0
# ωᵥ, Qᵥ = 0.0, 0.0
(((γaᵥ*ωᵥ - γbᵥ*Qᵥ)*ωᵥ - RLᵥ*Qᵥ - RTᵥ*(Qᵥ^2)) / Ifᵥ)*Δt
((30 - (γaᵥ*ωᵥ - γbᵥ*Qᵥ)*Qᵥ - 0.0*ωᵥ) / 0.02)*Δt

# Qᵥ = (((γaᵥ*ωᵥ - γbᵥ*Qᵥ)*ωᵥ - RLᵥ*Qᵥ - RTᵥ*(Qᵥ^2)) / Ifᵥ)*Δt
# ωᵥ = ((30 - (γaᵥ*ωᵥ - γbᵥ*Qᵥ)*Qᵥ - 0.1*ωᵥ) / 0.02)*Δt


# --------------------------------------------------------------------------
# Complex

@parameters ρ, T, γa, γb, Is, If
@variables ω(t), Q(t)

# Fluid
# ρ = 998;     # kg/m^3 - specific mass
# Centrifugal pump P100
# de = 0.137
dᵥ = 0.108   # m - impeller external diameter
rᵥ = dᵥ/2;   # m - impeller external radiu
k1ᵥ = 2.18016;
γ = k1ᵥ/(2*rᵥ);
# γa = ρ*re^2;
# γb = ρ*γ1;

gpump = γa*ω + γb*Q

# Shaft
# T = 10
@named sT = Se(T)
@named sI = Mass(m=1.0)
@named sR1 = Damper(c=10.0)
@named sJ1 = Junction1(sT, -sI, -sR1, sgn=-1)

# Impeller
function hΔPt(ρ, ω, Q, d, k4, k5)
    # ΔP_turb = ρ*(ω^2)*(d^2)*k4*(1 - k5*Q/(ω*d^3))^2
    if ω != 0.0
        return ρ*(ω^2)*(d^2)*k4*(1 - k5*Q/(ω*d^3))^2
    else
        return 0.0
    end
end

function hΔPl(ρ, ω, Q, d, k6)
    # ΔP_loc = k6*ρ*(ω^2)*(d^2)*(Q/(ω*(d^3)))^2
    if ω != 0.0
        return k6*ρ*(ω^2)*(d^2)*(Q/(ω*(d^3)))^2
    else
        return 0.0
    end
end

@parameters d, k4, k5, k6, ΔPl, ΔPt
@register hΔPt(ρ, ω, Q, d, k4, k5)
@register hΔPl(ρ, ω, Q, d, k6)

@named iI = Mass(m=1.0)
@named iRL = GenericDamper(ΔPl)
@named iRT = GenericDamper(ΔPt)
@named iJ1 = Junction1(-iI, -iRL, -iRT, sgn=-1)

# mGY Connection
@named gy = mGY(sJ1, iJ1, g = gpump)

# iJ1.iRL.Q ~ Q, iJ1.iRT.Q ~ Q, iJ1.iRL.ω ~ ω, iJ1.iRT.ω ~ ω
equality = [sJ1.f ~ ω, iJ1.f ~ Q]
@named sys = ODESystem(vcat(equations(gy), equality)) 
@named mdl = reducedobs(structural_simplify(sys))

equations(mdl)
# iJ1.iRT.d => d, iJ1.iRT.k4 => k4, iJ1.iRT.k5 => k5, iJ1.iRL.d => d, iJ1.iRL.k6 => k6, iJ1.iRT.ρ => ρ, iJ1.iRL.ρ => ρ
human = Dict(sJ1.sI.f => ω, iJ1.iI.f => Q, sJ1.sI.I => Is, iJ1.iI.I => If, sJ1.sT.T => T, iJ1.iRT.ΔPt => hΔPt(ρ, ω, Q, d, k4, k5), 
iJ1.iRL.ΔPl => hΔPl(ρ, ω, Q, d, k6))
mdl = renamevars(mdl, human)

equations(mdl)

τᵥ = 1;
ρᵥ = 998;
k4ᵥ = 0.14011;
k5ᵥ = 8.63530;
k6ᵥ = 16.03615;
γaᵥ = ρᵥ*rᵥ^2;
γbᵥ = ρᵥ*k1ᵥ;

vals = [T => τᵥ, d => dᵥ, ρ => ρᵥ, k4 => k4ᵥ, k5 => k5ᵥ, k6 => k6ᵥ, γa => γaᵥ, γb => γbᵥ, sJ1.sR1.R => 0.1, Is => 0.02, If => 0.0001]
prob = ODEProblem(mdl, [0.0, 0.0], (0.0, 1.0), vals)

sol = solve(prob, reltol=1e-8, abstol=1e-8)

plot(sol)

A = π*(dᵥ/2)^2
plot(sol.t, sol[2, :]./A)

0.0001/A

ωᵥ, Qᵥ = (0.1, 0.001)

(τᵥ - (γaᵥ*ωᵥ - γbᵥ*Qᵥ)*Qᵥ - 100*ωᵥ) / 3000

(-hΔPl(ρᵥ, ωᵥ, Qᵥ, dᵥ, k6ᵥ) - hΔPt(ρᵥ, ωᵥ, Qᵥ, dᵥ, k4ᵥ, k5ᵥ) - (γaᵥ*ωᵥ - γbᵥ*Qᵥ)*ωᵥ)/3000


@parameters d, q, Ω

expr = d^5*ρ*Ω^3*(1/2 - 2*q/(d^3*Ω))^2*(μ/(d^2*ρ*Ω) + (1/2 - 2*q/(d^3*Ω)))
simplify(expand(expand(expr)))


0.125*ρ*d^5*ω^3 + 0.25μ*d^3*ω^2 + μ*d^3*((-2*Q)^2/d^6) + ρ*(d^5)*(-2*Q)^3/(d^9)+ 1.5*ρ*(d^5)*((-2Q)^2 * ω / (d^6)) - 2μ*Q*ω - 1.5ρ*(d^2)*(ω^2)*Q