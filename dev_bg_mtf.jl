using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Symbolics, Symbolics.Latexify

include("lib_bg.jl")

# =============================================================================
# Transformer

function mTF(subsys...; name, r = 1.0, couple = true)
    # Get connections
    c = collect(Base.Flatten([subsys]))

    # Get variables
    # vars = extract_vars(r)
    # if isnothing()

    if couple

        pos = .!isnothing.(c)
        @assert length(sum(pos)) == 1

        @named power = Power()
        @unpack e, f = power

        # Effort and flow 
        print(pos)
        if pos[1]
            eqs = [
                c[1].f * r ~ f,
                e * r ~ c[1].e,
            ]
        else
            eqs = [
                f * r ~ c[2].f,
                c[2].e * r ~ e,
            ]
        end
        # Build subsystem
        # ODESys
        filter!(!isnothing, c)
        compose(extend(ODESystem(eqs, t, [], []; name = name), power), c...)
    else

        @assert length(c) == 2

        # Effort and flow 
        eqs = [
            c[1].f * r ~ c[2].f,
            c[2].e * r ~ c[1].e,
        ]
        # Build subsystem
        compose(ODESystem(eqs, t, [], []; name = name), c...)
    end
end


# =============================================================================
# System dynamics (palm2014) P4.65 - Solutions
# I modified the example by adding a bearing friction but the equations should
# not change significantly

# -----------------------------------------------------------------------------
# Definitions

m = 1.0
X = 1.0
k = 10.0
c = 1.0
@named bm = Mass(m = m)
@named bs = Spring(k = k, x = X)
@named bd = Damper(c = c * 2)

@named rm = Mass(m = m * 2)
@named rd = Damper(c = c * 5)

# -----------------------------------------------------------------------------
# Using Junctions and mTR

@named r = Junction1(rm, rd)
@named b = Junction1(bm, bs, bd)
@named mdl = mTF(r, b, r = 2, couple = false)

sys = structural_simplify(mdl)
equations(sys)
parameters(sys)
states(sys)

prob = ODEProblem(sys, [], (0.0, 10.0))
sol = solve(prob)
plot(sol)
plot(sol.t, sol[r.e])
plot(sol.t, sol[b.bs.q])

@named sysr = reducedobs(sys)
equations(sysr)

# -----------------------------------------------------------------------------
# Setting equations by hand

@variables e₁, e₂, f₁, f₂
# @variables R  -> It generates problem with more variable than equations
R = 2.0

eqs = [
    e₁ ~ rd.e + rm.e,
    rd.f ~ rm.f,
    f₁ ~ rd.f,
    e₁ ~ e₂ * R,
    f₁ * R ~ f₂,
    e₂ ~ bm.e + bd.e + bs.e,
    f₂ ~ bm.f,
    f₂ ~ bd.f,
    f₂ ~ bs.f,
]

mdl = compose(ODESystem(eqs, t, [], []; name = :tst), rm, rd, bm, bd, bs)
sys = structural_simplify(mdl)
dae_index_lowering(ODESystem(eqs, t, [], []; name = :tst))
structural_simplify(ODESystem(eqs, t, [], []; name = :tst))


# -----------------------------------------------------------------------------
# Setting equations by hand simplified

ps = @parameters R
ps = [2.0]

eqs = [
    rd.e + rm.e ~ (bm.e + bd.e + bs.e) * ps[1],
    rd.f ~ rm.f,
    rd.f * ps[1] ~ bm.f,
    bd.f ~ bm.f,
    bs.f ~ bd.f,
]

mdl = compose(ODESystem(eqs, t, [], ps; name = :tst), bm, bd, bs, rm, rd)
mdl = compose(ODESystem(eqs, t, [], []; name = :tst), bm, bd, bs, rm, rd)
sys = structural_simplify(mdl)
equations(sys)

prob = ODEProblem(sys, [], (0.0, 20.0))
sol = solve(prob)
plot(sol)

states(sys)

@named sysr = reducedobs(sys)
equations(sysr)

# The equations obtained by hand and by the mTR and junctions functions were 
# the same for some reason the DAE obtained is not in the canonical form. They
# do not depend on the states variables but on effort and flow :/


# =============================================================================
# Pendulum

# Specifying a time variable forcing function -> https://mtk.sciml.ai/stable/tutorials/ode_modeling/#Specifying-a-time-variable-forcing-function

# -----------------------------------------------------------------------------
# Definitions

m = 1.0
X = pi / 2
k = 0.01
c = 2.0
l = 1.5
@named lm = Mass(m = m)

@named js = Spring(k = k, x = X)
@named jd = Damper(c = c)

# -----------------------------------------------------------------------------
# Using Junctions and mTR
@variables θ
θ = GlobalScope(θ)

@named j = Junction1(js, jd)
@named lx = Junction0(lm)
@named ly1 = Junction1(lm, se = [-9.81])
@named ly = Junction0(ly1)
@named mtfx = mTF(lx, r = l * sin(θ))
@named mtfy = mTF(ly, r = l * cos(θ))
@named mdl = Junction0(mtfx, mtfy, j, couple = false)


# -----------------------------------------------------------------------------
# Setting equations by hand

@variables e₁, e₂, f₁, f₂
# @variables R  -> It generates problem with more variable than equations
R = 2.0

eqs = [
    e₁ ~ rd.e + rm.e,
    rd.f ~ rm.f,
    f₁ ~ rd.f,
    e₁ ~ e₂ * R,
    f₁ * R ~ f₂,
    e₂ ~ bm.e + bd.e + bs.e,
    f₂ ~ bm.f,
    f₂ ~ bd.f,
    f₂ ~ bs.f,
]

# -----------------------------------------------------------------------------

equations(mtfx)

eqs = [θ ~ j.js.q]
mdl = extend(ODESystem(eqs, t, [], []; name = :mdl), mdl)

sys = structural_simplify(mdl)

equations(sys)
parameters(sys)
states(sys)

prob = ODEProblem(sys, [], (0.0, 4.0))
sol = solve(prob)
plot(sol)
plot(sol.t, sol[r.e])
plot(sol.t, sol[b.bs.q])

# -----------------------------------------------------------------------------
# Setting equations by hand

@variables e[1:11](t), f[1:11](t)
@variables θ = 0.0
@variables e₇(t) = 0.0
@variables f₈(t) = 0.0

@parameters l = 1.5
θ = GlobalScope(θ)
R = l * sin(θ)

eqs = [
    # Sf_x = 0
    # 0 ~ e[1] + e[2],
    f[1] ~ f[2],
    # Mass x
    0 ~ e[3] - lm.e,
    f[3] ~ lm.f,
    # Junction0x
    0 ~ f[2] - f[3] - f[4],
    e[3] ~ e[2],
    e[2] ~ e[4],
    # MTFx
    f[4] * R ~ f[5],
    e[5] * R ~ e[4],
    # Junction1 θ
    0 ~ e[6] - jd.e - js.e,
    f[6] ~ jd.f,
    jd.f ~ js.f,
    # Junction0
    0 ~ f[5] - f[6] + f[7],
    e[5] ~ e[6],
    e[6] ~ e₇,
    # MTFy
    f₈ * R ~ f[7],
    e₇ * R ~ e[8],
    # Junction0y
    0 ~ f[9] - f₈ - f[11],
    e[8] ~ e[9],
    e[9] ~ e[11],
    # Mass y
    0 ~ e[11] - lm.e + lm.I * 9.81,
    f[3] ~ lm.f,
    # Sf_y = 0
    # 0 ~ e[10] + e[9],
    f[10] ~ f[9],
    # General equations
    θ ~ js.q,
    f[1] ~ 0,
    f[10] ~ 0,
]

mdl = compose(ODESystem(eqs, t, [e₇, f₈], [l]; name = :mdl), jd, js, lm)

parameters(mdl)

sys = structural_simplify(mdl)
sys

equations(sys)
states(sys)
parameters(sys)

states(sys)

prob = ODEProblem(sys, [], (0.0, 4.0))
sol = solve(prob)
plot(sol)



equations(sys)
structural_simplify(sys)