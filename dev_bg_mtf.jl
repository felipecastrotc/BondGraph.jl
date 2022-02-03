using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Symbolics, Symbolics.Latexify

include("lib_bg.jl")

# =============================================================================
# Transformer

function mTF(subsys...; name, r = 1.0)

    # Get connections
    c = subsys isa ODESystem ? [subsys, nothing] : collect(subsys)

    # Remove nothing from c array
    pos = .!isnothing.(c)
    c = c[pos]
    @assert !isempty(c)

    # If only one subsys is passed it automatically generates an open  
    # connection
    if sum(pos) == 1
        @named power = Power()
        @unpack e, f = power
        # Set variables according to the position
        if pos[1]
            e₁, f₁ = e, f
            e₂, f₂ = c[1].e, c[1].f
        else
            e₂, f₂ = e, f
            e₁, f₁ = c[1].e, c[1].f
        end
    else
        @assert length(c) == 2
        e₁, f₁ = c[1].e, c[1].f
        e₂, f₂ = c[2].e, c[2].f
    end

    # Transformer equation
    eqs = [
        f₂ ~ f₁ * r,
        e₁ ~ e₂ * r,
    ]

    # Check if it is a modulated transformer
    if isvariable(r)
        sts, ps = [], [r]
    elseif istree(unwrap(r))
        sts = []
        ps = collect(Set(ModelingToolkit.get_variables(r)))
    else
        sts, ps = [], []
    end

    sys = ODESystem(eqs, t, sts, ps; name = name)

    if @isdefined power
        sys = extend(sys, power)
    end

    compose(sys, c...)
end


# =============================================================================
# System dynamics (palm2014) P4.65 - Solutions
# I modified the example by adding a bearing friction but the equations should
# not change significantly.
# Bond graph example: https://youtu.be/gxjbtkLANOE?t=893

# -----------------------------------------------------------------------------
# Definitions

@named rm = Mass(m = 1.0)
@named rk = Spring(k = 10.0, x = 0.1)
@named rf = Damper(c = 1.0)

@named pj = Mass(m = 1.0)
@named pf = Damper(c = 0.01)
@variables T
@named pT = Se(3)

# -----------------------------------------------------------------------------
# Using Junctions and mTR
# Bond graph example 

# Gear radius
@variables R
R = GlobalScope(R)

@named p = Junction1(pT, -pj, -pf, sgn = -1)
@named r = Junction1(-rm, -rk, -rf)
@named tf = mTF(p, r, r = 10)

equations(tf)
mdl = structural_simplify(tf)
equations(mdl)

@named sys = reducedobs(structural_simplify(tf))
eqs = equations(sys)

latexify(eqs[4])
latexify(simplify(eqs[4]))

# The last equation generated it basically states the relationship between the 
# linear and angular acceleration. I am guessing that it could be used to 
# change the coordinates. The book solutions obtain a single equation by 
# changing the coordinates. 
# It is easier to observe before replace the varaibles.

# -----------------------------------------------------------------------------
# Setting equations by hand

@variables e₁(t), e₂(t), f₁(t), f₂(t)
# @variables R  -> It generates problem with more variable than equations
R = 2.0

eqs = [
    0 ~ 3.0 - pj.e - pf.e - e₁,
    pj.f ~ pf.f,
    pf.f ~ f₁,
    e₁ ~ e₂ * R,
    f₂ ~ f₁ * R,
    0 ~ -rm.e - rf.e - rk.e + e₂,
    rm.f ~ rf.f,
    rf.f ~ rk.f,
    rk.f ~ f₂,
]

mdl = compose(ODESystem(eqs, t, [], []; name = :tst), rm, rf, rk, pj, pf)
equations(mdl)
sys = structural_simplify(mdl)
equations(sys)

latexify(equations(sys)[4])

# The manual equations after simplification resulted in the same set of equation
# obtained from the usage of the functions. So the functions are behavouring as 
# expected

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