using BondGraph, Plots, Symbolics.Latexify
import BondGraph: t, D, mTF, tp
import ModelingToolkit: isvariable, istree
import Symbolics: unwrap, wrap

using DifferentialEquations

struct tptf end

# =============================================================================
# Transformer

function mTF(subsys...; name, r = 1.0)

    # Get connections
    c = subsys isa ODESystem ? [subsys, nothing] : collect(subsys)

    # Remove nothing from c array
    pos = .!isnothing.(c)
    c = c[pos]
    @assert !isempty(c)

    # Generate the in and out connection
    @named pin = Power(type = tptf)
    @named pout = Power(type = tptf)

    # Alias for simpler view about the mTF port
    e₁, f₁ = pin.e, pin.f
    e₂, f₂ = pout.e, pout.f

    # Transformer equation
    eqs = [f₁ * r ~ f₂, e₂ * r ~ e₁]

    # Check if it is a modulated transformer
    if isvariable(r)
        sts, ps = [], [r]
    elseif istree(unwrap(r))
        vars = collect(Set(ModelingToolkit.get_variables(r)))
        sts = filter(x -> ~isindependent(Num(x)), vars)
        ps = filter(x -> isindependent(Num(x)), vars)
    else
        sts, ps = [], []
    end

    # Apply connections
    if sum(pos) == 1
        if pos[1]
            push!(eqs, connect(c[1].power, pin))
        else
            push!(eqs, connect(pout, c[1].power))
        end
    else
        @assert length(c) == 2
        push!(eqs, connect(c[1].power, pin), connect(pout, c[2].power))
    end

    return compose(ODESystem(eqs, t, sts, ps; name = name), pin, pout)
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
@variables g(t)
R = GlobalScope(R)

@named p = Junction1(pT, pj, pf)
@named r = Junction1(rm, rk, rf)
@named tf = mTF(p, r, r = 10)
equations(tf)

mdl = compose(tf, r, p)
equations(mdl)
generate_graph(mdl)
emdl = expand_connections(mdl)
equations(emdl)
@named sys = reducedobs(structural_simplify(emdl))

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
    0 ~ 3.0 - pj.power.e - pf.power.e - e₁,
    pj.power.f ~ pf.power.f,
    pf.power.f ~ f₁,
    e₁ ~ e₂ * R,
    f₂ ~ f₁ * R,
    0 ~ -rm.power.e - rf.power.e - rk.power.e + e₂,
    rm.power.f ~ rf.power.f,
    rf.power.f ~ rk.power.f,
    rk.power.f ~ f₂,
]

mdl = compose(ODESystem(eqs, t, [], []; name = :tst), rm, rf, rk, pj, pf)
equations(mdl)
@named sys = reducedobs(structural_simplify(mdl))
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

l = 1.5
@named m = Mass(m = 1.0)
@named J = Mass(m = 1.0)
@named js = Spring(k = 0.01, x = -pi / 6)
@named jd = Damper(c = 0.01)
@named fix = Sf(0)
@named g = Se(9.81 * 1.0)

# -----------------------------------------------------------------------------
# Yellow blue Using Junctions and mTR
@variables θ(t) = 0.999 * pi
θ = GlobalScope(θ)

@named x = Junction1(-m)
@named xtf = mTF(x, r = (1.5 * cos(θ)))

@named y = Junction1(-m, g)
@named ytf = mTF(y, r = -(1.5 * sin(θ)))
@named mdl = Junction1(-J, -ytf, -xtf, couple = false)


eqs = [D(θ) ~ J.f]
mdl = extend(ODESystem(eqs, t, [θ], []; name = :mdl), mdl)

sys = structural_simplify(mdl)
equations(sys)

@named sys = reducedobs(sys)
equations(sys)

prob = ODEProblem(sys, [], (0.0, 30.0), [θ => pi / 2])
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
plot(sol.t, sol[θ])
# plot(sol)


# -----------------------------------------------------------------------------
# Book Using Junctions and mTR WRONGG!!!
@variables qx(t) = 0.0 qy(t) = 0.0
@variables θ(t) = 0.0
θ = GlobalScope(θ)

# @named xs = Junction1(fix, sgn = -1)
@named xm = Junction1(-m)
@named x0 = Junction0(-fix, -xm, sgn = -1)
@named xtf = mTF(x0, r = -(1.5 * cos(θ)))

# @named ys = Junction1(fix, sgn = -1)
@named ym = Junction1(-m, g)
@named y0 = Junction0(-fix, -ym, sgn = -1)
@named ytf = mTF(y0, r = -(1.5 * sin(θ)))

# @named j = Junction1(-js, -jd)
# @named j = Junction1(-js)
# @named j = Junction1(-m)
@named mdl = Junction1(-m, ytf, xtf, couple = false)
# @named mdl = Junction0(-j, ytf, xtf, couple = false)
# @named mdl = Junction0(ytf, xtf, couple = false)

# eqs = [D(θ) ~ ytf.f]
# eqs = [θ ~ j.js.q]
eqs = [D(θ) ~ mdl.m.f]
mdl = extend(ODESystem(eqs, t, [θ], []; name = :mdl), mdl)

sys = structural_simplify(mdl)
equations(sys)

@named sys = reducedobs(sys)
equations(sys)

prob = ODEProblem(sys, [], (0.0, 4.0), check_length = false)
sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
plot(sol)

plot(sol.t, sol[j.js.q])


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
