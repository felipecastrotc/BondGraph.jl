# One-port elements. Their laws are in the variables of their bond as drawn
# (`power.e`, `power.f`); the docstrings give the conventional drawing (Karnopp §3.1).

"""
    Se(E; name)
    Se(; name)

Effort source, `e = E`. `E` is a number or any expression of `t` and system variables
(an expression of another subsystem's variables needs `@named`). Without `E`, the
effort comes from the signal connector `input::RealInput` (a modulated source, MSe).
Conventionally its bond points away from the source.
"""
@oneport function Se(E = nothing; name)
    if E === nothing
        @named input = RealInput()
        E = input.u
    end
    [power.e ~ E]
end

"""
    Sf(F; name)
    Sf(; name)

Flow source, `f = F`. `F` is a number or any expression, as for [`Se`](@ref). Without
`F`, the flow comes from `input::RealInput` (MSf). Conventionally its bond points away
from the source.
"""
@oneport function Sf(F = nothing; name)
    if F === nothing
        @named input = RealInput()
        F = input.u
    end
    [power.f ~ F]
end

"""
    De(; name)

Effort detector: outputs the effort of its bond as `output::RealOutput` and carries no
power (`f = 0`), so it doesn't change the model. It sits on a 0-junction, where every
bond shares the effort it reads (Borutzky 2011, Ch. 3 p. 116).
"""
@oneport function De(; name)
    @named output = RealOutput()
    [power.f ~ 0, output.u ~ power.e]
end

"""
    Df(; name)

Flow detector: outputs the flow of its bond as `output::RealOutput` and carries no
power (`e = 0`). It sits on a 1-junction, where every bond shares the flow it reads.
"""
@oneport function Df(; name)
    @named output = RealOutput()
    [power.e ~ 0, output.u ~ power.f]
end

"""
    Resistance(; name, R = 1.0)

Resistance, `e = R f`. Conventionally its bond points into it, so it dissipates
`e f = R f² ≥ 0`. Drawn the other way it keeps its law and delivers power.
"""
@oneport function Resistance(; name, R = 1.0)
    @parameters R = R
    [power.e ~ R * power.f]
end

"""
    Compliance(; name, C = 1.0, q = 0.0)

Compliance with displacement state `q`: `q̇ = f`, `e = q / C`. `q` is the initial
displacement. For a spring of stiffness k, `C = 1/k`. Conventionally its bond points
into it.
"""
@oneport function Compliance(; name, C = 1.0, q = 0.0)
    @parameters C = C
    @variables q(t) = q
    [D(q) ~ power.f, power.e ~ q / C]
end

"""
    Inertia(; name, I = 1.0, p = 0.0)

Inertia with momentum state `p`: `ṗ = e`, `f = p / I`. `p` is the initial momentum.
Conventionally its bond points into it.
"""
@oneport function Inertia(; name, I = 1.0, p = 0.0)
    @parameters I = I
    @variables p(t) = p
    [D(p) ~ power.e, power.f ~ p / I]
end

"""
    Mass(; name, m = 1.0, u = 0.0)

Mechanical shorthand for [`Inertia`](@ref) with mass `m` and initial velocity `u`:
`I = m`, `p = m*u`. Its bond conventionally points into it (Karnopp §3.1).
"""
Mass(; name, m = 1.0, u = 0.0) = Inertia(; name, I = m, p = m * u)

"""
    Spring(; name, k = 10, x = 0.0)

Mechanical shorthand for [`Compliance`](@ref) with stiffness `k` and initial
displacement `x`: `C = 1/k`, `q = x`. Its bond conventionally points into it
(Karnopp §3.1).
"""
Spring(; name, k = 10, x = 0.0) = Compliance(; name, C = 1 / k, q = x)

"""
    Damper(; name, c = 10)

Mechanical shorthand for [`Resistance`](@ref) with damping coefficient `c`:
`R = c`. Its bond conventionally points into it (Karnopp §3.1).
"""
Damper(; name, c = 10) = Resistance(; name, R = c)
