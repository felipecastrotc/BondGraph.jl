# One-port elements. Their laws are in the variables of their bond as drawn
# (`power.e`, `power.f`); the docstrings give the conventional drawing (Karnopp §3.1).

"""
    Se(E; name)

Effort source, `e = E`. `E` is a number or any expression of `t` and system variables
(an expression of another subsystem's variables needs `@named`). Conventionally its
bond points away from the source.
"""
@oneport function Se(E; name)
    [power.e ~ E]
end

"""
    Sf(F; name)

Flow source, `f = F`. `F` is a number or any expression, as for [`Se`](@ref).
Conventionally its bond points away from the source.
"""
@oneport function Sf(F; name)
    [power.f ~ F]
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
