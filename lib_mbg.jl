using ModelingToolkit, DifferentialEquations, LinearAlgebra
using Symbolics

include("lib_bg.jl")

# =============================================================================
# Types

struct ODESystem2D{T}
    x::T
    y::T
    θ::T
end


# =============================================================================
# Ports

function Junction12D(ps2d...; name, subsys2d = [], couple = true, sgn = 1)
    # Check subsys type
    sys2d = subsys2d isa ODESystem2D ? [subsys2d] : subsys2d

    labels = ["x", "y", "θ"]
    sys, psucouple = [], []
    for i in 1:length(labels)
        namei = Symbol((name |> String) * labels[i])

        ps = [getproperty(p, Symbol(labels[i])) for p in ps2d]
        subsys = [getproperty(s, Symbol(labels[i])) for s in sys2d]

        push!(sys,
            Junction1(ps...; name = namei, subsys = subsys, couple = couple, sgn = sgn)
        )
    end

    ODESystem2D(sys...)

end

function Junction02D(ps2d...; name, subsys2d = [], couple = true, sgn = 1)
    # Check subsys type
    sys2d = subsys2d isa ODESystem2D ? [subsys2d] : subsys2d

    labels = ["x", "y", "θ"]
    sys, psucouple = [], []
    for i in 1:length(labels)
        namei = Symbol((name |> String) * labels[i])

        ps = [getproperty(p, Symbol(labels[i])) for p in ps2d]
        subsys = [getproperty(s, Symbol(labels[i])) for s in sys2d]

        push!(sys,
            Junction0(ps...; name = namei, subsys = subsys, couple = couple, sgn = sgn)
        )
    end

    ODESystem2D(sys...)

end

# =============================================================================
# Modulated Transformers

function mTF2Dtrans(a, θ, d)
    [
        1 0 a*[-sin(θ) cos(θ)]*d
        0 1 a*[-cos(θ) -sin(θ)]*d
        0 0 1
    ]
end

function mTF2Dtrans(p::Dict)
    M2Dtrans(p["a"], p["θ"], p["d"])
end

# =============================================================================
# Elements

# -----------------------------------------------------------------------------
# Planar elements
function Mass2D(; name, m::Vector = [1.0, 1.0, 1.0], u::Vector = [0.0, 0.0, 0.0])
    labels = ["x", "y", "ω"]
    sys = []
    for i in 1:length(m)

        namei = Symbol((name |> String) * labels[i])

        @named power = Power(flow = u[i])
        @unpack e, f = power

        ps = @parameters I = m[i]

        eqs = [
            D(f) ~ e / I
        ]

        push!(sys, extend(ODESystem(eqs, t, [], ps; name = namei), power))
    end

    ODESystem2D(sys...)
end

function Spring2D(; name, k::Vector = [1.0, 1.0, 1.0], x::Vector = [0.0, 0.0, 0.0])
    labels = ["x", "y", "ω"]
    sys = []
    for i in 1:length(k)
        namei = Symbol((name |> String) * labels[i])

        @named power = Power()
        @unpack e, f = power

        @variables q(t) = x[i]
        ps = @parameters C = 1 / k[i]

        eqs = [
            e ~ q / C
            D(q) ~ f
        ]

        push!(sys, extend(ODESystem(eqs, t, [q], ps; name = namei), power))
    end

    ODESystem2D(sys...)
end

function Damper2D(; name, c::Vector = [1.0, 1.0, 1.0])

    labels = ["x", "y", "ω"]
    sys = []
    for i in 1:length(c)
        namei = Symbol((name |> String) * labels[i])

        @named power = Power()
        @unpack e, f = power

        ps = @parameters R = c[i]
        eqs = [
            e ~ f * R
        ]

        push!(sys, extend(ODESystem(eqs, t, [], ps; name = namei), power))
    end

    ODESystem2D(sys...)

end
