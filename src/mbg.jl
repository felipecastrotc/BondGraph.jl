import ModelingToolkit: states, parameters

# TODO: This is was initial conception for multidimension bond graph. It still require to develop it. It is based on a older version of the code. So, for now it is not working.

# =============================================================================
# Types
const planar_labels = ["x", "y", "θ"]

# =============================================================================
# Types

struct ODESystem2D{T}
    x::T
    y::T
    θ::T
end

# =============================================================================
# Elements

# -----------------------------------------------------------------------------
# Planar elements
function Mass2D(; name, m::Vector = [1.0, 1.0, 1.0], u::Vector = [0.0, 0.0, 0.0])
    labels = planar_labels
    sys = []
    for i = 1:length(m)

        namei = Symbol((name |> String) * labels[i])

        @named power = Power(flow = u[i])
        @unpack e, f = power

        ps = @parameters I = m[i]

        eqs = [D(f) ~ e / I]

        push!(sys, extend(ODESystem(eqs, t, [], ps; name = namei), power))
    end

    ODESystem2D(sys...)
end

function Spring2D(; name, k::Vector = [1.0, 1.0, 1.0], x::Vector = [0.0, 0.0, 0.0])
    labels = planar_labels
    sys = []
    for i = 1:length(k)
        namei = Symbol((name |> String) * labels[i])

        @named power = Power()
        @unpack e, f = power

        @variables q(t) = x[i]
        ps = @parameters C = 1 / k[i]

        eqs = [
            e ~ q / C
            D(q) ~ f
            # D(e) ~ f/C
        ]

        push!(sys, extend(ODESystem(eqs, t, [q], ps; name = namei), power))
    end

    ODESystem2D(sys...)
end

function Damper2D(; name, c::Vector = [1.0, 1.0, 1.0])

    labels = planar_labels
    sys = []
    for i = 1:length(c)
        namei = Symbol((name |> String) * labels[i])

        @named power = Power()
        @unpack e, f = power

        ps = @parameters R = c[i]
        eqs = [e ~ f * R]

        push!(sys, extend(ODESystem(eqs, t, [], ps; name = namei), power))
    end

    ODESystem2D(sys...)

end

# =============================================================================
# Ports
function Junction12D(ps2d...; name, subsys2d = [], couple = true, sgn = 1)
    # Check subsys type
    sys2d = subsys2d isa ODESystem2D ? [subsys2d] : subsys2d

    labels = planar_labels
    sys = []
    for i = 1:length(labels)
        namei = Symbol((name |> String) * labels[i])

        ps = getcomponent(ps2d, labels[i])
        subsys = getcomponent(sys2d, labels[i])

        push!(
            sys,
            Junction1(ps...; name = namei, subsys = subsys, couple = couple, sgn = sgn),
        )
    end

    ODESystem2D(sys...)

end

function Junction02D(ps2d...; name, subsys2d = [], couple = true, sgn = 1)
    # Check subsys type
    sys2d = subsys2d isa ODESystem2D ? [subsys2d] : subsys2d

    labels = planar_labels
    sys = []
    for i = 1:length(labels)
        namei = Symbol((name |> String) * labels[i])

        ps = getcomponent(ps2d, labels[i])
        subsys = getcomponent(sys2d, labels[i])

        push!(
            sys,
            Junction0(ps...; name = namei, subsys = subsys, couple = couple, sgn = sgn),
        )
    end

    ODESystem2D(sys...)

end

# =============================================================================
# Modulated Transformers

function M2Dtrans(a, θ, d)
    [
        1 0 a*[-sin(θ) cos(θ)]*d
        0 1 a*[-cos(θ) -sin(θ)]*d
        0 0 1
    ]
end

function M2Dtrans(p::Dict)
    M2Dtrans(p["a"], p["θ"], p["d"])
end

function mTF2Dtrans(subsys...; name, a = 1.0, θ = 0.0, d = [0.0, 0.0])

    labels = planar_labels

    # Get connections
    c = subsys isa ODESystem2D ? [subsys, nothing] : collect(subsys)

    # Remove nothing from c array
    pos = .!isnothing.(c)
    c = c[pos]
    @assert !isempty(c)

    # If only one subsys is passed it automatically generates an open  
    # connection
    if sum(pos) == 1
        sys_vec = ODESystem[]
        for l in labels
            namei = Symbol((name |> String) * l)

            @named power = Power()
            @unpack e, f = power

            push!(sys_vec, extend(ODESystem(Equation[], t, [], []; name = namei), power))
        end
        sys = ODESystem2D(sys_vec...)

        # Set variables according to the position
        if pos[1]
            E₁, F₁ = getvecvar(sys, "e"), getvecvar(sys, "f")
            E₂, F₂ = getvecvar(c[1], "e"), getvecvar(c[1], "f")
        else
            E₂, F₂ = getvecvar(sys, "e"), getvecvar(sys, "f")
            E₁, F₁ = getvecvar(c[1], "e"), getvecvar(c[1], "f")
        end
    else
        @assert length(c) == 2
        E₁, F₁ = getvecvar(c[1], "e"), getvecvar(c[1], "f")
        E₂, F₂ = getvecvar(c[2], "e"), getvecvar(c[2], "f")
    end

    # Translation matrix
    M = M2Dtrans(a, θ, d)

    # Transformer equation
    eqs = []
    push!(eqs, F₂ .~ M * F₁)
    push!(eqs, E₁ .~ M' * E₂)
    eqs = collect(zip(eqs...))

    sys_vec = []
    for i = 1:length(labels)
        l = labels[i]
        namei = Symbol((name |> String) * l)
        sym = Symbol(l)

        num_vars = [x for x in (a, θ, d[1], d[2]) if x isa Num]
        sts = [x for x in num_vars if !isindependent(x)]
        ps = [x for x in num_vars if isindependent(x)]

        if @isdefined sys
            updt_sys = getproperty(sys, sym)

            eq = collect(eqs[i])
            updt_sys = extend(ODESystem(eq, t, sts, ps; name = namei), updt_sys)
            updt_sys = compose(updt_sys, [getproperty(j, sym) for j in c]...)

            push!(sys_vec, updt_sys)
        else
            eq = collect(eqs[i])
            updt_sys = ODESystem(eq, t, sts, ps; name = namei)
            updt_sys = compose(updt_sys, [getproperty(j, sym) for j in c]...)
            push!(sys_vec, updt_sys)
        end
    end

    ODESystem2D(sys_vec...)
end

# =============================================================================
# Utils
struct SgnODESystem2D
    x::ODESystem
    y::ODESystem
    θ::ODESystem
    sgn::Num
end

-(sys::ODESystem2D) = SgnODESystem2D(sys.x, sys.y, sys.θ, -1)
+(sys::ODESystem2D) = SgnODESystem2D(sys.x, sys.y, sys.θ, 1)

function getvecvar(sys::ODESystem2D, var::String)
    sym = Symbol(var)
    [getproperty(getproperty(sys, Symbol(l)), sym) for l in planar_labels]
end

function getcomponent(sys, comp::String)

    sgn = Symbol("sgn")
    component = Symbol(comp)

    out = SgnODESystem[]
    for sub in sys
        c = getproperty(sub, component)

        if !hasproperty(sub, sgn)
            push!(out, +c)
        else
            if getproperty(sub, sgn) > 0
                push!(out, +c)
            else
                push!(out, -c)
            end
        end

    end

    out

end

function setveceq(sys::ODESystem2D, var::String)
    sym = Symbol(var)
    [getproperty(getproperty(sys, Symbol(l)), sym) for l in planar_labels]
end

function get2Dsystem(sys::ODESystem2D; name)

    sts = states(sys)
    ps = parameters(sys)
    eqs = equations(sys)

    ODESystem(eqs, t, sts, ps, name = name)
end

# -----------------------------------------------------------------------------
# Extending
function equations(sys::ODESystem2D)
    eqs = Equation[]
    for l in planar_labels
        push!(eqs, equations(getproperty(sys, Symbol(l)))...)
    end
    eqs
end

function states(sys::ODESystem2D)
    sts = []
    for l in planar_labels
        push!(sts, states(getproperty(sys, Symbol(l)))...)
    end
    sts
end

function parameters(sys::ODESystem2D)
    ps = []
    for l in planar_labels
        push!(ps, parameters(getproperty(sys, Symbol(l)))...)
    end
    ps
end
