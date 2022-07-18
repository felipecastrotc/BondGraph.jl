using ModelingToolkit

@variables t
D = Differential(t)

function Power(; effort = 0.0, flow = 0.0, n = 1)
    if n > 1
        sts = @variables e[1:n](t) = effort f[1:n](t) = flow
    elseif n <= 0
        sts = missing, missing
    else
        sts = @variables e(t) = effort f(t) = flow
    end
    sts
end

function Mass(; name, m = 1.0, u = 0.0)
    e, f = Power(flow = u)
    
    ps = @parameters I = m
    # @variables p(t)

    eqs = [
        # It(e) ~ f*I
        # D(p) ~ e,
        # p ~ I * f,
        # f ~ p/I,
        D(f) ~ e / I,
    ]
    ODESystem(eqs, t, [e, f], ps; name = name)
    # ODESystem(eqs, t, [e, f, p], ps; name = name)
    # extend(ODESystem(eqs, t, [], ps; name = name), power)
    # ODESystem(eqs, t, [], ps; name = name)
    # extend(ODESystem(eqs, t, [p], ps; name = name), power)
end

function Damper(; name, c = 10, u=1.0)
    e, f = Power(flow = u)
    
    ps = @parameters R = c
    eqs = [e ~ f * R]
    ODESystem(eqs, t, [e, f], ps; name = name)
end

struct BgODESystem
    ode::ODESystem
    sgn::Num
    subsys::Array{ODESystem, 1}
    couple::Num
end

function BgODESystem(ode::ODESystem; sgn=1, subsys=ODESystem[], couple=0)
    BgODESystem(ode, sgn, subsys, couple)
end

function BgODESystem(bg::BgODESystem, couple)
    BgODESystem(bg.ode, bg.sgn, bg.subsys, couple)
end

function Base.getindex(object::BgODESystem, i)
    BgODESystem(object, i)
end

function Base.getproperty(object::BgODESystem, item::Symbol)
    bg_fields = [:sgn, :subsys, :couple, :ode]
    if item in bg_fields
        getfield(object, item)
    elseif item in [:e, :f]
        if object.couple == 0
            getproperty(getproperty(object, :ode), item)
        else
            getproperty(getproperty(object, :ode), item)[object.couple]
        end
    else
        getproperty(getproperty(object, :ode), item)
    end
end

# bg = BgODESystem(m, :j1)
sumvar(con, var) = sum(c -> getproperty(c, var), con)

function equalityeqs(con, sym, signs, vars)
    C = con
    eqs = Equation[]
    for i = 1:(length(C)-1)
        f₁ = getproperty(C[i], sym)
        f₂ = getproperty(C[i+1], sym)
        push!(eqs, f₁ ~ f₂)
    end

    if length(signs) > 1
        f₁ = getproperty(C[end], sym)
        for (sgn, v) in zip(signs, vars)
            push!(eqs, f₁ ~ sgn*v)
        end

    elseif length(signs) == 1
        f₁ = getproperty(C[end], sym)
        f₂ = signs[1] * vars
        push!(eqs, f₁ ~ f₂)
    end

    eqs
end


function Junction1(ps...; name, subsys = [], coupling = [])
    # Get connections
    ps = collect(Base.Flatten([ps]))
    con = addsgn.(collect(Base.Flatten([ps, subsys])))

    e, f = Power(n = length(coupling))

    # Σ efforts
    eqs = [0 ~ sumvar(con, :e) + sumcoupling(e, coupling)]
    # f₁ = f₂ = f₃
    eqs = vcat(eqs, equalityeqs(con, :f, coupling, f))

    # Build subsystem
    sts = vcat([[e[i], f[i]] for (i, v) in enumerate(coupling)]...)
    sys = ODESystem(eqs, t, sts, [], name = name)
    sys = compose(sys, rmsgn.(ps)...)

    return BgODESystem(sys)
end

@named m = Mass()
@named d = Damper()

ps = [m, d]
subsys = [b1[1]]
coupling = []
coupling = [1,1]

@named b1 = Junction1(m, d, coupling = [1, 1])
@named b2 = Junction1(m, d, subsys=[b1[1]])
@named b3 = Junction1(m, d, subsys=[b1[2]])


equations(compose(b1.ode, b2.ode, b3.ode))