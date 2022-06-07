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
    @variables p(t)

    eqs = [
        # It(e) ~ f*I
        # D(p) ~ e,
        # p ~ I * f,
        # f ~ p/I,
        D(f) ~ e / I,
    ]
    ODESystem(eqs, t, [e, f, p], ps; name = name)
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


mutable struct BgODESystem
    ode::ODESystem
    sgn::Num
    subsys::Array{BgODESystem, 1}
    con::Symbol
    name::Symbol
end


function BgODESystem(ode::ODESystem, con::Symbol)
    BgODESystem(ode, 1, ODESystem[], con, ode.name)
end

function Base.getproperty(object::BgODESystem, item::Symbol)
    bg_fields = [:sgn, :subsys, :con, :ode, :name]
    if item in bg_fields
        getfield(object, item)
    elseif item == :e
        "wow"
    elseif item == :f
        "wowf"
    else
        getproperty(getproperty(object, :ode), item)
    end
end

# bg = BgODESystem(m, :j1)
sumvar(con, var) = sum(c -> getproperty(c, var), con)

function equalityeqs(C, sym)
    eqs = Vector{Equation}(undef, length(C) - 1)
    for i = 1:(length(C)-1)
        # f₁ = C[i].sgn * getproperty(C[i].ode, sym)
        # f₂ = C[i+1].sgn * getproperty(C[i+1].ode, sym)
        f₁ = getproperty(C[i], sym)
        f₂ = getproperty(C[i+1], sym)
        eqs[i] = f₁ ~ f₂
    end
    eqs
end


test(1)

# function Junction1(ps...; name, subsys = [], coupling = [])
    ps = [m, d]
    # Get connections
    ps = collect(Base.Flatten([ps]))
    # con = addsgnODE.(collect(Base.Flatten([ps, subsysv])))
    con = vcat(ps, subsys)

    coupling = []
    e, f = Power(n = length(coupling))
    
    sum(con, con)
    # Σ efforts

    ?
    eqs = [0 ~ sum(con, :e) + sum(e, coupling)]
    # f₁ = f₂ = f₃
    eqs = vcat(eqs, equalityeqs(con, :f, couple = couple))
    # Remove empty equations
    if couple
        filter!(x -> filterexpr(x, ignore = [e, f]), eqs)
    end

    # Build subsystem
    if couple
        sys = ODESystem(eqs, t, [e, f], []; name = name)
    else
        sys = ODESystem(eqs, t, [], []; name = name)
    end

    out = compose(sys, rmsgnODE.(ps)...)
    # out = addsubsys(out, ps)
    if length(subsysv) == 0
        return BgODESystem(out, :j1)
    else
        BgODESystem(out, 1, rmsgnODE.(subsysv), :j1)
    end
# end





@named m = Mass()
@named d = Damper()

con = [m, d];

eqs = [0 ~ sumvar(con, :e) ]
eqs = vcat(eqs, equalityeqs(con, :f))

j1 = ODESystem(eqs, t, [], []; name = :j1)

subsys = []


@named m2 = Mass()
@named d2 = Damper()

con = [m2, d2];
subsys = [j1]
foreign =

j1
eqs = [0 ~ sumvar(con, :e) ]
eqs = vcat(eqs, equalityeqs(con, :f))

j2 = ODESystem(eqs, t, [], []; name = :j2)

subsys = []








function Power(; effort = 0.0, flow = 0.0)
    sts = @variables e(t) = effort f(t) = flow
    sts
end

@connector function J1(; name, effort = 0.0, flow = 0.0)
    sts = @variables e(t) = effort [connect = Flow] f(t) = flow
    ODESystem(Equation[], t, sts, []; name=name)
end

@connector function J0(; name, effort = 0.0, flow = 0.0)
    sts = @variables e(t) = effort f(t) = flow [connect = Flow]
    ODESystem(Equation[], t, sts, []; name=name)
end


@named pj1 = J1()
@named pj0 = J0()

@named m = Mass(m = 10)
@named d = Mass(m = 10)

def = ModelingToolkit.getdefault(m.f)
@variables e(t) = def [connect = Flow]


ODESystem(connect(m, d), t; name = :a)



connect(pj1.e, m.e)



function Mass(; name, m = 1.0, u = 0.0)
    e, f = Power(flow = u)
    
    ps = @parameters I = m
    @variables p(t)

    eqs = [
        # It(e) ~ f*I
        # D(p) ~ e,
        # p ~ I * f,
        # f ~ p/I,
        D(f) ~ e / I,
    ]
    ODESystem(eqs, t, [e, f, p], ps; name = name)
    # extend(ODESystem(eqs, t, [], ps; name = name), power)
    # ODESystem(eqs, t, [], ps; name = name)
    # extend(ODESystem(eqs, t, [p], ps; name = name), power)
end


function Spring(; name, k = 10, x = 0.0)
    @named power = Power(flow = u)
    @unpack e, f = power

    @variables q(t) = x
    ps = @parameters C = 1 / k

    eqs = [
        e ~ q / C
        D(q) ~ f
        # D(e) ~ f/C
    ]
    extend(ODESystem(eqs, t, [q], ps; name = name), power)
end

function Damper(; name, c = 10)
    @named power = Power(flow = u)
    @unpack e, f = power

    ps = @parameters R = c
    eqs = [e ~ f * R]
    extend(ODESystem(eqs, t, [], ps; name = name), power)
end


@connector function junction1(;name)
    sts = @variables e(t)=1.0 [connect = Flow] f(t)=1.0
    ODESystem(Equation[], t, sts, []; name=name)
end

@named j1 = junction1()

function Junction1(ps...; name)
     @named j1 = junction1()

     eqs = Equation[]
     for p in ps
        push!(eqs, connect(p.e, j1.e))
        push!(eqs, connect(p.f, j1.f))
     end

     compose(ODESystem(eqs, t, [], []; name = name), j1)
end



@named m = Mass(m = 10)
@named d = Damper(c = 10)

@named test = Junction1(m, d)


@named j1 = junction1()
eqs = Equation[]

ps = [m, d]

p = ps[1]
push!(eqs, connect(p.e, j1.e))
push!(eqs, connect(p.f, j1.f))

for p in ps
    push!(eqs, connect(p.e, j1.e))
    push!(eqs, connect(p.f, j1.f))
end






@connector function TwoPhaseFluidPort(; name, P = 0.0, m_flow = 0.0, h_outflow = 0.0)
    vars = @variables h_outflow(t)=h_outflow [connect = Stream] m_flow(t)=m_flow [
        connect = Flow,
    ] P(t)=P
    ODESystem(Equation[], t, vars, []; name = name)
end