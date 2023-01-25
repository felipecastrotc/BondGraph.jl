struct ConvConst
    g::Float64      # Gravity
    H::Float64      # Convert Pascal to head
    function ConvConst()
        g = 9.81        # Gravity
        H = 9.81*1000   # MMCA
        new(g, H)
    end
end

struct Fluid
    a::Float64      # (m/s) Speed of sound
    ρ::Float64      # (kg/m^3) Specific mass
    μ::Float64      # (Pa*s) Viscosity
end

struct Impeller
    # Impeller geometry
    γa::Float64
    γb::Float64
    d::Float64         # (m) Impeller diameter
    l::Float64        # (m) Equivalent length
    # These are an estimation thinking the impeller as a channel
    d_m::Float64       # (m) Diâmetro médio do impelidor
    h_m::Float64       # (m) Altura média do impelidor
    l_m::Float64         # (m) Comprimento médio do impelidor
    I::Float64         # (m) Impeller Inertia
    function Impeller(γa, γb, d, l, d_m, h_m, l_m, fluid::Fluid)
        # d (mm)
        # l (m)
        # d_m (mm)
        # h_m (mm)
        # l_m (mm)
        # da (mm)
        # la (m)

        # Convert mm to m
        di = d / 1000
        di_m = d_m / 1000
        h_m = h_m / 1000
        l_m = l_m / 1000

        # Impeller inertia
        I = fluid.ρ * l_m / (π * d_m * h_m)
        new(γa, γb, d, l, d_m, h_m, l_m, I)
    end
end

struct Shaft
    # Axis
    c::Float64         # Axis damping coefficient
    d::Float64         # (m) Axis diameter
    l::Float64         # (m) Axis length
    ρ::Float64         # (kg/m^3) Specific mass for steel
    I::Float64         # Axis Inertia
    function Shaft(c, d, l, ρ)
        # l (m)
        # d (mm)

        # Convert mm to m
        d = d / 1000
        # Axis inertia
        I = ((ρ * π * (d / 2)^2 * l) * (d / 2)^2) / 2
        new(c, d, l, ρ, I)
    end
end

struct Pipe
    d::Float64        # (m) Pipe diameter
    l::Float64       # (m) Pipe length
    ϵ::Float64        # (m) Steel rugosity
    A::Float64        # (m^2) Pipe area
    # Pipeline dynamics -> using Tanaka pipe
    # https://sci-hub.se/10.1109/IECON.2000.972506
    I::Float64
    C::Float64
    function Pipe(d, l, ϵ, fluid::Fluid)
        # d (mm)
        # l (m)

        d = d / 1000
        A = π * (d / 2)^2
        I = fluid.ρ * l / A
        C = A * l / (fluid.ρ * fluid.a^2)

        new(d, l, ϵ, A, I, C)
    end
end

struct System
    fluid::Fluid
    pipe::Pipe
    shaft::Shaft
    pump::Impeller
    conv::ConvConst
    function System(fluid, pipe, shaft, pump)
        new(fluid, pipe, shaft, pump, ConvConst())
    end
end

struct Inputs
    p::Vector
    func::Function
end

struct Sim
    sys::System
    u0::Vector
    tspan::Tuple
    model::Function
    input::Inputs
end

function type2dict(x)
    return Dict(key => getfield(x, key) for key ∈ fieldnames(typeof(x)))
end

function sim2dict(x::System)
    return Dict(
        :fluid => type2dict(x.fluid),
        :pipe => type2dict(x.pipe),
        :shaft => type2dict(x.shaft),
        :pump => type2dict(x.pump),
        :g => g,
        :H => 1000*g,
    )
end

function sim2dict(x::Sim)
    return Dict(
        :sys => sim2dict(x.sys),
        :u0 => x.u0,
        :tspan => x.tspan,
        :model => string(x.model),
        :input => Dict(:func => string(x.input.func), :p => x.input.p),
    )
end

function get_unpack(x::Sim)
    return eval(Symbol(split(string(sim.model), "!")[1] * "_unpack"))
end

function get_sim_ps(x::Sim)
    return eval(Symbol(split(string(sim.model), "!")[1] * "_ps"))
end
