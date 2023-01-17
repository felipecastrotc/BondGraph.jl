g = 9.81        # Gravity
H = 9.81*1000   # MMCA

struct Fluid
    a::Float64      # (m/s) Speed of sound
    ρ::Float64      # (kg/m^3) Specific mass
    μ::Float64      # (Pa*s) Viscosity
end

struct Impeller
    # Impeller geometry
    γa::Float64
    γb::Float64
    di::Float64         # (m) Impeller diameter
    Δsi::Float64        # (m) Equivalent length
    # These are an estimation thinking the impeller as a channel
    di_m::Float64       # (m) Diâmetro médio do impelidor
    hi_m::Float64       # (m) Altura média do impelidor
    li::Float64         # (m) Comprimento médio do impelidor
    Ii::Float64         # (m) Impeller Inertia
    function Impeller(γa, γb, di, Δsi, di_m, hi_m, li, fluid::Fluid)
        # di (mm)
        # Δsi (m)
        # di_m (mm)
        # hi_m (mm)
        # li (mm)
        # da (mm)
        # la (m)

        # Convert mm to m
        di = di / 1000
        di_m = di_m / 1000
        hi_m = hi_m / 1000
        li = li / 1000

        # Impeller inertia
        Ii = fluid.ρ * li / (π * di_m * hi_m)
        new(γa, γb, di, Δsi, di_m, hi_m, li, Ii)
    end
end

struct Shaft
    # Axis
    ca::Float64         # Axis damping coefficient
    da::Float64         # (m) Axis diameter
    la::Float64         # (m) Axis length
    ρa::Float64         # (kg/m^3) Specific mass for steel
    Ia::Float64         # Axis Inertia
    function Shaft(ca, da, la, ρa)
        # la (m)
        # da (mm)

        # Convert mm to m
        da = da / 1000
        # Axis inertia
        Ia = ((ρa * π * (da / 2)^2 * la) * (da / 2)^2) / 2
        new(ca, da, la, ρa, Ia)
    end
end

struct Pipe
    d::Float64        # (m) Pipe diameter
    Δs::Float64       # (m) Pipe length
    ϵ::Float64        # (m) Steel rugosity
    A::Float64        # (m^2) Pipe area
    # Pipeline dynamics -> using Tanaka pipe
    # https://sci-hub.se/10.1109/IECON.2000.972506
    L::Float64
    C::Float64
    function Pipe(d, Δs, ϵ, fluid::Fluid)
        # d (mm)
        # Δs (m)

        d = d / 1000
        A = π * (d / 2)^2
        L = fluid.ρ * Δs / A
        C = A * Δs / (fluid.ρ * fluid.a^2)

        new(d, Δs, ϵ, A, L, C)
    end
end

struct System
    fluid::Fluid
    pipe::Pipe
    shaft::Shaft
    pump::Impeller
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
