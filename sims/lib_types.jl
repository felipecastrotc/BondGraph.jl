using YAML

struct ConvConst
    g::Float64      # Gravity
    H::Float64      # Convert Pascal to head
    function ConvConst()
        g = 9.81        # Gravity
        H = 9.81 * 1000   # MMCA
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
    # Shock loss
    ks::Float64         # Shock loss constant
    q_star::Float64     # (m^3/s) Design flow rate
    # Parameters from fit
    K::Vector{Float64}        # Fit parameters
    c_f::Float64

    function Impeller(γa, γb, d, l, d_m, h_m, l_m, ks, q_star, fluid::Fluid)
        # d (mm)
        # l (m)
        # d_m (mm)
        # h_m (mm)
        # l_m (mm)
        # d (mm)
        # l (m)
        # ks (-)
        # q_star (mˆ3/s)

        # Convert mm to m
        d = d / 1000
        d_m = d_m / 1000
        h_m = h_m / 1000
        l_m = l_m / 1000

        # Impeller inertia
        I = fluid.ρ * l_m / (π * d_m * h_m)
        new(γa, γb, d, l, d_m, h_m, l_m, I, ks, q_star, [0.0], 0.0)

    end

    function Impeller(γa, γb, d, l, d_m, h_m, l_m, fluid::Fluid, c_f::Number=0.0)
        Impeller(γa, γb, d, l, d_m, h_m, l_m, 0.0, 1.0, fluid, c_f)
    end

    function Impeller(p::Impeller, fluid::Fluid)
        # Update the Impeller inertia
        l_m, d_m, h_m = p.l_m, p.d_m, p.h_m
        I = fluid.ρ * l_m / (π * d_m * h_m)
        new(p.γa, p.γb, p.d, p.l, p.d_m, p.h_m, p.l_m, I, p.ks, p.q_star, p.K, p.c_f)
    end

    function Impeller(path::String, key::String, d_m, h_m, l_m, fluid::Fluid, c_f::Number=0.0)

        # Load constants
        K = (YAML.load_file(path))[key]

        # Convert mm to m
        d_m = d_m / 1000
        h_m = h_m / 1000
        l_m = l_m / 1000
        I = fluid.ρ * l_m / (π * d_m * h_m)

        new(1.0, 1.0, 1.0, 1.0, d_m, h_m, l_m, I, 1.0, 1.0, K, c_f)
    end
    function Impeller(path::String, key::Vector{String}, d_m, h_m, l_m, fluid::Fluid, c_f::Number=0.0)

        # Load constants
        K = (YAML.load_file(path))[key[1]][key[2]]["K"]

        # Convert mm to m
        d_m = d_m / 1000
        h_m = h_m / 1000
        l_m = l_m / 1000
        I = fluid.ρ * l_m / (π * d_m * h_m)

        new(1.0, 1.0, 1.0, 1.0, d_m, h_m, l_m, I, 1.0, 1.0, K, c_f)
    end
end

struct Shaft
    # Axis
    c::Float64         # Axis damping coefficient
    d::Float64         # (m) Axis diameter
    l::Float64         # (m) Axis length
    ρ::Float64         # (kg/m^3) Specific mass for steel
    I::Float64         # Axis Inertia
    # Parameters from fit
    K::Vector{Float64}        # Fit parameters
    function Shaft(c, d, l, ρ)
        # l (m)
        # d (mm)

        # Convert mm to m
        d = d / 1000
        # Axis inertia
        I = ((ρ * π * (d / 2)^2 * l) * (d / 2)^2) / 2
        new(c, d, l, ρ, I, [0.0])
    end

    function Shaft(d, l, ρ, path::String, key::Vector{String})
        # l (m)
        # d (mm)

        # Load constants
        K = (YAML.load_file(path))[key[1]][key[2]]["K"]

        # Convert mm to m
        d = d / 1000
        # Axis inertia
        I = ((ρ * π * (d / 2)^2 * l) * (d / 2)^2) / 2
        new(0.0, d, l, ρ, I, K)
    end
end

struct Valve
    C_v::Function
    C_v_ps::Vector
    F_l::Float64
    F_p::Float64
    F_d::Function
    F_d_ps::Vector
    N_1::Float64
    N_2::Float64
    N_4::Float64
    function Valve(C_v, C_v_ps, N_1)
        new(C_v, C_v_ps, 1.0, 1.0, (x) -> x, [1.0], N_1, 1.0, 1.0)
    end
    function Valve(path::String, F_p::Number, key::String="datasheet")
        # path: directory to a YAML file containing the valve parameters
        ps = (YAML.load_file(path))[key]
        new(eval(Symbol(ps["C_v_func"])), ps["C_v"], ps["F_l"], F_p, eval(Symbol(ps["F_d_func"])), ps["F_d"], ps["N_1"], ps["N_2"], ps["N_4"])
    end
end

struct Twin
    ω::Float64
    k_1::Float64
    k_2::Float64
    h::Float64

    function Twin(ω::Float64, path::String, key::String)
        # Load constants
        ps = (YAML.load_file(path))[key]
        new(ω, ps["kd"], ps["ks"], ps["h"])
    end
    function Twin(ω::Float64, path::String, key::Vector{String})
        # Load constants
        ps = (YAML.load_file(path))[key[1]][key[2]]
        new(ω, ps["kd"], ps["ks"], ps["h"])
    end
end

mutable struct Pipe{T}
    d::Float64        # (m) Pipe diameter
    l::Float64       # (m) Pipe length
    ϵ::Float64        # (m) Steel rugosity
    A::Float64        # (m^2) Pipe area
    # Pipeline dynamics -> using Tanaka pipe
    # https://sci-hub.se/10.1109/IECON.2000.972506
    I::Float64
    C::Float64
    # Localized loss
    k::Float64      # (?) K loss
    others::T
    function Pipe(d, l, ϵ, fluid::Fluid, k=0.0)
        # d (mm)
        # l (m)

        d = d / 1000
        A = π * (d / 2)^2
        I = fluid.ρ * l / A
        C = A * l / (fluid.ρ * fluid.a^2)

        new{Nothing}(d, l, ϵ, A, I, C, k, nothing)
    end

    function Pipe(d, l, ϵ, fluid::Fluid, other::Union{Valve,Vector{Valve}}, k=0.0)
        # d (mm)
        # l (m)

        d = d / 1000
        A = π * (d / 2)^2
        I = fluid.ρ * l / A
        C = A * l / (fluid.ρ * fluid.a^2)

        new{typeof(other)}(d, l, ϵ, A, I, C, k, other)
    end
    function Pipe(d, l, ϵ, fluid::Fluid, other::Union{Twin,Vector{Twin}}, k=0.0)
        # d (mm)
        # l (m)

        d = d / 1000
        A = π * (d / 2)^2
        I = fluid.ρ * l / A
        C = A * l / (fluid.ρ * fluid.a^2)

        new{typeof(other)}(d, l, ϵ, A, I, C, k, other)
    end

end

struct System{T}
    fluid::Fluid
    pipe::T
    shaft::Shaft
    pump::Impeller
    conv::ConvConst
    function System(fluid, pipe, shaft, pump)
        new{typeof(pipe)}(fluid, pipe, shaft, pump, ConvConst())
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
    return Dict(key => nestedfield(x, key) for key ∈ fieldnames(typeof(x)))
end

function nestedfield(x, key)
    types = [ConvConst, Fluid, Impeller, Shaft, Valve, Twin, Pipe, Inputs]
    # Get the field value
    val = getfield(x, key)
    # Check if the val is a simulation type
    if typeof(val) in types
        return type2dict(val)
    else
        return val
    end
end


function type2dict(x::Vector)
    return [type2dict(i) for i in x]
end

function sim2dict(x::System)
    return Dict(
        :fluid => type2dict(x.fluid),
        :pipe => type2dict(x.pipe),
        :shaft => type2dict(x.shaft),
        :pump => type2dict(x.pump),
        :const => type2dict(ConvConst()),
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
