using MAT, YAML

include("lib_types.jl")
include("lib_models.jl")
include("lib_in_func.jl")

PATH_SIM_STORE = "/Users/fctc/me/doutorado/dev/pinn-jax/cfg-sim/"

function store_sim(name, sol, sim, t, offset=0.0)

    # Build a dictionary with all settings used
    sim_dct = sim2dict(sim)

    # Store info
    f = open(PATH_SIM_STORE * name * "_cfg.yml", "w")
    YAML.write(f, Dict(Symbol(name) => sim_dct))
    close(f)

    # Store simulation
    nsol = sol(t)
    dnsol = sol(t, Val{1})
    t .-= offset

    file = matopen(PATH_SIM_STORE * name * "_sim.h5", "w")

    write(file, "sol", Matrix(nsol))
    write(file, "dsol", Matrix(dnsol))
    write(file, "t", collect(t))
    write(file, "dt", diff(t)[1])

    close(file)

end