# Patch path: ./src/systems/diffeqs/odesystem.jl

# var_to_name functions
import ModelingToolkit: scalarize, value, process_variables!
import ModelingToolkit: collect_var_to_name!, RefValue, nameof
import ModelingToolkit: SymbolicContinuousCallbacks, ODESystem
import ModelingToolkit: _merge, todict

# var_to_name -> modifications
const EMPTY_TGRAD = Vector{Num}(undef, 0)
const EMPTY_JAC = Matrix{Num}(undef, 0, 0)


function ODESystem(
    deqs::AbstractVector{<:Equation}, iv, dvs, ps;
    controls = Num[],
    observed = Equation[],
    systems = ODESystem[],
    name = nothing,
    default_u0 = Dict(),
    default_p = Dict(),
    defaults = _merge(Dict(default_u0), Dict(default_p)),
    connector_type = nothing,
    preface = nothing,
    continuous_events = nothing,
    checks = true,
    var_to_name = Dict()
)
    name === nothing && throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
    deqs = scalarize(deqs)
    @assert all(control -> any(isequal.(control, ps)), controls) "All controls must also be parameters."

    iv′ = value(scalarize(iv))
    dvs′ = value.(scalarize(dvs))
    ps′ = value.(scalarize(ps))
    ctrl′ = value.(scalarize(controls))

    if !(isempty(default_u0) && isempty(default_p))
        Base.depwarn("`default_u0` and `default_p` are deprecated. Use `defaults` instead.", :ODESystem, force = true)
    end
    defaults = todict(defaults)
    defaults = Dict{Any,Any}(value(k) => value(v) for (k, v) in pairs(defaults))

    var_to_name_all = Dict()
    process_variables!(var_to_name_all, defaults, dvs′)
    process_variables!(var_to_name_all, defaults, ps′)
    isempty(observed) || collect_var_to_name!(var_to_name_all, (eq.lhs for eq in observed))

    var_to_name = merge(var_to_name_all, var_to_name)

    tgrad = RefValue(EMPTY_TGRAD)
    jac = RefValue{Any}(EMPTY_JAC)
    ctrl_jac = RefValue{Any}(EMPTY_JAC)
    Wfact = RefValue(EMPTY_JAC)
    Wfact_t = RefValue(EMPTY_JAC)
    sysnames = nameof.(systems)
    if length(unique(sysnames)) != length(sysnames)
        throw(ArgumentError("System names must be unique."))
    end
    cont_callbacks = SymbolicContinuousCallbacks(continuous_events)
    ODESystem(deqs, iv′, dvs′, ps′, var_to_name, ctrl′, observed, tgrad, jac, ctrl_jac, Wfact, Wfact_t, name, systems, defaults, nothing, connector_type, nothing, preface, cont_callbacks, checks = checks)
end
