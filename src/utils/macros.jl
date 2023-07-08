using MacroTools


# TODO: Add the support for initial states
macro oneport(ex)
    block_initial = quote
        @named power = Power(type=op)
    end

    block_midle = MacroTools.striplines(ex.args[2])

    # Generate the equation variable
    block_eqs = :(eqs = vcat([$(block_midle.args[end])]...))

    block_end = quote

        # Get all variables from the equations
        vars = vcat([collect(Set(ModelingToolkit.get_variables(eq))) for eq in eqs]...)
        # Remove the power variables
        filter!(x -> !isequal(x, power.e) && !isequal(x, power.f), vars)
        # Remove duplicates
        vars = collect(Set(vars))

        # Get the parameters and states of the system
        if length(vars) > 0
            sts = filter(x -> ~isindependent(Num(x)), vars)
            ps = filter(x -> isindependent(Num(x)), vars)
        else
            sts, ps = [], []
        end
        # Generate the compatible system
        compose(ODESystem(eqs, t, sts, ps; name=name), power)
    end

    # Create the function body
    body = Expr(:block, block_initial.args..., block_midle.args[1:end-1]..., block_eqs.args..., block_end.args...)
    # Cleanup
    body = MacroTools.striplines(body)
    # Build function
    ex.args[2] = body
    return ex
end

function setinitialval(power; effort=nothing, flow=nothing)

    if !isa(effort, Number)
        effort = ModelingToolkit.getdefault(power.e)
    end
    if !isa(flow, Number)
        flow = ModelingToolkit.getdefault(power.f)
    end

    type = get_bg_junction(power.e)[1]

    return Power(name=power.name, effort=effort, flow=flow, type=type)
end