@connector function Power(; name, effort = 0.0, flow = 0.0, rename_var = Dict())
    sts = @variables e(t) = effort f(t) = flow
    sts_name = Dict()
    if haskey(rename_var, :e)
        merge!(sts_name, Dict(:e => Sym{FnType{Tuple{Number},Number}}(rename_var[:e])(t)))
    end
    if haskey(rename_var, :f)
        merge!(sts_name, Dict(:f => Sym{FnType{Tuple{Number},Number}}(rename_var[:f])(t)))
    end
    if !isempty(sts_name)
        ODESystem(Equation[], t, sts, []; name = name, var_to_name = sts_name)
    else
        ODESystem(Equation[], t, sts, []; name = name)
    end
end


function Mass(; name, m = 1.0, u = 0.0, rename_var = Dict())
    @named power = Power(flow = u, rename_var = rename_var)
    @unpack e, f = power
    ps = @parameters I = m
    eqs = [D(f) ~ e / I]
    keepparameters!(rename_var)
    extend(ODESystem(eqs, t, [], ps; name = name, var_to_name = rename_var), power)
end
