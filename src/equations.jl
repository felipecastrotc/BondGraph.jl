"""
    bg_equations(sys; rename = Dict())

The equations of the model `sys` in terms of its states, as `full_equations` of the
compiled model (`sys` is compiled with `mtkcompile` unless it already is). `rename` maps
model variables or parameters to the names to show instead, e.g.
`Dict(model.C.q => x, D(model.C.q) => ẋ, model.R.R => b)`. Renaming a variable also
renames its derivative, unless the derivative has its own entry. The keys may come from
`sys` itself or from its compiled form. The model is not modified: for a renamed
variable in a solution, add an alias equation such as `x ~ model.C.q` to the model.

```julia
@variables x(t) ẋ(t)
latexify(bg_equations(model; rename = Dict(model.C.q => x, D(model.C.q) => ẋ)))
```
"""
_isderivative(x) = iscall(x) && operation(x) isa Differential
_underived(x) = _isderivative(x) ? only(arguments(x)) : x

function bg_equations(sys; rename = Dict())
    compiled = ModelingToolkit.iscomplete(sys)
    eqs = full_equations(compiled ? sys : mtkcompile(sys))
    # A key from the uncompiled model is namespaced by its name: `model₊C₊q` is `C₊q`.
    root = compiled ? "" : string(nameof(sys), "₊")
    vars = Dict(string(getname(v)) => v for eq in eqs
                for v in _underived.(get_variables(eq.lhs - eq.rhs)))
    find(k) = get(vars, chopprefix(string(getname(k)), root)) do
        throw(ArgumentError("`$k` does not appear in the model's equations"))
    end
    keys_ = unwrap.(collect(keys(rename)))
    rules = Dict{Any, Any}(operation(k)(find(_underived(k))) => unwrap(rename[k])
                           for k in filter(_isderivative, keys_))
    for k in filter(!_isderivative, keys_)
        rules[find(k)] = unwrap(rename[k])
        get!(rules, unwrap(D(find(k))), unwrap(D(rename[k])))
    end
    [substitute(eq, rules) for eq in eqs]
end
