# Transformer and gyrator, through-oriented: port 1 in, port 2 out (Karnopp §3.2).
# A number gives the plain two-port, with the modulus as a parameter; an expression
# gives the modulated one (MTF, MGY). An expression of another subsystem's variables
# needs `@named`, which scopes those variables to the parent model.

"""
    TF(; name, m)

Transformer with modulus `m`: `e₁ = m e₂`, `f₂ = m f₁` (Karnopp Eq. 3.9). Port `p1` is
the bond into it, `p2` the bond out. `m` can be an expression (modulated transformer,
MTF); build it with `@named` when the expression uses other subsystems' variables.
"""
@twoport function TF(; name, m)
    if !(m isa Num)   # a constant modulus is a parameter
        @parameters m = m
    end
    [e1 ~ m * e2, f2 ~ m * f1]
end

"""
    GY(; name, r)

Gyrator with modulus `r`: `e₁ = r f₂`, `e₂ = r f₁` (Karnopp Eq. 3.11). Port `p1` is the
bond into it, `p2` the bond out. `r` can be an expression (modulated gyrator, MGY);
build it with `@named` when the expression uses other subsystems' variables.
"""
@twoport function GY(; name, r)
    if !(r isa Num)
        @parameters r = r
    end
    [e1 ~ r * f2, e2 ~ r * f1]
end
