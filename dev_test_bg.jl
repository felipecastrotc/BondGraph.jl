using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Symbolics, Symbolics.Latexify

include("lib_bg.jl")


@variables t
D = Differential(t)


@connector function Power(; name, effort = 0.0, flow = 0.0)
    sts = @variables e(t) = effort f(t) = flow
    ODESystem(Equation[], t, sts, []; name = name)
end

# =============================================================================
# Ports

function Junction1(ps...; name, subsys = [], couple = true)
    # Check subsys type
    sys = subsys isa ODESystem ? [subsys] : subsys

    # Get connections
    con = collect(Base.Flatten([ps, sys]))
    if couple
        @named power = Power()
        @unpack e, f = power

        # Σ efforts
        eqs = [e ~ sum(c -> c.e, con)]
        # f₁ = f₂ = f₃
        if length(con) > 1
            for i in 1:(length(con)-1)
                push!(eqs, con[i].f ~ con[i+1].f)
            end
        end
        push!(eqs, con[1].f ~ f)
        # Build subsystem
        ex = extend(ODESystem(eqs, t, [], []; name = name), power)
        if length(ps) > 0
            compose(ex, ps...)
        else
            ex
        end
    else
        # Σ efforts
        eqs = [0 ~ sum(p -> p.e, con)]
        # f₁ = f₂ = f₃
        if length(con) > 1
            for i in 1:(length(con)-1)
                push!(eqs, con[i].f ~ con[i+1].f)
            end
        end
        # Build system
        ex = ODESystem(eqs, t, [], []; name = name)
        if length(ps) > 0
            compose(ex, ps...)
        else
            ex
        end
    end
end

function Junction0(ps...; name, subsys = [], couple = true)
    # Check subsys type
    sys = subsys isa ODESystem ? [subsys] : subsys

    # Get connections
    con = collect(Base.Flatten([ps, sys]))
    if couple
        @named power = Power()
        @unpack e, f = power

        # Σ flows
        eqs = [f ~ sum(c -> c.f, con)]
        # e₁ = e₂ = e₃
        if length(con) > 1
            for i in 1:(length(con)-1)
                push!(eqs, con[i].e ~ con[i+1].e)
            end
        end
        push!(eqs, con[1].e ~ e)
        # Build subsystem
        ex = extend(ODESystem(eqs, t, [], []; name = name), power)
        if length(ps) > 0
            compose(ex, ps...)
        else
            ex
        end
    else
        # Σ flows
        eqs = [0 ~ sum(p -> p.f, con)]
        # e₁ = e₂ = e₃
        if length(con) > 1
            for i in 1:(length(con)-1)
                push!(eqs, con[i].e ~ con[i+1].e)
            end
        end
        # Build system
        ex = ODESystem(eqs, t, [], []; name = name)
        if length(ps) > 0
            compose(ex, ps...)
        else
            ex
        end
    end
end

function keepparameters!(rename_var)
    filter!(p -> p.first != :e && p.first != :f, rename_var)
end

# =============================================================================
# Elements

function Mass(; name, m = 1.0, u = 0.0)
    @named power = Power(flow = u)
    @unpack e, f = power
    ps = @parameters I = m
    eqs = [
        D(f) ~ e / I
    ]
    extend(ODESystem(eqs, t, [], ps; name = name), power)
end

function Spring(; name, k = 10, x = 0.0)
    @named power = Power()
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

function Spring3(; name, k = 10, x = 0.0)
    @named power = Power()
    @unpack e, f = power

    @variables q(t) = x
    ps = @parameters C = 1 / k

    eqs = [
        e ~ q^3 / C
        D(q) ~ f
        # D(e) ~ f/C
    ]
    extend(ODESystem(eqs, t, [q], ps; name = name), power)
end

function Damper(; name, c = 10)
    @named power = Power()
    @unpack e, f = power

    ps = @parameters R = c
    eqs = [
        e ~ f * R
    ]
    extend(ODESystem(eqs, t, [], ps; name = name), power)
end

# =============================================================================
# Definitions

m = 1.0
X = 1.0
k = 10.0
c = 1.0
@named mass = Mass(m = m)
@named spring = Spring(k = k, x = X)
@named damper = Damper(c = c)


# =============================================================================
# 2DOF

@named sd = Junction1(spring, damper)

M = [Junction1(spring, damper, mass; name = :b1)];
push!(M, Junction0(sd; name = :b2j0, subsys = M[1]));
push!(M, Junction1(mass; name = :b2j1, subsys = M[2], couple = false));

model = ODESystem(equations(M), t; name = :tst)
sys = structural_simplify(model)
equations(sys)

prob = ODEProblem(sys, [], (0.0, 20.0))
sol = solve(prob)
plot(sol)

# Alternative scheme - REVISAR PQ ESTÀ DANDO ERRADO!!!
M2 = [Junction1(spring, damper, mass; name = :b1)];
push!(M2, Junction1(mass; name = :b2j1));
push!(M2, Junction0(sd, M2[1], M2[2]; name = :b2j0, couple = false));

model2 = ODESystem(equations(M2[3]), t; name = :tst)
sys2 = structural_simplify(model2)
equations(sys2)

equations(sys)

prob2 = ODEProblem(sys2, [], (0.0, 20.0))
prob2
sol2 = solve(prob2)
plot(sol2)

sys2.states
sys.states

Array(sol2)[1, :] - Array(sol)[3, :]
Array(sol2)[2, :] - Array(sol)[1, :]
Array(sol2)[3, :] - Array(sol)[2, :]
Array(sol2)[4, :] + Array(sol)[4, :]

# =============================================================================
# Dev
# Debug



sys.states
sys.ps

parameters(sys)
states(sys)

Symbol(getname())


# =============================================================================
# Dev
# replace names

@named m = Mass(m = 1.0)
@named k = Spring(k = 10, x = 1.0)
@named c = Damper(c = 1.0)
@named s = Junction1(c, k)

@named sd = Junction1(k, c)
@named a = Junction1(k, c, m);
@named b1 = Junction1(m);
@named b2 = Junction0(s, a, b1, couple = false);

sys = structural_simplify(b2)

equations(mdl)

# Add alias names

# Mass
u = 1.0
m = 10.0
name = :tst

@named m = Mass(m = 1.0, rename_var = Dict(:f => :ẋ₁))

equations(m)

@named mdl = Junction1(m, k, couple = false)

# =============================================================================
# Dev
# rename by Dict()

eqs = equations(sys)
O = eqs[2].rhs


# function renamexpr(eq::Equation, sys)
#     _lhs = handleterms(eq.lhs, sys)
#     _rhs = handleterms(eq.rhs, sys)
#     _lhs ~ _rhs
# end

# function handleterms(O, sys) where {T}
# based on namespace_expr function from MTK
ivs = independent_variables(sys)
O = unwrap(O)
if any(isequal(O), ivs)
    return O
elseif isvariable(O)
    renamespace(sys, O)
elseif istree(O)
    renamed = map(a -> handleterms(a, sys), arguments(O))
    if symtype(operation(O)) <: FnType
        renamespace(sys, O)
    else
        similarterm(O, operation(O), renamed)
    end
elseif O isa Array
    map(Base.Fix2(namespace_expr, sys), O)
else
    O
end
# end


vars = vcat(parameters(sys), states(sys))
cvt = []
dct =
    x = arguments(O)[1]
x = O
# function renamespace(sys, x)
x = unwrap(x)
if x isa Symbolic
    let scope = getmetadata(x, SymScope, LocalScope())
        if scope isa LocalScope
            sys_name = getname(sys)
            var_name = getname(x)

            rename(x, renamespace(sys_name, var_name))
        elseif scope isa ParentScope
            setmetadata(x, SymScope, scope.parent)
        else # GlobalScope
            x
        end
    end
else

end
# end

Symbol()

