using BondGraph
#  Plots, Symbolics.Latexify
import BondGraph: t, D
# using DifferentialEquations, JLD2

function addsubsys(sys, pss)

    ps = pss isa BgODESystem ? [pss] : collect(pss)

    if BgODESystem in typeof.(ps)

        eqs, names = Equation[], Symbol[]

        rename = Dict()

        filter!(x -> x isa BgODESystem, ps)

        for p in ps
            for s in p.systems
                # if !(s.name in names)
                push!(names, s.name)
                eqs = vcat(eqs, equations(s))
                # end
                pname = string(p.name) * "₊"
                sname = string(s.name) * "₊"
                pname = pname * sname

                for st in s.states
                    str_st = string(st)
                    sym_st = split(str_st, "(")[1]
                    rename[pname*str_st] = Symbol(sname * sym_st)
                end
            end
        end

        sys = renamevars(sys, rename)
        # Remove duplicate equations
        eqs = collect(Set(vcat(equations(sys), eqs)))
        return ODESystem(eqs; name = sys.name)
    else
        return sys
    end
end

# Pipe
@named pI = Mass(m = 1)
@named pR = Damper(c = 1)
@named pJ1 = Junction1(-pR, -pI)

equations(pJ1)
# Impeller
@parameters RL, RT

# Impeller system
@named iI = Mass(m = 1.0)
@named iRL = Damper(c = 10)

# Impeller recirculation
@named iRR = Damper(c = 1000)
@named iRes = Junction1(-iRR)
@named iIn = Junction0(pJ1, subsys = [iRes])
# @named iIn = ODESystem(vcat(equations(iIn), equations(iRes)))
@named iOut = Junction0(-pJ1, subsys = [-iRes])
# @named iOut = ODESystem(vcat(equations(iOut), equations(iRes)))
@named sys = Junction1(-iI, -iRL, iIn, -iOut, couple = false);


pss = [-iI, -iRL, iIn, -iOut]

ps = pss isa BgODESystem ? [pss] : collect(pss)

eqs, con, subsys = Equation[], Dict{Symbol,Vector{Any}}(), Dict{Symbol,Any}()

rename = Dict()

filter!(x -> x isa BgODESystem, ps);

# p = ps[4]
# s = p.subsys[1]

for p in ps
    for s in p.subsys

        if !(s.name in keys(subsys))
            con[s.name] = [p]
            subsys[s.name] = s
            eqs = vcat(eqs, equations(s))
        else
            push!(con[s.name], p)
        end
    end
end

k = :iRes
v = con[k];
# for (k, v) in con

if length(v) > 1


else

end

# end
rename




pname = string(p.name) * "₊"
sname = string(s.name) * "₊"
pname = pname * sname

for st in s.states
    str_st = string(st)
    sym_st = split(str_st, "(")[1]
    rename[pname*str_st] = Symbol(sname * sym_st)
end
