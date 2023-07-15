using BondGraph
using BondGraph: t, D
using DifferentialEquations
using ModelingToolkit
using Plots
using Symbolics.Latexify

import BondGraph: adjmtx2eqs, mGY, get_var, get_bg_junction, tpgy, tptf, op, j1, j0, equalityeqs
import BondGraph: get_sm_mtx!, bgflow, bgeffort, sumvar, get_var

@named L = Mass(m=0.5)        # (H) Electric inductance
@named R = Damper(c=1.0)    # (Ohm) Electric resistance
@named V = Se(12.0)           # (V) Voltage

@named J = Mass(m=0.01)     # (kg⋅m^2) Shaft moment of inertia
@named c = Damper(c=0.1)    # (N⋅m⋅s) Shaft viscous friction
@named τ = Se(1.0)            # (N⋅m) Applied torque on the shaft

g = 0.01    # DC motor torque constant

@named jm = Junction1([-1, c], [-1, J])    # Mechanical part
@named jm = Junction1(τ, [-1, c], [-1, J])    # Mechanical part
@named je = Junction1(V, [-1, R], [-1, L])    # Electical part

@named gy, eqs = mGY(je, jm, g=g)

equations(gy)

@named mdl = ODESystem(eqs, t)
mdl = compose(mdl, gy, je, jm)
equations(mdl)

generate_graph(mdl)

equations(expand_connections(mdl))

@named sys = simplifysys(mdl)
equations(sys)

emdl = expand_connections(mdl)
equations(emdl)
equations(alias_elimination(emdl))
sys = structural_simplify(emdl)
observed(sys)
equations(structural_simplify(emdl))

tspan = (0, 10)

# Generate an `ODEProblem` as in the ModelingToolkit, passing the simulation interval.
prob = ODEProblem(sys, [], tspan)

# Use the `solve` to simulate the system.
sol = solve(prob)
plot(sol, xlabel="Time", ylabel="Amplitude")


import ModelingToolkit: ConnectionSet, ConnectionElement
import ModelingToolkit: namespaced_var, generate_connection_set!
import BondGraph: get_bg_connection_set!, csets2dict, csets2adjmtx,get_sm_mtx!

connectionsets = ConnectionSet[]

# Generate the bond graph connection sets
sys = generate_connection_set!(connectionsets, mdl, nothing, nothing)

# Get the bond graph connection sets
bgconnectionsets = get_bg_connection_set!(connectionsets)

# Create a dictionary mapping connection set names to connection elements
str2con = csets2dict(bgconnectionsets)

am = csets2adjmtx(bgconnectionsets, str2con)


function adjmtx2eqs(am, str2con)
    # Create a dictionary mapping indices to connection element keys
    idx2k = Dict(i => k for (i, k) in enumerate(keys(str2con)))

    # Generate the signal matrix
    sm = get_sm_mtx!(am, idx2k, str2con)

    # Calculate the input and output vectors
    in_vec = sum(am, dims=1)[:]
    out_vec = sum(am, dims=2)[:]

    idx_con = 1:size(am, 1)
    in_con = idx_con[in_vec.>0]

    eqs = Equation[]
    
    for i in in_con
        # i = 16
        v = get_var(i, idx2k, str2con)

        # Get variables leaving the node
        vout = Num[]
        if out_vec[i] == 1
            push!(vout, -get_var(i, idx2k, str2con))
        elseif out_vec[i] > 1
            for j = 1:out_vec[i]
                push!(vout, -add_idx(v, j))
            end
        end

        # Get variables being added to the node
        vin = Num[]

        # Iterate over the connections
        for j in idx_con[am[:, i].>0]
            # j = 1
            # Check if the input has multiple outputs
            msk = am[j, :] .> 0
            chk = sum(msk) > 1 ? findfirst(sort(idx_con[msk]) .== i) : nothing
            jtype = get_bg_junction(get_var(j, idx2k, str2con))[1]

            if !isnothing(chk)
                # Create the variable for the i node
                push!(vin, add_idx(get_var(j, idx2k, str2con), chk))
            elseif jtype === tpgy || jtype === tptf
                push!(vin, get_var(j, idx2k, str2con))
            elseif jtype === op
                # Apply the signal matrix
                sgn = sm[j, i]
                var = get_var(j, idx2k, str2con)
                # The negative direction in BG refers to the product e*f
                sgn = get_bg_junction(var)[2] === bgflow ? abs(sgn) : sgn
                push!(vin, sgn * var)
            else
                push!(vin, get_var(j, idx2k, str2con))
            end
        end

        # Apply the junction type
        jtype = get_bg_junction(v)[1]
        vtype = get_bg_junction(v)[2]

        if (jtype === j0 && vtype === bgeffort) || (jtype === j1 && vtype === bgflow)
            push!(eqs, equalityeqs(vcat(vout, vin))...)
        elseif (jtype === j0 && vtype === bgflow) || (jtype === j1 && vtype === bgeffort)
            push!(eqs, sumvar(vcat(vout, vin)))
        elseif (jtype === tpgy) || (jtype === tptf)
            lvin, lvout = length(vin), length(vout)

            if (lvin <= 1) && (lvout <= 1) && (lvout + lvout == 1)
                throw(
                    DomainError(
                        "The $jtype type only allows one connection in and one connection out",
                    ),
                )
            end

            if lvin == 1
                # The negative direction in BG refers to the product e*f
                sgn = get_bg_junction(vin[1])[2] === bgflow ? -1 : 1
                push!(eqs, equalityeqs(vcat(sgn*vin, v))...)
            elseif lvout == 1
                push!(eqs, equalityeqs(vcat(vout, v))...)
            end
        end
    end

    return eqs
end

# ----------------------------------------------------------------------
# Check simple

m = 10
@named mass = Mass(m=m)
@named damper = Damper(c=m)
@named spring = Spring(k=m)
@named f = Se(cos(t))

@named mdl = Junction1([-1, mass], [-1, damper], [-1, spring], f)
# @named mdl = Junction1([1, mass], [1, damper], [1, spring], [-1, f])
generate_graph(mdl)
equations(mdl)
equations(expand_connections(mdl))

@named sys = simplifysys(mdl)
equations(sys)

# ----------------------------------------------------------------------
# Check Wiley
# https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781119958376.app1

using DelimitedFiles

# Read the CSV file into a matrix
matrix = readdlm("ref_wiley.csv", ',')
sim20 = convert(Array{Float64}, matrix[2:end, :])

@parameters V, g, n, ra, J_1, J_2, B, τ

V = 1.0
g = 1.0
n = 1.0
ra = 1.0
J_1 = 1.0
J_2 = 1.0
B = 1.0
τ = 1.0

# Electric
@named v = Se(V)
@named Ra = Damper(c=ra)
# M1
@named J1 = Mass(m=J_1)
@named Rb = Damper(c=B)
# M2
@named J2 = Mass(m=J_2)
@named T = Se(τ)
# MP
@named je = Junction1(v, [-1, Ra])
@named jm1 = Junction1([-1, Rb], [-1, J1])
@named jm2 = Junction1([-1, J2], T)

# TP
@named gy, eqs = mGY(je, jm1, g=g)
@named tf = mTF(jm1, jm2, r=n, coneqs=eqs)

@named mdl = ODESystem(eqs, t)
mdl = compose(mdl, tf, gy, je, jm1, jm2)

generate_graph(mdl)

equations(alias_elimination(expand_connections(mdl)))

@named sys = simplifysys(mdl)

equations(sys)

prob = ODEProblem(sys, [], (0, 10))
sol = solve(prob)

plot(sol, label="BGToolkit")
plot!(sim20[:, 1], sim20[:, 2:end], label="20-sim", line=(4, 0.5, :dash), linewidth=2)



# ----------------------------------------------------------------------
# Check TF simple Wikipedia

m = 10
@named i3 = Mass(m=m)
@named r2 = Damper(c=m)
@named c6 = Spring(k=m)
@named r7 = Damper(c=m)
@named f = Se(cos(t))
@parameters r

@named jl = Junction1([-1, i3], [-1, r2], f)
@named jr = Junction0([-1, c6], [-1, r7])

@named tf, eqs = mTF(jl, jr, r=1.0)

@named mdl = ODESystem(eqs, t)
mdl = compose(mdl, jl, jr, tf)

generate_graph(mdl)
equations(expand_connections(mdl))


# @named mdl = Junction1([-1, mass], [-1, damper], [-1, spring], f)
equations(mdl)

@named sys = simplifysys(mdl)
equations(sys)

latexify(sys)