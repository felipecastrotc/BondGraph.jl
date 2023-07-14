using BondGraph
using BondGraph: t, D
using DifferentialEquations
using ModelingToolkit
using Plots
using Symbolics.Latexify

import BondGraph: adjmtx2eqs, mGY, get_var, get_bg_junction, tpgy, tptf, op, j1, j0, equalityeqs

@named L = Mass(m=0.5)        # (H) Electric inductance
@named R = Damper(c=1.0)    # (Ohm) Electric resistance
@named V = Se(12.0)           # (V) Voltage

@named J = Mass(m=0.01)     # (kg⋅m^2) Shaft moment of inertia
@named c = Damper(c=0.1)    # (N⋅m⋅s) Shaft viscous friction
@named τ = Se(1.0)            # (N⋅m) Applied torque on the shaft

g = 0.01    # DC motor torque constant

@named jm = Junction1(τ, [-1, c], [-1, J])    # Mechanical part
@named je = Junction1(V, [-1, R], [-1, L])    # Electical part

eqs = []
@named gy = mGY(je, jm, g=g, coneqs=eqs)

equations(gy)

@named mdl = ODESystem(eqs, t)
mdl = compose(mdl, gy, je, jm)
equations(mdl)

generate_graph(mdl)

equations(expand_connections(mdl))

@named sys = simplifysys(mdl)

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
    # 6 e 16
    idx2k
    # Generate the signal matrix
    sm = get_sm_mtx!(am, idx2k, str2con)

    # Calculate the input and output vectors
    in_vec = sum(am, dims=1)[:]
    out_vec = sum(am, dims=2)[:]

    idx_con = 1:size(am, 1)
    in_con = idx_con[in_vec.>0]

    eqs = Equation[]
    # in_con
    # [get_var(i, idx2k, str2con) for i in in_con]
    for i in in_con
        # i = 7
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
            # j = idx_con[am[:, i].>0][3]
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
                push!(vin, sgn * get_var(j, idx2k, str2con))
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
                push!(eqs, equalityeqs(vcat(vin, v))...)
            elseif lvout == 1
                push!(eqs, equalityeqs(vcat(-vout, v))...)
            end
        end
    end

    return eqs
end


m = 10
@named mass = Mass(m=m)
@named f = Se(10)
@named g = Se(9.81*m)

@named mdl = Junction1(mass, f, [-1, g])
equations(expand_connections(mdl))

@named left = Junction1(mass, [-1, g])
@named right = Junction1(f)

@named tf, eqs = mTF(left, right, r = 1.0)

isnothing(tf)