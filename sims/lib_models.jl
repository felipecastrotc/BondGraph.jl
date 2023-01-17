include("lib_flow.jl")

# -------------------------------------------------------------------
# Model manual -  1 upstream 1 downstream - Source of effort
# (H)ead Q-flow-rate (P)ipe-1 (I)mpeller-1
function hqp1i1!(du, u, p, t)
    Pi = 0.0
    Po = 0.0
    Qi, ω, Q1, Q3, P2, P4 = u[:]
    ϵ, Δs, d, L, C, ρ, μ, γa, γb, di, Δsi, Ii, ca, Ia, H = hqp1i1_unpack(p)
    τ, τp = p.input.func, p.input.p
    # Impeller
    du[1] = (P2 - P4 + (ρ / H) * (γa * Qi - γb * ω) * ω - darcyq_lam(Δs, d, μ) * Qi / H) * H / Ii
    # Shaft
    du[2] = (τ(t, τp...) - ca * ω - ρ * (γa * Qi - γb * ω) * Qi) / Ia
    # Pipeline
    # Upstream
    # Q1
    du[3] = (Pi - P2 - darcyq_lam(Δs, d, μ) * Q1 / H) * H / L
    # Downstream
    # Q3
    du[4] = (P4 - Po - darcyq_lam(Δs, d, μ) * Q3 / H) * H / L
    # Pressures
    du[5] = (Q1 - Qi) / (C * H)           # P2 Upstream inlet
    du[6] = (Qi - Q3) / (C * H)           # P4 Downstream coupling pump
end

function hqp1i1_unpack(sim)
    return [sim.sys.pipe.ϵ,
            sim.sys.pipe.Δs,
            sim.sys.pipe.d,
            sim.sys.pipe.L,
            sim.sys.pipe.C,
            sim.sys.fluid.ρ,
            sim.sys.fluid.μ,
            sim.sys.pump.γa,
            sim.sys.pump.γb,
            sim.sys.pump.di,
            sim.sys.pump.Δsi,
            sim.sys.pump.Ii,
            sim.sys.shaft.ca,
            sim.sys.shaft.Ia,
            H,
    ]
end

function hqp1i1_ps()
    return [:ϵ, :Δs, :d, :L, :C, :ρ, :μ, :γa, :γb, :di, :Δsi, :Ii, :ca, :Ia, :H]
end


