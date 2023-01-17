using StructuralIdentifiability

ode = @ODEmodel(
	x1'(t) = k * (1 - x1(t) - x2(t)),
	x2'(t) = r * (1 - x1(t) - x2(t)),
	y1(t) = x1(t),
	y2(t) = x2(t) # new output function!
)

local_id = assess_local_identifiability(ode, 0.99) # this is a different result!


ode = @ODEmodel(
    qi'(t) = ( - 2*P4(t) + ρ* (γa * Qi(t) - γb * ω(t)) * ω(t) - (128 * Δs * μ / (3.14159 * d^4)) * Qi)/ Ii, 
    # qi'(t) = (P2(t) - P4(t) + ρ* (γa * Qi(t) - γb * ω(t)) * ω(t) - (128 * Δs * μ / (3.14159 * d^4)) * Qi)/ Ii, 
    ω'(t) = (τ(t) - ca * ω(t) - ρ * (γa * Qi(t) - γb * ω(t)) * Qi(t)) / Ia,
    # Q1'(t) = (Pi - P2(t) - (128 * Δs * μ / (3.14159 * d^4)) * Q1(t)) / L, 
    Q3'(t) = (P4(t) - Po - (128 * Δs * μ / (3.14159 * d^4)) * Q3(t)) / L,
    # P2'(t) = (Q1(t) - Qi(t)) / C,
    P4'(t) = (Qi(t) - Q3(t)) / C,
    y1(t) = P4(t),
    y2(t) = τ(t)
)

local_id = assess_local_identifiability(ode, 0.99)