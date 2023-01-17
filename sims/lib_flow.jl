# Support functions 

function AEq(d)
    return π * (d / 2)^2
end

function ReEq(ρ, q, d, μ)
    A = π * (d / 2)^2
    return (ρ * (q / A) * d) / μ
end

function darcyq_lam(Δs, d, μ)
    # using SymPy
    return 128 * Δs * μ / (pi * d^4)
end
