using SymPy
using Latexify

ω, Q, Q_s, L, A, d, μ, ρ, γa, γb, kf, ks, kl = symbols("ω, Q, Q⁺, L, A, d, μ, ρ, γ_a, γ_b, k_f, k_s, k_l")

πs, β = symbols("π, β")

# Support variables
A = πs * (d / 2)^2
Re = ρ * Q * d / (μ * A)
darcy = ((ρ * L) / (2 * d))

# ------------------------------------------------------------------------------
# Friction
dp_f = kf * (64 / Re) * darcy * ((Q / A)^2)

# Simplify
k_f1 = symbols("k_f1")

# with sub
dp_f2 = subs(dp_f, 128 * L * k_f / (πs * d^4) => k_f1)
string(dp_f2)

# ------------------------------------------------------------------------------
# Euler
dp_eu = expand((ρ) * (γa * Q - γb * ω) * ω)
string(dp_eu)

# ------------------------------------------------------------------------------
# Q shock
dp_s = expand(-ks * ρ * (d / 2)^2 * ω^2 * (1 - Q / Q_s)^2)

# Simplify
k_s1, k_s2, k_s3 = symbols("k_s1, k_s2, k_s3")

# with sub
dp_s2 = subs(dp_s, ks * d^2 / (4 * Q_s^2) => k_s1, ks * d^2 / (2 * Q_s) => k_s2, ks * d^2 / 4 => k_s3)

# ------------------------------------------------------------------------------
# Loc
dp_l = kl * ρ * ((Q / A)^2)

# Simplify
k_l1 = symbols("k_l1")

# with sub
dp_l2 = subs(dp_l, 16 * kl / (πs^2 * d^4) => k_l1)

# ------------------------------------------------------------------------------
# New shock

dp_s = expand(ρ * ks * ((ω * (d / 2) - Q / (A * β))^2))

# Simplify
k_s1, k_s2, k_s3 = symbols("k_s1, k_s2, k_s3")

# with sub
dp_s2 = subs(dp_s, (16 * ks) / (d^4 * β^2 * πs^2) => k_s1, (4 * ks) / (d * β * πs) => k_s2, ks * d^2 / 4 => k_s3)

string(dp_s2)


# ------------------------------------------------------------------------------
# Disk Friction loss

# Reynolds number used on disk friction loss
Rew = ρ * (ω * d / 2) * d / μ

krr, r, s, δ, R_n = symbols("k_rr, r, s, δ, R_n")

# For cylinder case
P_df = krr * ρ * ω^3 * (d / 2)^4 * L      # Eq. 3.64 pag 136 Gullich 2010
krr_eq = πs * (d / 2) / (2 * Rew * s)     # Eq. 3.68 pag 136 Gullich 2010

P_dff = subs(P_df, krr => krr_eq) / ω
latexify(P_dff)

-ρ * Q^2 + ρ * Q * ω - P_dff

# For disk case
P_df = (krr/cos(δ))* ρ * ω^3 * (d/2)^5*(1 - (R_n/(d/2))^5)      # Eq. 3.64 pag 136 Gullich 2010
krr_eq = πs * (d / 2) / (2 * Rew * s)     # Eq. 3.68 pag 136 Gullich 2010
P_dff = subs(P_df, krr => krr_eq)
latexify(P_dff)



