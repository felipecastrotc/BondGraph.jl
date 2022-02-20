using Plots, LinearAlgebra

function getviews(Xᵢ)
    return view(Xᵢ, 1:length(Xᵢ)-2), view(Xᵢ, 3:length(Xᵢ))
end

function mocqhsteady!(Q, H, tᵢ, A, a)
    # Create a view from the current simulation time
    Qᵢ, Hᵢ = view(Q, tᵢ, :), view(H, tᵢ, :)
    Q₁, H₁ = view(Q, tᵢ + 1, :), view(H, tᵢ + 1, :)

    # Get A and B views
    Qa, Qb = getviews(Qᵢ)
    Ha, Hb = getviews(Hᵢ)
    
    R = f*(2 * d * A)
    Ca = g * A / a
    # Positive characteristic line
    Cp = Qa + Ca * Ha - R * Δt * Qa .* abs.(Qa)
    # Negative characteristic line
    Cn = Qb - Ca * Hb - R * Δt * Qb .* abs.(Qb)
    # Estimate Qp
    Qp = 0.5 .* (Cp + Cn)
    # From Qp = Cp - Ca*Hp or Qp = Cn + Ca*Hp estimate Hp
    Hp = (Cp - Qp) ./ Ca

    # Update the Q and H matrices
    Q₁[2:length(Q₁)-1] = Qp
    H₁[2:length(H₁)-1] = Hp
end

function bc_up_clevel!(Q, H, tᵢ, Hr, A, a; k = 0.0)
    # Create a view from the current simulation time
    Qᵢ, Hᵢ = view(Q, tᵢ, :), view(H, tᵢ, :)
    # Get Q and H of the B point
    Qb, Hb = Qᵢ[2], Hᵢ[2]
    R = f * (2 * d * A)

    # Get the boundary Q
    Ca = g * A / a
    Cn = Qb - Ca * Hb - R * Δt * Qb .* abs.(Qb)

    if k == 0.0
        # Negligible loss at entrance
        Qp = Cn + Ca * Hr
        Hp = Hr
    else
        # Non-egligible loss at entrance
        k1 = Ca * (1 + k) / (2 * g * (A^2))
        Qp = (-1 + sqrt(1 + 4 * k1 * (Cn + Ca * Hr))) / (2 * k1)
        Hp = Hr - (1 + k) * (Qp^2) / (2 * g * (A^2))
    end
    # Apply the BC
    Qᵢ[1] = Qp
    Hᵢ[1] = Hp
end

function bc_down_clevel!(Q, H, tᵢ, Hr, A, a; k = 0.0)
    # Create a view from the current simulation time
    Qᵢ, Hᵢ = view(Q, tᵢ, :), view(H, tᵢ, :)
    # Get Q and H of the B point
    Qa, Ha = Qᵢ[length(Qᵢ)-1], Hᵢ[length(Qᵢ)-1]
    R = f * (2 * d * A)

    # Get the boundary Q
    Ca = g * A / a
    Cp = Qa + Ca * Ha - R * Δt * Qa .* abs.(Qa)

    if k == 0.0
        # Negligible loss at entrance
        Qp = Cp - Ca * Hr
        Hp = Hr
    else
        # Non-egligible loss at entrance
        k1 = Ca * (1 - k) / (2 * g * (A^2))
        Qp = (1 - sqrt(1 - 4 * k1 * (Cp - Ca * Hr))) / (2 * k1)
        Hp = Hr - (1 - k) * (Qp^2) / (2 * g * (A^2))
    end
    # Apply the BC
    Qᵢ[length(Qᵢ)] = Qp
    Hᵢ[length(Qᵢ)] = Hp
end

function bc_down_clevel!(Q, H, tᵢ, Hr, A, a; k = 0.0)
    # Create a view from the current simulation time
    Qᵢ, Hᵢ = view(Q, tᵢ, :), view(H, tᵢ, :)
    # Get Q and H of the B point
    Qa, Ha = Qᵢ[length(Qᵢ)-1], Hᵢ[length(Qᵢ)-1]
    
    R = f * (2 * d * A)
    # Get the boundary Q
    Ca = g * A / a
    Cp = Qa + Ca * Ha - R * Δt * Qa .* abs.(Qa)

    if k == 0.0
        # Negligible loss at entrance
        Qp = Cp - Ca * Hr
        Hp = Hr
    else
        # Non-egligible loss at entrance
        k1 = Ca * (1 - k) / (2 * g * (A^2))
        Qp = (1 - sqrt(1 - 4 * k1 * (Cp - Ca * Hr))) / (2 * k1)
        Hp = Hr - (1 - k) * (Qp^2) / (2 * g * (A^2))
    end
    # Apply the BC
    Qᵢ[length(Qᵢ)] = Qp
    Hᵢ[length(Qᵢ)] = Hp
end


# Pipe
d = 0.1
L = 1
A = π * (d^2) / 4
# Fluid properties
# R = 9e4
f = 0.01
# General properties
g = 9.81
# Mesh properties
a = 1300
Δx = 0.01
Δt = Δx /(2*a)

# a = Δx / Δt
a * Δt / Δx
# a = 1395.0
ts = 0.5

x = 0:Δx:L
t = 0:Δt:ts-Δt

# Initialize volumetric flow-rate and head matrices
Q = zeros(length(t) + 1, length(x))
H = zeros(length(t) + 1, length(x))

# HR = vcat(range(10, 3, Int(length(t) / 2)), ones(Int(length(t) / 2)) .* 3)
# HR = HR .* exp.(-100 * t)
# HR = HR .* exp.(-10 * t)
# HR = 1 .+ exp.(-300.0 .* t) + 1.0 .* cos.(2 * π * 500.0 * t) .* exp.(-300.0 .* t)
# HR = 1 .+ exp.(-30.0 .* t) + 1.0 .* cos.(2 * π * 5.0 * t) .* exp.(-3.0 .* t) .+ 10
HR = zeros(length(t)) .+ 5
plot(HR)

HF = 3
H[1, :] .= HF;
for (tᵢ, v) in enumerate(t)
    mocqhsteady!(Q, H, tᵢ, A, a)
    bc_up_clevel!(Q, H, tᵢ + 1, HR[tᵢ], A, a; k =0.1)
    bc_down_clevel!(Q, H, tᵢ + 1, HF, A, a; k = 0)
end
plot(t * a / L, H[1:end-1, Int((size(H, 2) - 1) / 2)], xlabel = "Time (s)", ylabel = "Head (m)")

# ylim = (min(minimum(H) * 0.9, 0), maximum(H) * 1.1)
# args = Dict(
#     :ylim => ylim,
#     :label => "Pressure (m)",
#     :xlabel => "Length (m)",
#     :ylabel => "Head (m)",
# )
# @gif for h in eachrow(H)
#     plot(x, h; args...)
# end every 1000

H[end-1, end-1]
# plot(t * a / L, H[1:end-1, Int((size(H, 2) - 1) / 2)], xlabel = "Time (s)", ylabel = "Head (m)")

# ylim = (min(minimum(Q) * 0.9, 0), maximum(Q) * 1.1)
# args = Dict(
#     :ylim => ylim,
#     :label => "Pressure (m)",
#     :xlabel => "Length (m)",
#     :ylabel => "Head (m)",
# )
# @gif for q in eachrow(Q)
#     plot(x, q; args...)
# end every 50

# plot(t * a / L, Q[1:end-1, Int((size(Q, 2) - 1) / 2)] ./ A, xlabel = "Time (t*a/L)", ylabel = "V (m/s)")

# Q2 = copy(Q)
# t2 = copy(t)
plot(t * a / L, Q[1:end-1, end-1] ./ A, xlabel = "Time (t*a/L)", ylabel = "V (m/s)")
plot!(t2 * a / L, Q2[1:end-1, end-1] ./ A, xlabel = "Time (t*a/L)", ylabel = "V (m/s)")

maximum(Q ./ A)

Q[end-1, end-1] / A

hf = HR[end] - HF
v = sqrt(hf * 2 * g * d /(L*f))
1000*v*d/1e-5


mul!