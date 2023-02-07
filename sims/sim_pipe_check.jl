using DifferentialEquations
using Plots, Printf
using MAT, YAML

include("lib_types.jl")
include("lib_sim.jl")
include("lib_in_func.jl")
include("utils.jl")

conv = ConvConst()

oil = Fluid(1290, 882.32, 0.208)

d_pipe = 3 * 25.4         # (mm) Pipe diameter
ϵ_pipe = 4.6e-05          # (m) Pipe rugosity
n_segs = 4
pipe_seg = Pipe(d_pipe, 1 / n_segs, ϵ_pipe, oil)

function pipe8!(du, u, p, t)
    pipe, fluid = p
    P_i, P_o = 2 * fluid.ρ * 9.81, 0.0

    n_seg = n_segs * 2
    Q = view(u, 1:2:n_seg)
    P = view(u, 2:2:n_seg)
    dQ = view(du, 1:2:n_seg)
    dP = view(du, 2:2:n_seg)

    dQ[1] = pipe_cheng(t, [P_i, P[1], Q[1]], pipe, fluid, conv)
    dP[1] = pipe_p(t, [Q[1], Q[2]], pipe, fluid, conv)
    for i in 2:(length(P)-1)
        dQ[i] = pipe_cheng(t, [P[i-1], P[i], Q[i]], pipe, fluid, conv)
        dP[i] = pipe_p(t, [Q[i], Q[i+1]], pipe, fluid, conv)
    end
    dQ[end] = pipe_cheng(t, [P[end], P_o, Q[end]], pipe, fluid, conv)
    dP[end] = pipe_p(t, [Q[end-1], Q[end]], pipe, fluid, conv)
end

tspan = (0.0, 0.1)

# u0 = zeros(2 * 2)
u0 = zeros(n_segs * 2)

prob = ODEProblem(pipe8!, u0, tspan, [pipe_seg, oil])
sol = solve(prob, reltol=1e-8, abstol=1e-8)

t_0 = 0.0
Δt = 1 / 1290 / 100
t = collect(t_0:Δt:tspan[2])
t = collect(t_0:Δt:2/1290)
t = collect(t_0:Δt:40/1290)
nsol = sol(t)

p = plot(t, nsol[1, :])
for i in 1:2:n_segs*2
    plot!(t, nsol[i, :])
end
plot!([1 / 1290, 1 / 1290], [0, 1.5], color=:black, legend=false)

p = plot(t, nsol[2, :])
for i in 2:2:n_segs*2
    plot!(t, nsol[i, :])
end
plot!([1 / 1290, 1 / 1290], [0, 2e4], color=:black, legend=false)
p

anim = @animate for i ∈ 1:size(nsol, 2)
    plot([1/n_segs:1/n_segs:1], nsol[1:2:2*n_segs, i], legend=false)
    ylims!(0, 6.5)
    ylabel!("Velocity (m/s)")
    xlabel!("Pipe (m)")
    title!(@sprintf("Time: %f s", t[i]))
end
gif(anim, "anim_4-seg.gif", fps=15)

i = 50
plot([1/n_segs:1/n_segs:1], nsol[1:2:2*n_segs, i], legend=false)
ylims!(0, 1)
ylabel!("Velocity (m/s)")
xlabel!("Pipe (m)")
title!(@sprintf("Time: %f s", t[i]))
savefig("vel_4-seg_i-50.png")