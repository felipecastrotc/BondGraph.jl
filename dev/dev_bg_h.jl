using DifferentialEquations, DiffEqSensitivity, ModelingToolkit
using Symbolics, Symbolics.Latexify
using LinearAlgebra, Statistics
using Flux
using Flux.Optimise: ADAM, update!
using Plots, Printf

include("lib_bg.jl")

# =============================================================================
# Definitions

m = 1.0
X = 1.
k = 5.
c = 0.2
@named mass = Mass(m=m)
@named spring = Spring(k=k, x=X)
@named spring3 = Spring3(k=k, x=0.2)
@named damper = Damper(c=c)

# =============================================================================
# 1DOF

@named model_1dof = Junction1(spring, mass, damper, couple=false)
@named model_1dof = Junction1(spring3, mass, damper, couple=false)

sys_1dof = structural_simplify(model_1dof)
equations(sys_1dof)
states(sys_1dof)
parameters(sys_1dof)

latexify(sys_1dof)

prob_1dof = ODEProblem(sys_1dof, [], (0., 30.))
sol_1dof = solve(prob_1dof)

plot(sol_1dof)

# =============================================================================
# 2DOF

@named sd = Junction1(spring, damper)

M2 = [Junction0(sd, mass; name=:b3j0)];
push!(M2, Junction1(mass; name=:b2j1, subsys=M2[1]));
push!(M2, Junction0(sd; name=:b2j0, subsys=M2[2], couple=false));

@named model_2dof = ODESystem(equations(M2), t)
sys_2dof = structural_simplify(model_2dof)
equations(sys_2dof)

prob_2dof = ODEProblem(sys_2dof, [], (0., 20.))
sol_2dof = solve(prob_2dof)

plot(sol_2dof)

# =============================================================================
# 4DOF

@named sd = Junction1(spring, damper)

M4 = [Junction0(sd, mass; name=:b3j0)];
push!(M4, Junction1(mass; name=:b2j1, subsys=M4[1]));
push!(M4, Junction0(sd; name=:b2j0, subsys=M4[2]));
push!(M4, Junction1(mass; name=:b1j1, subsys=M4[3]));
push!(M4, Junction0(sd; name=:b1j0, subsys=M4[4]));
push!(M4, Junction1(mass, spring, damper; name=:b0, subsys=M4[5], couple=false));

@named model_4dof = ODESystem(equations(M4), t)

sys_4dof = structural_simplify(model_4dof)
equations(sys_4dof)
states(sys_4dof)

prob_4dof = ODEProblem(sys_4dof, [], (0., 20.))
sol_4dof = solve(prob_4dof)

plot(sol_4dof)


# =============================================================================
# Neural

# Not currently working on automatic mode only on explicit function definition
# m = Chain(Dense(1, 30), Dense(30, 1))
# m = Flux.f64(m)
# p, re = Flux.destructure(m);
# h(x, p) = dot(re(p)([x]), [1])

h(x, p) = x*p
@register h(x, p)

function SpringU(;name, x = 0., P=1.)
    @named power = Power()
    @unpack e, f = power

    @variables q(t)=x
    ps = @parameters p=P
    
    eqs = [
            e ~ h(q, p)
            D(q) ~ f
        ]
    extend(ODESystem(eqs, t, [q], ps; name=name), power)
end

@named springU = SpringU(x=1.)

# =============================================================================
# 1DOF

@named model_1dofu = Junction1(mass, springU, damper, couple=false)

sys_1dofu = structural_simplify(model_1dofu)
equations(sys_1dofu)
parameters(sys_1dofu)

prob_1dofu = ODEProblem(sys_1dofu, [], (0., 20.))
sol_1dofu = solve(prob_1dofu)
plot(sol_1dofu)

# Set the initial conditions
u0 = prob_1dofu.u0

# Declare a function with p
function pred_ode(p, array=true)
    ps = [1, p, 1]
    sol = solve(prob_1dofu,Tsit5(),u0=u0,p=ps,saveat=t_rng)
    if array
        Array(sol)
    else
        sol
    end
end

t_rng = 0:0.1:5
y = Array(sol_1dof(t_rng))
p = 0.1

plot(y')
plot!(pred_ode(1)')

# The x was changed due the position of the variables
loss(x, y) = mean((circshift(x, 1) .- y).^2);

# Train
a = Animation()
it = 300;
hist = zeros(it);
for i in 1:it
    g = gradient((p) -> loss(pred_ode(p), y), p);
    p += -g[1];     # Simple gradient descent
    hist[i] = loss(pred_ode(p), y);
    @printf("It: %d - loss: %.3e \n", i, hist[i]);
    
    # Animation
    plt = plot(pred_ode(p, false))
    plt = plot(plt, sol_1dof(t_rng), xlabel="Time (s)", ylabel="Amplitude ()")
    frame(a, plt)
end

gif(a)

# =============================================================================
# 4DOF

@named sd = Junction1(springU, damper)

MU = [Junction0(sd, mass; name=:b3j0)];
push!(MU, Junction1(mass; name=:b2j1, subsys=MU[1]));
push!(MU, Junction0(sd; name=:b2j0, subsys=MU[2]));
push!(MU, Junction1(mass; name=:b1j1, subsys=MU[3]));
push!(MU, Junction0(sd; name=:b1j0, subsys=MU[4]));
push!(MU, Junction1(mass, springU, damper; name=:b0, subsys=MU[5], couple=false));

@named model_4dofu = ODESystem(equations(MU), t)

sys_4dofu = structural_simplify(model_4dofu)
equations(sys_4dofu)
parameters(sys_4dofu)

prob_4dofu = ODEProblem(sys_4dofu, [], (0., 20.))
sol_4dofu = solve(prob_4dofu)
plot(sol_4dofu)

# Set the initial conditions
u0 = prob_4dofu.u0

# Declare a function with p
function pred_ode(p, array=true)
    # ps = [p[1], 1.0, 1.0, 1.0, p[2], 1.0, 1.0, p[3], 1.0, 1.0, p[4], 1.0]
    ps = [p[1], 0.2, 1.0, 1.0, p[2], 0.2, 1.0, p[3], 0.2, 1.0, p[4], 0.2]
    # ps = [p, 1.0, 1.0, 1.0, p, 1.0, 1.0, p, 1.0, 1.0, p, 1.0]
    sol = solve(prob_4dofu,Tsit5(),u0=u0,p=ps,saveat=t_rng)
    if array
        Array(sol)
    else
        sol
    end
end

t_rng = 0:0.1:5
y = Array(sol_4dof(t_rng))
p = [1, 8, 15, 3]

plot(y')
plot!(pred_ode(p)')

loss(x, y) = mean((x .- y).^2);

# Train
a = Animation()
it = 600;
hist = zeros(it);
for i in 1:it
    g = gradient((p) -> loss(pred_ode(p), y), p);
    p += -2*g[1];     # Simple gradient descent
    hist[i] = loss(pred_ode(p), y);
    @printf("It: %d - loss: %.3e \n", i, hist[i]);
    
    # Animation
    plt = plot((pred_ode(p)')[:, 3:4])
    plt = plot(plt, (y')[:, 3:4], xlabel="Time (s)", ylabel="Amplitude ()")
    frame(a, plt)
end

gif(a)


# =============================================================================
# Manual 1DOF

# Declare the differential equation as a function
function pODE(du,u,p,t)
    m = 1.0
    c = 0.2
    du[1] = (m^-1)*(-h(u[2], p) - (c*u[1]))
    du[2] = u[1]
end


# Set the neural network
m = Chain(Dense(1, 50, relu), Dense(50, 1))
p, re = Flux.destructure(m);
h(x, p) = re(p)([x])[1]
# h(x, p) = p[1]*x^3 + p[2]*x^2 + p[3]*x
# p = rand(3)

# Configure the simulation
t_f = 10.
tspan = (0., t_f*2);
u0 = [0.0, 1.0];

prob_1dofm = ODEProblem(pODE, u0, tspan)

function predict_nn_ode()
  Array(solve(prob_1dofm,Tsit5(),u0=u0,p=p,saveat=t_rng))
end

# Data settings
t_rng = 0:0.2:10
y = circshift(Array(sol_1dof(t_rng)), 1);
plot(predict_nn_ode()')
plot!(y')


ps = Flux.params(p);

# Initialise optimiser
# opt = AMSGrad();
opt = ADAM(0.1);

loss(x, y) = mean((x .- y).^2);

# Train
a = Animation()
it = 100;
hist = zeros(it);
for i in 1:it
    ∇loss = gradient(() -> loss(predict_nn_ode(), y), ps);
    
    update!(opt, ps, ∇loss);
    hist[i] = loss(predict_nn_ode(), y);
    @printf("It: %d - loss: %.3e \n", i, hist[i]);

    # Animation
    plt = plot(predict_nn_ode()')
    plt = plot(plt, y', xlabel="Time (s)", ylabel="Amplitude ()")
    frame(a, plt)
end

gif(a)

plot(solve(prob_1dofm,Tsit5(),u0=u0,p=p,saveat=0:0.1:20))
plot!(sol_1dof)