using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Symbolics, Symbolics.Latexify
using Flux, Statistics
using DiffEqSensitivity
using Flux.Optimise: ADAM, update!
using Printf

# =============================================================================
# Generate data

function msd(du, u, p, t)
    du[1] = (m^-1) * (-k*u[2] - (c * u[1]))
    du[2] = u[1]
end

function ms3d(du, u, p, t)
    du[1] = (m^-1) * (-k3*(u[2]^3) - (c * u[1]))
    du[2] = u[1]
end


# Parameters
m = 1.0
k = 10.0
k3 = 60.0
c = 0.9

# Simulation
tf = 10
u0 = [0.0, 1.0]
prob1 = ODEProblem(msd, u0, (0.0, tf))
prob3 = ODEProblem(ms3d, u0, (0.0, tf))

sol1 = solve(prob1, Tsit5(), u0 = u0);
sol3 = solve(prob3, Tsit5(), u0 = u0);

plot(sol1)
plot(sol3)


# =============================================================================
# Generate data

function pODE(du, u, p, t)
    du[1] = (m^-1) * (-h(u[2], p) - (c * u[1]))
    du[2] = u[1]
end

nn = Chain(Dense(1, 10, relu), Dense(10, 10, relu), Dense(10, 10, relu), Dense(10, 1))
# m = Flux.f64(m)
p, re = Flux.destructure(nn);
h(x, p) = re(p)([x])[1]


trng = 0:0.1:tf
prob = ODEProblem(pODE, u0, (0.0, tf))

function predict_n_ode()
    Array(solve(prob, Tsit5(), u0 = u0, p = p, saveat = trng))
end

plot(predict_n_ode()')

y = Array(sol3(trng));
# y = circshift(y, 1);


ps = Flux.params(p);

# Initialise optimiser
# opt = AMSGrad();
opt = ADAM(0.05);

# loss(x, y) = mean((x .- y).^2);
loss(x, y) = sum(abs2, y .- x);

# Iterate
# it = 19;
it = 1000;
hist = zeros(it);
for i = 1:it
    ∇loss = gradient(() -> loss(predict_n_ode(), y), ps)
    update!(opt, ps, ∇loss)
    hist[i] = loss(predict_n_ode(), y)
    cb()
    @printf("It: %d - loss: %.3e \n", i, hist[i])
end


cb = function (; doplot = false) # callback function to observe training
    pred = predict_n_ode()
    #   display(loss(pred, y))
    # plot current prediction against data
    pl = scatter(trng, y', label = "data")
    scatter!(pl, trng, pred', label = "prediction")
    display(plot(pl))
end
cb()
