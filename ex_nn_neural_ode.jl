using DifferentialEquations, DiffEqFlux, Plots



# ===============================================
# Optimization of Ordinary Differential Equations

function lotka_volterra!(du, u, p, t)
  x, y = u
  α, β, δ, γ = p
  du[1] = dx = α*x - β*x*y
  du[2] = dy = -δ*y + γ*x*y
end

# Initial condition
u0 = [1.0, 1.0]

# Simulation interval and intermediary points
tspan = (0.0, 10.0)
tsteps = 0.0:0.1:10.0

# LV equation parameter. p = [α, β, δ, γ]
p = [1.5, 1.0, 3.0, 1.0]

# Setup the ODE problem, then solve
prob = ODEProblem(lotka_volterra!, u0, tspan, p)
sol = solve(prob, Tsit5())
?Tsit5

# Plot the solution
plot(sol)
savefig("LV_ode.png")

function loss(p)
#   sol = solve(prob, Tsit5(), p=p, saveat = tsteps)
  sol = solve(prob, p=p, saveat = tsteps)
  loss = sum(abs2, sol.-1)
  return loss, sol
end

callback = function (p, l, pred)
  display(l)
  plt = plot(pred, ylim = (0, 6))
  display(plt)
  # Tell sciml_train to not halt the optimization. If return true, then
  # optimization stops.
  return false
end

result_ode = DiffEqFlux.sciml_train(loss, p,
                                    cb = callback,
                                    maxiters = 100)
# result_ode = DiffEqFlux.sciml_train(loss, p, ADAM(0.1), cb = callback)



# =======================================================
# Neural Ordinary Differential Equations with Flux.train!


u0 = Float32[2.; 0.]
datasize = 30
tspan = (0.0f0,1.5f0)

function trueODEfunc(du,u,p,t)
    true_A = [-0.1 2.0; -2.0 -0.1]
    du .= ((u.^3)'true_A)'
end
t = range(tspan[1],tspan[2],length=datasize)
prob = ODEProblem(trueODEfunc,u0,tspan)
ode_data = Array(solve(prob,Tsit5(),saveat=t))

dudt2 = Chain(x -> x.^3,
             Dense(2,50,tanh),
             Dense(50,2))

p, re = Flux.destructure(dudt2) # use this p as the initial condition!
dudt(u,p,t) = re(p)(u) # need to restructure for backprop!
prob = ODEProblem(dudt,u0,tspan)

function predict_n_ode()
  Array(solve(prob,Tsit5(),u0=u0,p=p,saveat=t))
end

function loss_n_ode()
    pred = predict_n_ode()
    loss = sum(abs2,ode_data .- pred)
    loss
end

loss_n_ode() # n_ode.p stores the initial parameters of the neural ODE

cb = function (;doplot=false) # callback function to observe training
  pred = predict_n_ode()
  display(sum(abs2,ode_data .- pred))
  # plot current prediction against data
  pl = scatter(t,ode_data[1,:],label="data")
  scatter!(pl,t,pred[1,:],label="prediction")
  display(plot(pl))
  return false
end

# Display the ODE with the initial parameter values.
cb()

data = Iterators.repeated((), 1000)
Flux.train!(loss_n_ode, Flux.params(u0,p), data, ADAM(0.05), cb = cb)
