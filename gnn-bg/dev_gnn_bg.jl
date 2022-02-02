using CSV
using LinearAlgebra
using Plots
using LaTeXStrings
using ScikitLearn
using Statistics
using DataFrames
using Serialization
using Random
using Symbolics
using LsqFit



# [f, e] * Aᵣ = [fᵣ, eᵣ]

k = 2;
m = 2.5;
b = 0.2;
A = [0 1; -k/m -b/m];

#   [ẋ    x]
X₀ = [0, 0.1];

n = 1000;
X = zeros(2, n);
Δt = 0.01;
for i in 1:n
    X[:, i] = X₀;
    Ẋ = A*X₀;
    X₀ += Ẋ*Δt;
end

t = 0:Δt:(n-1)*Δt;
plot(t, X[2, :])
plot!(t, X[1, :])


q = X[2, :];
f = X[1, :].*k/m;
e = k.*q;
p = m.*f;

Eₖ = 0.5.*q.*e;
Eₚ = 0.5.*f.*p;
plot(t, Eₖ)
plot!(t, Eₚ)
plot!(t, cumsum((b.*f).*f.*Δt))


plot(t, Eₖ .+ Eₚ)
plot!(t, Eₖ)
plot!(t, Eₚ)
plot!(t, cumsum((f.^2)*b*Δt))


plot(t, cumsum(f.^2*b*Δt) + (Eₖ .+ Eₚ))
plot(t, cumsum(f.^2*b*Δt))
plot!(t, (Eₖ[1] .+ Eₚ[1]) .- (Eₖ .+ Eₚ))



using DifferentialEquations


function msd!(du, u, p, t)
    A = [0 1; -k/m -b/m];
    A₃ = [0 0; -k₃/m 0];
    tmp = A*u + A₃*(u.^3);
    du[1] = tmp[1];
    du[2] = tmp[2];
end

k = 10;
k₃ = 3;
k₃ = 0;
m = 0.5;
b = 0;
b = 0.1;

u = [0, 1];
# u = [1, 0.2];
tₛ = (0.0, 100.0);

du = [0.0, 0.0]
prob = ODEProblem(msd!, u, tₛ);
sol = solve(prob);

Δt = 0.01
# t = 0:Δt:1;
t = 0:Δt:20;
t = 0:Δt:4;
X = sol(t);
plot(t, X')

q = X[2, :];
# f = X[1, :].*k/m;
f =  vcat([0], -diff(X[2, :])/Δt);

e = k₃.*(q.^3);
# e = k.*q;
p = m.*f;

Eₖ = 0.5.*q.*e;
Eₚ = 0.5.*f.*p;
plot(t, Eₖ)
plot!(t, Eₚ)
plot!(t, cumsum((b.*f).*f.*Δt))

plot(t, cumsum((b.*f).*f.*Δt))
plot!(t, (Eₖ + Eₚ)[1] .- (Eₖ + Eₚ))

plot(cumsum((b.*f).*f.*Δt) .+ (Eₖ + Eₚ))


# Analytical

x(t) = exp.(-0.1.*t).*(0.1.*cos.(1.0*t) .+ (0.01).*sin.(1.0.*t))
x(t) = exp.(-0.0500.*t).*(0.1.*cos.(.705*t) .+ (0.00709).*sin.(0.705.*t))

xx = x(t);
xv = vcat([0], -diff(xx)/Δt);

plot(t, xx)
plot!(t, X[2, :])

plot(t, xv)
plot!(t, X[1, :])


e = k.*xx;
p = m.*xv;
q = xx;
f = xv;

Eₖ = 0.5.*q.*e;
Eₚ = 0.5.*f.*p;

Eₖ = (k*q.^2)./2;
Eₚ = (m*f.^2)./2;

plot(t, Eₖ)
plot!(t, Eₚ)
plot(t, Eₖ .+ Eₚ)

(Eₖ .+ Eₚ)



# NN
using Flux
using Flux.Optimise: ADAM, update!
using Flux.Zygote
using Printf


module TT
    struct Compliance
        e::Float64
        f::Float64
        q::Float64
        p::Float64
    end
    struct Resistance
        e::Float64
        f::Float64
        q::Float64
        p::Float64
    end
    struct Inertance
        e::Float64
        f::Float64
        q::Float64
        p::Float64
    end
    mutable struct element{T, Y}
        port::T
        model::Y
        Δt::Float64
    end
end


function (C::TempTypes.Compliance)(f, q)
    out = C.m([C.e, f, q, C.p])
    C.e += out[1]*C.Δt
    C.f += out[2]*C.Δt
    C.q += out[3]*C.Δt
    C.p += out[4]*C.Δt
end

function (R::TempTypes.Resistance)(f, q)
    out = R.m([R.e, f, q, R.p])
    R.e += out[1]*R.Δt
    R.f += out[2]*R.Δt
    R.q += out[3]*R.Δt
    R.p += out[4]*R.Δt
end

function (I::TempTypes.Inertance)(f, q)
    out = I.m([I.e, f, q, I.p])
    I.e += out[1]*I.Δt
    I.f += out[2]*I.Δt
    I.q += out[3]*I.Δt
    I.p += out[4]*I.Δt
end


# Trial 1

mc = Chain(Dense(4, 16, relu, bias=true), Dense(16, 4, relu, bias=true));
mr = Chain(Dense(4, 16, relu, bias=true), Dense(16, 4, relu, bias=true));
mi = Chain(Dense(4, 16, relu, bias=true), Dense(16, 4, relu, bias=true));

# c = TT.element(TT.Compliance(0., 0., 0.1, 0.), mc, Δt)
# r = TT.element(TT.Resistance(0., 0., 0., 0.), mr, Δt)
# i = TT.element(TT.Inertance(0., 0., 0., 0.), mi, Δt)


function loss(dc, dr, di, d_loss)
    lc = zeros(4);
    lr = zeros(4);
    li = zeros(4);
    ctr_e = 0.0; 
    ctr_f = 0.0; 
    ctr_p = 0.0; 

    n = (size(dc, 1) - 1);
    for i in 1:n
        i = 100
        oc = mc(dc[i, :]);
        or = mr(dr[i, :]);
        oi = mi(di[i, :]);

        # Constraints
        λ = [0.075, 0.05, 0.05, 0.075];
        ctr_e += λ[1]*(oc[1] + or[1] + oi[1])^2;
        ctr_f += λ[2]*(oc[2] - or[2])^2  + λ[3]*(or[2] - oi[2])^2;
        ctr_p += λ[4]*(0.5*oc[3]*oc[1] + 0.5.*oi[4].*oi[2] - d_loss[i+1])^2;
        # Loss
        lc += 0.25*(oc - dc[i+1, :]).^2;
        lr += 0.25*(or - dr[i+1, :]).^2;
        li += 0.25*(oi - di[i+1, :]).^2;
    end

    return norm(sqrt.(lc./n) .+ sqrt.(lr./n) .+ sqrt.(li./n)) + ctr_e/n + ctr_f/n + ctr_p/n;
    # return norm(sqrt.(lc./n) .+ sqrt.(lr./n) .+ sqrt.(li./n) .+ ctr_e/n .+ ctr_f/n .+ ctr_p/n);

end

mc = Chain(Dense(4, 32, relu, bias=true), Dense(32, 4, relu, bias=true));
mr = Chain(Dense(4, 32, relu, bias=true), Dense(32, 4, relu, bias=true));
mi = Chain(Dense(4, 32, relu, bias=true), Dense(32, 4, relu, bias=true));

dc = hcat(k*q, f, q, p)
dr = hcat(b*f, f, q, p)
di = hcat(-(k*q + b*f), f, q, p)
d_loss = cumsum((b.*f).*f.*Δt);

ps = Flux.params(mc, mr, mi);

# Initialise optimiser
opt = ADAM(0.001);

# opt = RMSProp();


# Iterate
# it = 19;
it = 400;
e = zeros(it);
for i in 1:it
    ∇loss = gradient(() -> loss(dc, dr, di, d_loss), ps);
    update!(opt, ps, ∇loss);
    e[i] = loss(dc, dr, di, d_loss);
    @printf("It: %d - loss: %.3e \n", i, e[i]);
end

plot(e)

mc(dc[301, :])
dc[302, :]


# Trial 2

function loss(dc, dr, di, d_loss)
    lc = zeros(4);
    lr = zeros(4);
    li = zeros(4);
    ctr_e = 0.0; 
    ctr_f = 0.0; 
    ctr_p = 0.0; 

    n = (size(dc, 2) - 1);
    for i in 1:n
        oc = nn_c(dc[1:end-1,i]);
        oi = nn_i(di[1:end-1,i]);
        or = nn_i(dr[1:end-1,i]);

        # Constraints
        λ = [0.075, 0.05, 0.05, 0.075];
        ctr_e += λ[1]*(oc[1] + or[1] + oi[1])^2;
        ctr_f += λ[2]*(oc[2] - or[2])^2  + λ[3]*(or[2] - oi[2])^2;
        # ctr_p += λ[4]*(0.5*oc[3]*oc[1] + 0.5.*oi[4].*oi[2] - d_loss[i+1])^2;
        # Loss
        lc += (oc - dc[:, i+1]).^2;
        li += (oi - di[:, i+1]).^2;
        lr += (or - dr[:, i+1]).^2;
    end

    return mean(sqrt.(lc./n) .+ sqrt.(lr./n) .+ sqrt.(li./n)) + ctr_e/n + ctr_f/n;
    # return norm(sqrt.(lc./n) .+ sqrt.(lr./n) .+ sqrt.(li./n) .+ ctr_e/n .+ ctr_f/n .+ ctr_p/n);

end


function loss(dc, di, dr)
    
    out_c = nn_c(dc[1:end-1,1:end-1]);
    out_i = nn_i(di[1:end-1,1:end-1]);
    out_r = nn_r(dr[1:end-1,1:end-1]);
    
    # λ = [0.075, 0.05, 0.05, 0.075];
    λ = [0.1, 0.1, 0.1, 0.1];
    # λ = [1, 1, 1, 1];
    ctr_e = λ[1]*sqrt(mean((out_c[1, :] .+ out_r[1, :] .+ out_i[1, :]).^2));
    ctr_f = λ[2]*sqrt(mean((out_c[2, :] .- out_r[2]).^2))  + λ[3]*sqrt(mean((out_r[2, :] .- out_i[2, :]).^2));

    lc = sqrt(mean((nn_c(dc[1:end-1,1:end-1]) .- dc[:,2:end]).^2));
    li = sqrt(mean((nn_i(di[1:end-1,1:end-1]) .- di[:,2:end]).^2));
    lr = sqrt(mean((nn_i(dr[1:end-1,1:end-1]) .- dr[:,2:end]).^2));
    return li + lc + lr + ctr_e + ctr_f
end 

function loss_sep(dc, di, dr)
    out_c = nn_c(dc[1:end-1,1:end-1]);
    out_i = nn_i(di[1:end-1,1:end-1]);
    out_r = nn_r(dr[1:end-1,1:end-1]);
    
    lc = sqrt(mean((out_c .- dc[:,2:end]).^2));
    li = sqrt(mean((out_i .- di[:,2:end]).^2));
    lr = sqrt(mean((out_r .- dr[:,2:end]).^2));
    return li + lc + lr
end 

nn = 32
σ = leakyrelu;
nn_c = Chain(Dense(3, nn, σ, bias=true), Dense(nn, 4, identity, bias=true));
nn_i = Chain(Dense(3, nn, σ, bias=true), Dense(nn, 4, identity, bias=true));
nn_r = Chain(Dense(3, nn, σ, bias=true), Dense(nn, 4, identity, bias=true));

dc = hcat(e, f, q, p)';
di = hcat(-(e + b*f), f, q, p)';
dr = hcat(b*f, f, q, p)';
d_loss = cumsum((b.*f).*f.*Δt);

ps = Flux.params(nn_c, nn_i, nn_r);

# Initialise optimiser
opt = AMSGrad();
# opt = ADAM(0.001);

# Iterate
# it = 19;
it = 3000;
h = zeros(it);
for i in 1:it
    ∇loss = gradient(() -> loss_sep(dc, di, dr), ps);
    update!(opt, ps, ∇loss);
    h[i] = loss_sep(dc, di, dr);
    @printf("It: %d - loss: %.3e \n", i, h[i]);
end

plot(h)


# Initialise optimiser
opt = AMSGrad();
# opt = ADAM(0.001);
# Iterate
# it = 19;
it = 1000;
h = zeros(it);
for i in 1:it
    ∇losss = gradient(() -> loss(dc, di, dr), ps);
    update!(opt, ps, ∇losss);
    h[i] = loss(dc, di, dr);
    @printf("It: %d - loss: %.3e \n", i, h[i]);
end

plot(h)


i = 1
oc = nn_c(dc[1:end-1,i]);
oi = nn_i(di[1:end-1,i]);
or = nn_r(dr[1:end-1,i]);

oc[1] + oi[1] + or[1]
mean([oc[2], oi[2], or[2]])
dr[2, 2]
# mean([oc[3], oi[3], or[3]])

# Test

tₛ = 20
tmp = zeros(4, Int(tₛ/Δt));
tmp[:, 1] = hcat(e, f, q, p)'[:, 1];
for i in 2:size(tmp, 2)
    tmp[:, i] = nn_c(tmp[1:end-1, i-1]);
end

tt = 0:Δt:(tₛ-Δt)
plot(tt, tmp[3, :])
plot!(t, q)

plot(tt, tmp[1, :])
plot!(t, e)

plot(tt, tmp[2, :])
plot!(t, f)

plot(tt, tmp[4, :])
plot!(t, p)


# Test with only position

function loss(q)
    lc = sqrt(mean((nn_q(q[1:end-1]') .- q[2:end]').^2));
    return lc
end 


nn = 8
σ = relu
nn_q = Chain(Dense(1, nn, σ, bias=true), Dense(nn, 1, identity, bias=true));

ps = Flux.params(nn_q);

# Initialise optimiser
opt = AMSGrad();
# opt = ADAM(0.001);

# Iterate
# it = 19;
it = 1000;
h = zeros(it);
for i in 1:it
    ∇loss = gradient(() -> loss(q), ps);
    update!(opt, ps, ∇loss);
    h[i] = loss(q);
    @printf("It: %d - loss: %.3e \n", i, h[i]);
end


tmp = zeros(1, 200);
tmp[1, 1] = q[1];
for i in 2:size(tmp, 2)
    tmp[:, i] = nn_q(tmp[:, i-1]);
end

plot(tmp')
plot!(q)


# New model

function msd!(du, u, p, t)
    A = [0 1; -k/m -b/m];
    A₃ = [0 0; -k₃/m 0];
    tmp = A*u + A₃*(u.^3);
    du[1] = tmp[1];
    du[2] = tmp[2];
end

k = 0;
k₃ = 0;
k₃ = 6;
m = 0.5;
b = 0;
b = 0.02;

u = [0, 1];
# u = [1, 0.2];
tₛ = (0.0, 100.0);

du = [0.0, 0.0]
prob = ODEProblem(msd!, u, tₛ);
sol = solve(prob);

Δt = 0.01
# t = 0:Δt:1;
t = 0:Δt:4;
t = 0:Δt:20;
X = sol(t);
plot(X[1, :], X[2, :])
plot(t, X')


q = X[2, :];
# f = X[1, :].*k/m;
f =  vcat([0], -diff(X[2, :])/Δt);

e = k₃.*(q.^3);
# e = k.*q;
p = m.*f;

Eₖ = 0.5.*q.*e;
Eₚ = 0.5.*f.*p;
plot(t, Eₖ)
plot!(t, Eₚ)
plot!(t, cumsum((b.*f).*f.*Δt))

plot(t, cumsum((b.*f).*f.*Δt))



