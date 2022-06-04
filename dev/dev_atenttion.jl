using Symbolics

@variables a b c d t


x = cos(a) + sin(b) + exp(-2/b)

Symbolics.TreeViews(cos(a) + sin(b))

Symbolics.toexpr(x)



X = [a, b,]

W = ones(2, length(X))
W = rand(size(W)...) .> 0.3


X2 = vcat(sin.(W*X),
exp.(W*X),
log.(W*X),
1.0./(W*X))

W2 = ones(2, length(X2))
W2
W2*tmp

tmp = vcat(sin.(W2*X2),
exp.(W2*X2),
log.(W2*X2),
1.0./(W2*X2))


exp(-1*log(2*10))


exp(-3*log(10))

x = -1*log(a)
y = -1*log(b)
z = -1*log(c)

x + y + z = log(a) + log(b) + log(x) = log(a*b*c)

exp(x)*exp(y)*exp(z) = exp(x + y + z) 


exp(0.5*log(10) + 2*log(8) -3*log(2)) 

exp()
10^0.5*8^2/2^3

exp(3*log(complex(-1)))
exp(2.2*log(complex(-1)))

(Complex(-1))^(2.2)



20^-1


exp(2*log(10))

@variables reₛ, dₛ, ϵₛ, ρₛ, Lₛ

re(p)*abs(u[1]), d(p), ϵ(p), ρ(p), L(p)

u = [reₛ, dₛ, ϵₛ, ρₛ, Lₛ]

simplify(γ[1]*exp(sum(γ[2:end].*log.(u))))
simplify(γ[1]*exp(sum(γ[2:end].*log.(u))))


simplify(γ[1]*prod(u.^γ[2:end]))


0.0001^0.00048

2^0.7052

1/0.7502