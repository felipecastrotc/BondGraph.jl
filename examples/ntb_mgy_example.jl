### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 818e3cee-1df8-48b1-bd64-227265772680
begin	
	import Pkg
	Pkg.add("SymbolicUtils")
end

# ╔═╡ 9469b113-95f2-4e67-b599-ef9bab6865ae
Pkg.add("BondGraph")

# ╔═╡ 25edf63e-84fe-11ec-13a0-e93833c9b99d
begin
	using BondGraph, Plots, Symbolics.Latexify
	import BondGraph: t, D
	
	using DifferentialEquations
end

# ╔═╡ f8656632-5785-4b65-83be-ef0db325323d
begin
	@named L = Mass(m = 0.5);
	@named R = Damper(c = 1.0);
	@named Uₐ = Se(12.0);
	
end

# ╔═╡ c1bcea2c-eb82-4082-9b3b-b9be15ba143c
begin
	@named J = Mass(m = 0.01)
	@named b = Damper(c = 0.1)
	@named Tₗ = Se(1.0)
end

# ╔═╡ 6611c596-4909-4671-87bd-f7df1b76e7a9
g = 0.01

# ╔═╡ 861d9a17-a9ba-4de3-98b1-37c5dd2238f6
begin
	@named je = Junction1(Uₐ, -R, -L, sgn = -1)
	@named jm = Junction1(Tₗ, -b, -J)
	@named gy = mGY(je, jm, g = g)
end

# ╔═╡ b8ca3891-ab03-452f-8a41-60e752d2a0f5
equations(gy)

# ╔═╡ ae9d9340-3705-4779-ac84-d0b9db0fba0f
@named sys = reducedobs(structural_simplify(gy))

# ╔═╡ 15c818c5-2203-47f1-8c78-5bb8601dec0f
equations(sys)

# ╔═╡ 974a0ab1-3800-4de1-a26a-7625899d5e43
begin
	prob = ODEProblem(sys, [], (0.0, 4.0))
	sol = solve(prob)
	plot(sol)
end

# ╔═╡ Cell order:
# ╠═f8656632-5785-4b65-83be-ef0db325323d
# ╟─c1bcea2c-eb82-4082-9b3b-b9be15ba143c
# ╠═6611c596-4909-4671-87bd-f7df1b76e7a9
# ╟─861d9a17-a9ba-4de3-98b1-37c5dd2238f6
# ╟─b8ca3891-ab03-452f-8a41-60e752d2a0f5
# ╠═ae9d9340-3705-4779-ac84-d0b9db0fba0f
# ╟─15c818c5-2203-47f1-8c78-5bb8601dec0f
# ╟─974a0ab1-3800-4de1-a26a-7625899d5e43
# ╟─818e3cee-1df8-48b1-bd64-227265772680
# ╠═25edf63e-84fe-11ec-13a0-e93833c9b99d
# ╠═9469b113-95f2-4e67-b599-ef9bab6865ae
