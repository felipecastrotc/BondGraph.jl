### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 69e2c94d-876d-42b8-98ae-53bab44d72cf
begin
	import Pkg
    Pkg.activate(".")
	
	using DifferentialEquations, DiffEqSensitivity, ModelingToolkit
	using Symbolics, Symbolics.Latexify
	using LinearAlgebra, Statistics
	using Flux
	using Flux.Optimise: ADAM, update!
	using Plots, Printf
	
	include("lib_bg.jl")
end

# ╔═╡ 1da302c0-61eb-4f6b-a7b9-53329063a528
md"""## Bond graph modelling and Neural ODE"""

# ╔═╡ f3ebfb0c-8280-4871-950a-af1c2fe6775e
md"""
Declare the system components:
* Mass;
* Linear spring;
* Cubic spring;
* Damper.
"""

# ╔═╡ c9b5a317-7a74-4c82-bd05-d71adb3e3daf
@named mass = Mass(m=1.);

# ╔═╡ 8148306d-f5ac-4cc2-ba8a-f0c6c03e1802
@named spring = Spring(k=5., x=1.);

# ╔═╡ 1f0b94b4-6634-4881-b618-900cef590028
@named spring3 = Spring3(k=5., x=1.);

# ╔═╡ f5f5742b-cd28-40b6-954d-b8b24008b3a1
@named damper = Damper(c=0.2);

# ╔═╡ d16eda67-f782-42f4-992f-2879ab7dbba6
md"""
### 1DOF Mass-spring-damper system
#### Linear spring

Modeling using the bondgraph technique we obtain the following set of equations:
"""

# ╔═╡ 777a2ca6-dbf1-4318-973e-9a60ffd2ab52
@named model_1dof = Junction1(spring, mass, damper, couple=false)

# ╔═╡ 33d214c1-5f1b-4e8c-adaa-e3e11253b1db
md""" Simplifying the system"""

# ╔═╡ 4d24d91a-ac54-4f48-98cf-266e18c974dc
sys_1dof = structural_simplify(model_1dof)

# ╔═╡ fafc29b1-8612-49f4-8366-1955ed91ee8a
latexify(sys_1dof) * "sdas"

# ╔═╡ cfccc3f0-b2d8-46cc-9b77-ffa804551e0d
md""" Simulate the system:

Time: $(@bind t_sim_sys_1dof html"<input type=range min=1 max=50>") """

# ╔═╡ dc437255-5f00-406b-a986-4a65253c9428
(t_sim_sys_1dof,)

# ╔═╡ b3f287b2-d562-498f-9a06-8047a0279f5a
begin
	prob_1dof = ODEProblem(sys_1dof, [], (0., t_sim_sys_1dof))
	sol_1dof = solve(prob_1dof)
	
	plot(sol_1dof)
end

# ╔═╡ 107db558-1f28-4291-bec6-5d8029e98f45
md"""
#### Cubic spring

Modeling using the bondgraph technique we obtain the following set of equations:
"""

# ╔═╡ 0ce9d827-abe0-4c0d-996c-084796646dd3
@named model_1dof3 = Junction1(spring3, mass, damper, couple=false)

# ╔═╡ 276d732b-76e4-4f8d-84b3-1fb288494a3b
md""" Simplifying the system"""

# ╔═╡ 17ccb290-f3b6-4bb5-8c2b-3892703a3acd
sys_1dof3 = structural_simplify(model_1dof3)

# ╔═╡ 8d99d4df-ece1-4e24-b085-0054632ef485
md""" Simulate the system:

Time: $(@bind t_sim_sys_1dof3 html"<input type=range min=1 max=50>") """

# ╔═╡ 959f666e-4263-40bc-bcd7-9d064fd3ed63
(t_sim_sys_1dof,)

# ╔═╡ 05ba0c86-867c-4869-bc8a-8cf47551581f
begin
	prob_1dof3 = ODEProblem(sys_1dof3, [], (0., t_sim_sys_1dof3))
	sol_1dof3 = solve(prob_1dof3)
	
	plot(sol_1dof3)
end

# ╔═╡ ed83f3b3-b5e3-40a0-aff1-2b68498abd45
md"""
### 2DOF Mass-spring-damper system
#### Linear spring

Modeling using the bondgraph technique we obtain the following set of equations:
"""

# ╔═╡ be6368c3-868e-48db-859c-762c6640dea8
begin
	@named sd = Junction1(spring, damper)
	
	M2 = [Junction1(spring, damper, mass; name = :b1)];
	push!(M2, Junction0(sd; name = :b2j0, subsys = M2[1]));
	push!(M2, Junction1(mass; name = :b2j1, subsys = M2[2], couple = false));	
	
	@named model_2dof = ODESystem(equations(M2), t)
end;

# ╔═╡ 2c3d02f3-0d93-4ea6-8519-d27ad64a0e07
md""" Simplifying the system"""

# ╔═╡ 8bc1b7c3-82d3-4149-9ded-9c650faa07ba
sys_2dof = structural_simplify(model_2dof)

# ╔═╡ 02f0e9c0-ffce-46f3-9a57-7c99b907f6b4
latexify(sys_2dof)*"1"

# ╔═╡ 61dd3028-543f-40b8-98aa-79f5bc9e7020
md""" Simulate the system:

Time: $(@bind t_sim_sys_2dof html"<input type=range min=1 max=50>") """

# ╔═╡ 7d8f3b64-10a4-45fa-849b-aee554965e47
(t_sim_sys_2dof,)

# ╔═╡ 5f100395-5c1e-45d1-bc3a-d2cb0296652c
begin
	
	prob_2dof = ODEProblem(sys_2dof, [], (0., t_sim_sys_2dof))
	sol_2dof = solve(prob_2dof)
	
	plot(sol_2dof)
	
end

# ╔═╡ 0262d087-75f2-4f02-99d9-978603b7ac37
md"""
### 4DOF Mass-spring-damper system
#### Linear spring

Modeling using the bondgraph technique we obtain the following set of equations:
"""

# ╔═╡ 1651d15b-fcdd-4102-9f8b-bd2f5df4c8b9
begin	
	M4 = [Junction0(sd, mass; name=:b3j0)];
	push!(M4, Junction1(mass; name=:b2j1, subsys=M4[1]));
	push!(M4, Junction0(sd; name=:b2j0, subsys=M4[2]));
	push!(M4, Junction1(mass; name=:b1j1, subsys=M4[3]));
	push!(M4, Junction0(sd; name=:b1j0, subsys=M4[4]));
	push!(M4, Junction1(mass, spring, damper; name=:b0, subsys=M4[5], couple=false));
	
	@named model_4dof = ODESystem(equations(M4), t);
end;

# ╔═╡ 4501ede9-5e01-44ee-a5cd-2f610b1364d6
md"""Simplifying"""

# ╔═╡ c69023cd-6052-4555-9ba9-bcfd90125e8e
sys_4dof = structural_simplify(model_4dof)

# ╔═╡ 27435a45-2924-4428-a27b-7827f551b601
md""" Simulate the system:

Time: $(@bind t_sim_sys_4dof html"<input type=range min=1 max=50>") """

# ╔═╡ 5a476d63-07e6-40dd-8067-e11aa422938e
(t_sim_sys_4dof,)


# ╔═╡ 2c93e3aa-6254-49f2-bc72-2dff4965313b
begin
	
prob_4dof = ODEProblem(sys_4dof, [], (0., t_sim_sys_4dof))
sol_4dof = solve(prob_4dof)

plot(sol_4dof)

	
end

# ╔═╡ d51e75dc-77a5-49db-a1bd-ac9bbe176571
md"""
### NeuralODE + Bondgraph
#### 1DOF

The system is modelled using the bondgraph technique but, for instance, unknown spring function is approximated using a neural network or any function.

Firstly, it is required to create a new spring model that has unknown function, $h(x, p)$.

"""

# ╔═╡ 416cd952-bba7-4b04-bb3e-1cd13cccf560
@register h(x, p)

# ╔═╡ c6f0c12f-96fd-4d28-8b12-09ba31cb7cce
md"""Create new spring model, SpringU:"""

# ╔═╡ a2bbe472-cde7-44fc-93fa-5eaf42bd2c19
md"""Declare the new spring:"""

# ╔═╡ 9137c9fb-906f-4e6f-8470-8a01096ee60b
md""" Guess a function type for $h(x,p)$."""

# ╔═╡ 7a4ca281-91e0-435e-865e-7555a44efa71
h(x, p) = x*p

# ╔═╡ 0f7a2c98-fb1a-444c-9798-1276a47b0324
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

# ╔═╡ 94b11ce2-ff6f-4fea-8a50-675117b0fbf5
@named springU = SpringU(x=1.)

# ╔═╡ e6763942-e37a-47b8-bb80-e3100f58febd
md"""Now, model the mass-spring-damper system using the generic SpringU, instead of the linear spring."""

# ╔═╡ c6df6aef-7faa-43d0-a8e0-5125225616c1
@named model_1dofu = Junction1(mass, springU, damper, couple=false)

# ╔═╡ f3c7935f-0ee2-42af-ac72-1e2f9afb4642
md"""Simplifying the system:"""

# ╔═╡ 306e69cb-d997-4613-ac8d-f77b0c383307
sys_1dofu = structural_simplify(model_1dofu)

# ╔═╡ 61c7bbfc-f6c7-4897-802b-229f96b680d8
latexify(sys_1dofu)*"1"

# ╔═╡ b9585204-8336-407b-96dc-e6a5647a7041
md""" Simulate the system:

Time: $(@bind t_sim_sys_1dofu html"<input type=range min=1 max=50>")
p: $(@bind p_1dofu html"<input type=range min=0 max=20 step=0.1>") 
"""

# ╔═╡ 7b1b5a63-a39a-4de7-8a3b-2cb7110d7f16
["t" => t_sim_sys_1dofu, "p" => p_1dofu]

# ╔═╡ 45b5ec93-47b9-448d-9530-9ba283b7a6e9
begin
	prob_1dofu = ODEProblem(sys_1dofu, [], (0., t_sim_sys_1dofu), [springU.p => p_1dofu])
	sol_1dofu = solve(prob_1dofu)
	plot(sol_1dofu)
end

# ╔═╡ 4343e21d-1233-4336-8664-05b6226826c6
md""" Now, suppose you know the truth spring stiffness is $5$, but you do not know it previously and you try to find a $p$ that matches it:"""

# ╔═╡ c75d18b4-d56f-49c5-ba9d-c023786571aa
begin
	# Declare a function with p
	function pred_ode(p, array=true)
		# ps = Array(prob_1dofu.p)
		ps = [1.0, p, 0.2];
		# Set the initial conditions
		u0 = prob_1dofu.u0
	    sol = solve(prob_1dofu,Tsit5(),u0=u0,p=ps,saveat=t_rng)
	    if array
	        Array(sol)
	    else
	        sol
	    end
	end
	
	function train(p, y, it=300)

        # The x was changed due the position of the variables
        loss(x, y) = mean((circshift(x, 1) .- y).^2);
		# loss(x, y) = mean((x .- y).^2);
        
        # Train
        a = Animation()
        hist = zeros(it);
        for i in 1:it
            g = gradient((p) -> loss(pred_ode(p), y), p);
            p += -0.5*g[1];     # Simple gradient descent
            hist[i] = loss(pred_ode(p), y);
            @printf("It: %d - loss: %.3e \n", i, hist[i]);
            
			if hist[i] < 1e-5
				break
			end
			
            # Animation
            plt = plot(pred_ode(p, false))
            plt = plot(plt, sol_1dof(t_rng), xlabel="Time (s)", ylabel="Amplitude")
            frame(a, plt)
        end

        return p, a
    end
	
	t_rng = 0:0.1:10
	y = Array(sol_1dof(t_rng))
	
	p, a = train(p_1dofu, y, 600);
	
	gif(a)
end

# ╔═╡ 30c83012-cc69-44d8-a451-dfe0e857df34
md"""
#### 4DOF

Extending for 4DOF, the system's equation modelled is shown below:

"""

# ╔═╡ 631ec632-2f26-4512-8438-56c5ef138e20
begin
	
	@named sdU = Junction1(springU, damper)
	
	MU = [Junction0(sdU, mass; name=:b3j0)];
	push!(MU, Junction1(mass; name=:b2j1, subsys=MU[1]));
	push!(MU, Junction0(sdU; name=:b2j0, subsys=MU[2]));
	push!(MU, Junction1(mass; name=:b1j1, subsys=MU[3]));
	push!(MU, Junction0(sdU; name=:b1j0, subsys=MU[4]));
	push!(MU, Junction1(mass, springU, damper; name=:b0, subsys=MU[5], couple=false));
	
	@named model_4dofu = ODESystem(equations(MU), t)
	
	sys_4dofu = structural_simplify(model_4dofu)
end

# ╔═╡ 6dcc9425-1ce2-47ed-985c-2a1f73ab344f
md""" Simulate the system:

Time: $(@bind t_sim_sys_4dofu html"<input type=range min=1 max=50>")
"""

# ╔═╡ 0a97dc09-a651-4ea5-96a0-781ebdabcbdc
p_4dofu = [2., 8., 10, 7.]

# ╔═╡ e8be615b-f1bc-4064-a8c9-715a7f095741
begin
	prob_4dofu = ODEProblem(sys_4dofu, [], (0., t_sim_sys_4dofu))
	sol_4dofu = solve(prob_4dofu)
	plot(sol_4dofu)
end

# ╔═╡ 0deae790-cf20-42e8-a2e0-160b445bf4e7
md""" Now trying to the spring stiffness:"""

# ╔═╡ 07c559fe-c877-4363-a139-a6fa1570e6d2
begin
	# Declare a function with p
	function pred_ode4(p, array=true)
		# ps = Array(prob_1dofu.p)
		ps = [p[1], 0.2, 1.0, 1.0, p[2], 0.2, 1.0, p[3], 0.2, 1.0, p[4], 0.2]
		# Set the initial conditions
		u0 = prob_4dofu.u0
	    sol = solve(prob_4dofu,Tsit5(),u0=u0,p=ps,saveat=t_rng)
	    if array
	        Array(sol)
	    else
	        sol
	    end
	end
	
	function train4(p, y, it=300)

        # The x was changed due the position of the variables
		loss(x, y) = mean((x .- y).^2);
        
        # Train
        a = Animation()
        hist = zeros(it);
        for i in 1:it
            g = gradient((p) -> loss(pred_ode4(p), y), p);
            p += -2*g[1];     # Simple gradient descent
            hist[i] = loss(pred_ode4(p), y);
            @printf("It: %d - loss: %.3e \n", i, hist[i]);
            
			if hist[i] < 1e-5
				break
			end
			
            # Animation
            plt = plot(pred_ode4(p, false)[3:4, :]')
            plt = plot(plt, sol_4dof(t_rng)[3:4, :]', xlabel="Time (s)", ylabel="Amplitude ()")
            frame(a, plt)
        end

        return p, a
    end
	
	y4 = Array(sol_4dof(t_rng))
	
	p4, a4 = train4(p_4dofu, y4, 600);
	
	gif(a4)
end

# ╔═╡ f4a4d755-20ef-4d69-b7fb-f82be2ba5fb6


# ╔═╡ 1a0442d6-18d5-4923-83ca-6418c0a5a055
md"""### Import functions and libraries"""

# ╔═╡ Cell order:
# ╟─1da302c0-61eb-4f6b-a7b9-53329063a528
# ╟─f3ebfb0c-8280-4871-950a-af1c2fe6775e
# ╠═c9b5a317-7a74-4c82-bd05-d71adb3e3daf
# ╠═8148306d-f5ac-4cc2-ba8a-f0c6c03e1802
# ╠═1f0b94b4-6634-4881-b618-900cef590028
# ╠═f5f5742b-cd28-40b6-954d-b8b24008b3a1
# ╟─d16eda67-f782-42f4-992f-2879ab7dbba6
# ╠═777a2ca6-dbf1-4318-973e-9a60ffd2ab52
# ╟─33d214c1-5f1b-4e8c-adaa-e3e11253b1db
# ╟─4d24d91a-ac54-4f48-98cf-266e18c974dc
# ╠═fafc29b1-8612-49f4-8366-1955ed91ee8a
# ╟─cfccc3f0-b2d8-46cc-9b77-ffa804551e0d
# ╟─dc437255-5f00-406b-a986-4a65253c9428
# ╠═b3f287b2-d562-498f-9a06-8047a0279f5a
# ╟─107db558-1f28-4291-bec6-5d8029e98f45
# ╠═0ce9d827-abe0-4c0d-996c-084796646dd3
# ╟─276d732b-76e4-4f8d-84b3-1fb288494a3b
# ╟─17ccb290-f3b6-4bb5-8c2b-3892703a3acd
# ╟─8d99d4df-ece1-4e24-b085-0054632ef485
# ╟─959f666e-4263-40bc-bcd7-9d064fd3ed63
# ╟─05ba0c86-867c-4869-bc8a-8cf47551581f
# ╟─ed83f3b3-b5e3-40a0-aff1-2b68498abd45
# ╠═be6368c3-868e-48db-859c-762c6640dea8
# ╟─2c3d02f3-0d93-4ea6-8519-d27ad64a0e07
# ╟─8bc1b7c3-82d3-4149-9ded-9c650faa07ba
# ╠═02f0e9c0-ffce-46f3-9a57-7c99b907f6b4
# ╟─61dd3028-543f-40b8-98aa-79f5bc9e7020
# ╟─7d8f3b64-10a4-45fa-849b-aee554965e47
# ╟─5f100395-5c1e-45d1-bc3a-d2cb0296652c
# ╟─0262d087-75f2-4f02-99d9-978603b7ac37
# ╠═1651d15b-fcdd-4102-9f8b-bd2f5df4c8b9
# ╟─4501ede9-5e01-44ee-a5cd-2f610b1364d6
# ╠═c69023cd-6052-4555-9ba9-bcfd90125e8e
# ╟─27435a45-2924-4428-a27b-7827f551b601
# ╟─5a476d63-07e6-40dd-8067-e11aa422938e
# ╟─2c93e3aa-6254-49f2-bc72-2dff4965313b
# ╟─d51e75dc-77a5-49db-a1bd-ac9bbe176571
# ╠═416cd952-bba7-4b04-bb3e-1cd13cccf560
# ╟─c6f0c12f-96fd-4d28-8b12-09ba31cb7cce
# ╠═0f7a2c98-fb1a-444c-9798-1276a47b0324
# ╟─a2bbe472-cde7-44fc-93fa-5eaf42bd2c19
# ╠═94b11ce2-ff6f-4fea-8a50-675117b0fbf5
# ╟─9137c9fb-906f-4e6f-8470-8a01096ee60b
# ╠═7a4ca281-91e0-435e-865e-7555a44efa71
# ╟─e6763942-e37a-47b8-bb80-e3100f58febd
# ╠═c6df6aef-7faa-43d0-a8e0-5125225616c1
# ╟─f3c7935f-0ee2-42af-ac72-1e2f9afb4642
# ╟─306e69cb-d997-4613-ac8d-f77b0c383307
# ╠═61c7bbfc-f6c7-4897-802b-229f96b680d8
# ╟─b9585204-8336-407b-96dc-e6a5647a7041
# ╟─7b1b5a63-a39a-4de7-8a3b-2cb7110d7f16
# ╟─45b5ec93-47b9-448d-9530-9ba283b7a6e9
# ╟─4343e21d-1233-4336-8664-05b6226826c6
# ╠═c75d18b4-d56f-49c5-ba9d-c023786571aa
# ╟─30c83012-cc69-44d8-a451-dfe0e857df34
# ╟─631ec632-2f26-4512-8438-56c5ef138e20
# ╟─6dcc9425-1ce2-47ed-985c-2a1f73ab344f
# ╟─0a97dc09-a651-4ea5-96a0-781ebdabcbdc
# ╟─e8be615b-f1bc-4064-a8c9-715a7f095741
# ╟─0deae790-cf20-42e8-a2e0-160b445bf4e7
# ╟─07c559fe-c877-4363-a139-a6fa1570e6d2
# ╠═f4a4d755-20ef-4d69-b7fb-f82be2ba5fb6
# ╠═1a0442d6-18d5-4923-83ca-6418c0a5a055
# ╠═69e2c94d-876d-42b8-98ae-53bab44d72cf
