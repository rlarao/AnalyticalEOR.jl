### A Pluto.jl notebook ###
# v0.19.9

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

# ╔═╡ d8d0af8e-0c5e-11ed-1370-cb8d2ab22ce7
# ╠═╡ show_logs = false
begin
	import Pkg
	Pkg.activate("/home/rick/.julia/dev/AnalyticalEOR/Project.toml")
	using AnalyticalEOR, PlutoUI, Plots, LaTeXStrings
end

# ╔═╡ 7f366edf-af92-4eb1-828f-1869760a1733
begin
	c_labels = [L"Na" L"Mg" L"Ca" L"Cl"]
	a_labels = [L"XNa" L"X_2Mg" L"X_2Ca"]
	c_colors = [:gold :green :blue :pink]
	a_colors = [:gold :green :blue]
end

# ╔═╡ 43b3e07f-09ba-4f91-afe8-0d676d2b1d58
md"""#### Experiment 3.2 Voegelin et al. (2000)"""

# ╔═╡ b198dd49-34ff-4176-9648-8226936ccce6
# begin
# 	ζᵢ = [	1e-9, # Na
# 			1e-9/2, # Mg
# 			0.01/2, # Ca
# 		]
# 	push!(ζᵢ, ζᵢ[1] + 2ζᵢ[2] + 2ζᵢ[3]) # Anion Cl conc. by charge balance
	
# 	# * Injected concentrations
# 	ζⱼ = [	0.095, 	# Na
# 			0.0052/2, # Mg
# 			1e-9/2, 	# Ca
# 		]
		
# 	push!(ζⱼ, ζⱼ[1] + 2ζⱼ[2] + 2ζⱼ[3] ) # Anion Cl conc. by charge balance
	
# end

# ╔═╡ b65fca11-04e6-4042-b7ea-79ae83e623cd
md"""#### Experiment 6.1 Voegelin et al. (2000)"""

# ╔═╡ 8753454d-6a9e-442f-81bf-108803c7d762


# ╔═╡ cbdc5f1f-6e7f-4565-a6a4-1301f508566b
# begin
# 	ζᵢ = [	0.0017, # Na
# 			0.00124/2, # Mg
# 			0.0062/2, # Ca
# 		]
# 	push!(ζᵢ, ζᵢ[1] + 2ζᵢ[2] + 2ζᵢ[3]) # Anion Cl conc. by charge balance
	
# 	# * Injected concentrations
# 	ζⱼ = [	0.47, 	# Na
# 			0.098/2, # Mg
# 			0.023/2, 	# Ca
# 		]
		
# 	push!(ζⱼ, ζⱼ[1] + 2ζⱼ[2] + 2ζⱼ[3] ) # Anion Cl conc. by charge balance
	
# end

# ╔═╡ e6f3b900-efdb-4143-9bc2-57965a5f1a21
# begin
# 	ζᵢ = [	0.47, # Na
# 			0.098/2, # Mg
# 			0.023/2, # Ca
# 		]
# 	push!(ζᵢ, ζᵢ[1] + 2ζᵢ[2] + 2ζᵢ[3]) # Anion Cl conc. by charge balance
	
# 	# * Injected concentrations
# 	ζⱼ = [	0.0017, 	# Na
# 			0.00124/2, # Mg
# 			0.0062/2, 	# Ca
# 		]
		
# 	push!(ζⱼ, ζⱼ[1] + 2ζⱼ[2] + 2ζⱼ[3] ) # Anion Cl conc. by charge balance
	
# end

# ╔═╡ f058fde4-864d-4363-8551-6d6ce508f1cc
md""" ### Experiment 6.2 Voegelin et al. (2000)"""

# ╔═╡ 7fde1209-3683-45db-9a6a-d5a0a8cacc3f
begin
	ζᵢ = [	0.0017, # Na
			0.00124/2, # Mg
			0.0062/2, # Ca
		]
	push!(ζᵢ, ζᵢ[1] + 2ζᵢ[2] + 2ζᵢ[3]) # Anion Cl conc. by charge balance
	
	# * Injected concentrations
	ζⱼ = [	0.47, 	# Na
			0.098/2, # Mg
			0.023/2, 	# Ca
		]
		
	push!(ζⱼ, ζⱼ[1] + 2ζⱼ[2] + 2ζⱼ[3] ) # Anion Cl conc. by charge balance
end

# ╔═╡ 590ad963-6e67-422a-9682-3e12ab2d69b7
md""" ### Equilibrium Constants"""

# ╔═╡ febdeadf-fff3-4b3d-ad18-4be6658d62b7
begin
	ν =  [1, 2, 2, 1] # charges
	
	# * Equilibrium constants and 
	K₂₁ = 10^1.67148
	K₃₁ = 10^1.83148     
	K₂₃ =  10^-0.16
	
	# * Cation Exchange Capacity
	ρ = 2.8082271
	ϕ = 0.59
	cec = 0.06
	Z = cec * ((1 - ϕ) / ϕ) * ρ # Conversion of cation exchange capacity into moles/liter
	
	ec = IonExchangeProblem(K₂₁, K₃₁, K₂₃, Z, ν)
end

# ╔═╡ 3ee83a07-85cb-4996-82bb-75ff654dee7c
Z

# ╔═╡ 99bb428f-6f32-41e2-997a-b9d9d34beb2d
it = solve_Ion_Transport(ζᵢ, ζⱼ, ec)

# ╔═╡ c0184f6b-6c80-4ea3-b6e2-2aceeabeafcc
@bind t Slider(0.01:0.01:5, default=0.1)

# ╔═╡ 0193979e-496c-41f7-94ed-4d64364cb6ee
md""" t = $t"""

# ╔═╡ 61881095-681d-45b3-808e-b93184e24f8b
begin
plot_adsorbed_conc(it, t, a_labels, a_colors)
plot!(yscale=:log10)
end

# ╔═╡ a466c05b-9e88-4153-ab62-342497ce5ed8
begin
plot_flowing_conc(it, t, c_labels, c_colors)
plot!(yscale=:log10)
end

# ╔═╡ Cell order:
# ╠═d8d0af8e-0c5e-11ed-1370-cb8d2ab22ce7
# ╟─7f366edf-af92-4eb1-828f-1869760a1733
# ╟─43b3e07f-09ba-4f91-afe8-0d676d2b1d58
# ╠═b198dd49-34ff-4176-9648-8226936ccce6
# ╟─b65fca11-04e6-4042-b7ea-79ae83e623cd
# ╠═8753454d-6a9e-442f-81bf-108803c7d762
# ╟─cbdc5f1f-6e7f-4565-a6a4-1301f508566b
# ╠═e6f3b900-efdb-4143-9bc2-57965a5f1a21
# ╟─f058fde4-864d-4363-8551-6d6ce508f1cc
# ╠═7fde1209-3683-45db-9a6a-d5a0a8cacc3f
# ╟─590ad963-6e67-422a-9682-3e12ab2d69b7
# ╠═febdeadf-fff3-4b3d-ad18-4be6658d62b7
# ╠═3ee83a07-85cb-4996-82bb-75ff654dee7c
# ╠═99bb428f-6f32-41e2-997a-b9d9d34beb2d
# ╟─0193979e-496c-41f7-94ed-4d64364cb6ee
# ╟─c0184f6b-6c80-4ea3-b6e2-2aceeabeafcc
# ╟─61881095-681d-45b3-808e-b93184e24f8b
# ╟─a466c05b-9e88-4153-ab62-342497ce5ed8
