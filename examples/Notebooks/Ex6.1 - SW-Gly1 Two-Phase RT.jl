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
	using AnalyticalEOR, PlutoUI, Plots, Roots, LaTeXStrings
end

# ╔═╡ 79399bed-95e7-43ea-bef4-dd56a3364e0c
begin
	c_labels = ["COO" "SO4" "Gly" "Na"]
	c_colors = [:brown :darkorchid :orangered :gold]
	a_labels = ["XCOO" "X2SO4" "XGly"]
	a_colors = [:brown :darkorchid :orangered]
end

# ╔═╡ febdeadf-fff3-4b3d-ad18-4be6658d62b7
begin
	ν = [1, 2, 1, 1] # charges

    ζᵢ = [9e-3,   # COO	
        7e-3, # S
        1e-9,  # Gly
    ]
	
	# Cation Na conc. by charge balance
    push!(ζᵢ, ν[1] * ζᵢ[1] + ν[2] * ζᵢ[2] + ν[3] * ζᵢ[3])

    # * Injected concentrations
    ζⱼ = [9e-3,   # COO
        1.961e-2,  # S
        0.1/5, # Gly
    ]
	
	# Cation Na conc. by charge balance
    push!(ζⱼ, ν[1] * ζⱼ[1] + ν[2] * ζⱼ[2] + ν[3] * ζⱼ[3])
    
    # * Cation Exchange Capacity
    Z = 0.07

    # * Equilibrium constants and 
    # K₂₁ = 10^(2.478 - log10(Z))
    # K₂₁ = 10^(2.8 - log10(2))
    # K₃₁ = 10^1.78
    K₂₁ = 10^(2.8 - log10(2Z))
    K₃₁ = 10^1.78
    K₂₃ = K₂₁ / K₃₁

	ec = IonExchangeProblem(K₂₁, K₃₁, K₂₃, Z, ν)

end

# ╔═╡ a72e76c9-16ba-4ecd-844d-77460f6333b0
# ╠═╡ show_logs = false
it = solve_Ion_Transport(ζᵢ, ζⱼ, ec)

# ╔═╡ c0184f6b-6c80-4ea3-b6e2-2aceeabeafcc
@bind t Slider(0.1:0.01:3, default=0.5)

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
	# plot!(yscale=:log10)
end

# ╔═╡ af1a7f63-9a1c-4953-9f73-56a0d357669e
md""" #### Two-Phase Flow Solution"""

# ╔═╡ 9eada2f0-d022-41f6-aa4e-b02691a7e251
it.λ

# ╔═╡ 37640243-9cfe-46b5-ab5a-1b9f92dd9d29
begin
	λᵢ = it.λ[27]
	λⱼ = it.λ[1]
	λₘ₁ = it.λ[23]
	λₘ₂ = it.λ[19]
end

# ╔═╡ b83d4f0b-f104-4cbd-8a50-17bb7ce51813


# ╔═╡ 48da0d27-b166-4765-8566-a3d6e22223cc
begin
	ĉᵢ = it.ĉ[27,:]
	ĉⱼ = it.ĉ[1,:]
	ĉₘ₁ = it.ĉ[23,:]
	ĉₘ₂ = it.ĉ[19,:]
end

# ╔═╡ 47fb48e2-c548-43b4-8727-a2ff6b6dcc06
ĉᵢ = it.ĉ[:,1]

# ╔═╡ Cell order:
# ╠═d8d0af8e-0c5e-11ed-1370-cb8d2ab22ce7
# ╠═79399bed-95e7-43ea-bef4-dd56a3364e0c
# ╠═febdeadf-fff3-4b3d-ad18-4be6658d62b7
# ╠═a72e76c9-16ba-4ecd-844d-77460f6333b0
# ╟─0193979e-496c-41f7-94ed-4d64364cb6ee
# ╟─c0184f6b-6c80-4ea3-b6e2-2aceeabeafcc
# ╟─61881095-681d-45b3-808e-b93184e24f8b
# ╠═a466c05b-9e88-4153-ab62-342497ce5ed8
# ╟─af1a7f63-9a1c-4953-9f73-56a0d357669e
# ╠═9eada2f0-d022-41f6-aa4e-b02691a7e251
# ╠═47fb48e2-c548-43b4-8727-a2ff6b6dcc06
# ╠═48da0d27-b166-4765-8566-a3d6e22223cc
# ╠═37640243-9cfe-46b5-ab5a-1b9f92dd9d29
# ╠═b83d4f0b-f104-4cbd-8a50-17bb7ce51813
