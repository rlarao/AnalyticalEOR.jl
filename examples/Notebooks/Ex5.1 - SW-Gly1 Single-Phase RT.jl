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
	using Revise, AnalyticalEOR, PlutoUI, Plots, Roots, LaTeXStrings
end

# ╔═╡ 79399bed-95e7-43ea-bef4-dd56a3364e0c
begin
    c_colors = ["brown" "darkorchid" "orangered" "gold"]
    c_labels = [L"COO^-" L"SO_4^{2-}" L"Gly^-" L"Na^+"]
    a_colors = ["brown" "darkorchid" "orangered"]
    a_labels = [L"COOX" L"SO_4X_2" L"GlyX"]
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

	ec = ExchangeConstants(K₂₁, K₃₁, K₂₃, Z, ν)

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

# ╔═╡ Cell order:
# ╠═d8d0af8e-0c5e-11ed-1370-cb8d2ab22ce7
# ╠═79399bed-95e7-43ea-bef4-dd56a3364e0c
# ╠═febdeadf-fff3-4b3d-ad18-4be6658d62b7
# ╠═a72e76c9-16ba-4ecd-844d-77460f6333b0
# ╟─0193979e-496c-41f7-94ed-4d64364cb6ee
# ╟─c0184f6b-6c80-4ea3-b6e2-2aceeabeafcc
# ╟─61881095-681d-45b3-808e-b93184e24f8b
# ╠═a466c05b-9e88-4153-ab62-342497ce5ed8
