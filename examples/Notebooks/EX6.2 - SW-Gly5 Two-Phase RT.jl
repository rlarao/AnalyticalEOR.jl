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

# ╔═╡ a06e0ff4-9540-4985-8a24-0d2f20b5d3ef
md""" #### Single Phase Reactive Transport"""

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
        0.1, # Gly
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
# plot!(ylims=(1e-5, 1e-3))
end

# ╔═╡ a466c05b-9e88-4153-ab62-342497ce5ed8
begin
	plot_flowing_conc(it, t, c_labels, c_colors)
	# plot!(yscale=:log10)
end

# ╔═╡ f56bcd18-29a1-4f0c-bc94-9a63f86cf750
begin
	ĉᵢ = it.ĉ[7,:]
	ĉⱼ = it.ĉ[2,:]
	ĉₘ₁ = it.ĉ[6,:]
	ĉₘ₂ = it.ĉ[3,:]
end

# ╔═╡ 3a3a2022-f96d-40b6-8457-8ce159ffdefb
md""" #### Two-Phase Flow Solution"""

# ╔═╡ 71d08b96-bb9d-4be5-85cc-85a542476a81
md"""Wave velocities"""

# ╔═╡ 5882ab82-e458-4fc8-ad69-387bb3d20f59
begin
	λᵢ = it.λ[7]
	λⱼ = it.λ[2]
	λₘ₁ = it.λ[5]
	λₘ₂ = it.λ[3]
end

# ╔═╡ 68ea7ef4-3867-4498-9be2-1c8a79a31904
md""" Eigenvalues"""

# ╔═╡ a8efdad5-4dbf-4f35-9c06-82f6a4445730
begin
	σᵢ = it.σ[7]
	σⱼ = it.σ[2]
	σₘ₁ = it.σ[5]
	σₘ₂ = it.σ[3]
end

# ╔═╡ 7563ec36-dadd-4294-aec3-7f330a424686
begin
	σ₁ = σₘ₁ + 1
	σ₂ = σₘ₂ + 1
end

# ╔═╡ 4e9bb37c-c505-42c5-9a76-39c694c68727
begin
	λ₁ = 1 / σ₁
	λ₂ = 1 / σ₂
end

# ╔═╡ e89c7dbe-f64f-4618-af85-50c6cbb132b1
md""" Calculate wettability index"""

# ╔═╡ a3f49e8e-fc91-4131-afd2-c44741beb01e
begin
	coo_ww = ĉⱼ[1]
	coo_ow = ĉᵢ[1]
	coo_mw = ĉₘ₂[1]

	ω = (coo_mw - coo_ww) / (coo_ow - coo_ww)
end

# ╔═╡ 7ffb2390-c065-43f4-bab1-c18eaf2c64f0
@bind no_ww Slider(1.0:0.1:3.5, default=2.)

# ╔═╡ 92b8430a-9a58-47c9-856e-703a8d8d558f
println("no_ww = ", no_ww)

# ╔═╡ d5e97d51-e14a-43cb-a408-80f42c74900c
@bind no_ow Slider(1.0:0.1:3.5, default=3.0)

# ╔═╡ eb310cd1-3143-4932-a37c-1d9d8a3a4ecc
begin
	kr_ww = RelPerms(swr=0.11,
	                sor=0.10,
	                krw0=0.2,
	                kro0=0.8,
	                nw=3.,
	                no=no_ww)
	
	kr_ow = RelPerms(swr=0.11,
	                sor=0.16,
	                krw0=0.4,
	                kro0=0.5,
	                nw=2.,
	                no=no_ow)
end

# ╔═╡ c38d2695-de63-417e-b94e-35b0fd8da214
begin
	sor(ω) = (1-ω) * kr_ww.sor + ω * kr_ow.sor
	krw0(ω) = (1-ω) * kr_ww.krw0 + ω * kr_ow.krw0
	kro0(ω) = (1-ω) * kr_ww.kro0 + ω * kr_ow.kro0
	no(ω) = (1-ω) * kr_ww.no + ω * kr_ow.no
	nw(ω) = (1-ω) * kr_ww.nw + ω * kr_ow.nw
	
	kr_mw = RelPerms(swr=0.11,
	                sor=sor(ω),
	                krw0=krw0(ω),
	                kro0=kro0(ω),
	                nw=nw(ω),
	                no=no(ω))
end

# ╔═╡ 25f1e16c-0eae-48fa-b78f-c2783fcf655e
println("no_ow = ", no_ow)

# ╔═╡ 4cdae795-360b-4f21-893c-d1932d995410
@bind μo Slider(5.:1.:300., default=5.)

# ╔═╡ 1febe860-956e-4d2d-afc5-400ac7ab0119
begin
   	μw = 1.0
    # μo = μo

    wf = WaterFlooding( kr=kr_ow,
                        μw=μw,
                        μo=μo,
                        )
end

# ╔═╡ f65d64fe-230b-49c6-85ba-0dca86c198ce
begin
    cf = ChemicalFlooding(
                    wf = wf,
                    krs = [kr_ww, kr_mw],
                    D = [λ₁, λ₂]
                    )
end

# ╔═╡ 64c2766a-cc40-4cf7-9506-873af608361e
md""" μo = $μo"""

# ╔═╡ b718c546-3a69-40a2-9f9f-eeb144ee2ba8
@bind si Slider(0.11:0.01:0.9, default=0.2)

# ╔═╡ 871b9f85-c90d-41e4-bc7c-0934fdff4fbf
md"""si = $si"""

# ╔═╡ 18dac1ae-cfe1-4341-9e7f-368537ea81e3
@bind sj Slider(0.11:0.01:0.9, default=0.8)

# ╔═╡ 2ea00e43-ec3a-4e59-a57a-768dd1e95e13
md"""sj = $sj"""

# ╔═╡ 3744cb12-ce8c-4dd1-8bd9-d6cdfe5c59d8
# ╠═╡ show_logs = false
begin
	sol = solve_cf(cf, si, sj)
    plot_fractional_flow(cf, sol)
end

# ╔═╡ 67a62549-2b20-4819-80ca-7ae0aa0f2dc1
@bind t2 Slider(0.01:0.01:1.5, default=0.2)

# ╔═╡ 66294f94-5be6-4df6-b992-03d3ec29a7e7
md"""t = $t2"""

# ╔═╡ ad7ca85f-43b7-471a-ac7f-fefb7ec9b4c0
    plot_sat_profile(sol, t2)

# ╔═╡ Cell order:
# ╟─d8d0af8e-0c5e-11ed-1370-cb8d2ab22ce7
# ╟─79399bed-95e7-43ea-bef4-dd56a3364e0c
# ╟─a06e0ff4-9540-4985-8a24-0d2f20b5d3ef
# ╠═febdeadf-fff3-4b3d-ad18-4be6658d62b7
# ╠═a72e76c9-16ba-4ecd-844d-77460f6333b0
# ╟─0193979e-496c-41f7-94ed-4d64364cb6ee
# ╟─c0184f6b-6c80-4ea3-b6e2-2aceeabeafcc
# ╟─61881095-681d-45b3-808e-b93184e24f8b
# ╟─a466c05b-9e88-4153-ab62-342497ce5ed8
# ╟─f56bcd18-29a1-4f0c-bc94-9a63f86cf750
# ╟─3a3a2022-f96d-40b6-8457-8ce159ffdefb
# ╟─71d08b96-bb9d-4be5-85cc-85a542476a81
# ╠═5882ab82-e458-4fc8-ad69-387bb3d20f59
# ╟─68ea7ef4-3867-4498-9be2-1c8a79a31904
# ╠═a8efdad5-4dbf-4f35-9c06-82f6a4445730
# ╠═7563ec36-dadd-4294-aec3-7f330a424686
# ╠═4e9bb37c-c505-42c5-9a76-39c694c68727
# ╟─e89c7dbe-f64f-4618-af85-50c6cbb132b1
# ╠═a3f49e8e-fc91-4131-afd2-c44741beb01e
# ╠═eb310cd1-3143-4932-a37c-1d9d8a3a4ecc
# ╠═c38d2695-de63-417e-b94e-35b0fd8da214
# ╠═1febe860-956e-4d2d-afc5-400ac7ab0119
# ╠═f65d64fe-230b-49c6-85ba-0dca86c198ce
# ╟─92b8430a-9a58-47c9-856e-703a8d8d558f
# ╟─7ffb2390-c065-43f4-bab1-c18eaf2c64f0
# ╟─25f1e16c-0eae-48fa-b78f-c2783fcf655e
# ╟─d5e97d51-e14a-43cb-a408-80f42c74900c
# ╟─64c2766a-cc40-4cf7-9506-873af608361e
# ╟─4cdae795-360b-4f21-893c-d1932d995410
# ╟─871b9f85-c90d-41e4-bc7c-0934fdff4fbf
# ╟─b718c546-3a69-40a2-9f9f-eeb144ee2ba8
# ╟─2ea00e43-ec3a-4e59-a57a-768dd1e95e13
# ╟─18dac1ae-cfe1-4341-9e7f-368537ea81e3
# ╟─3744cb12-ce8c-4dd1-8bd9-d6cdfe5c59d8
# ╟─66294f94-5be6-4df6-b992-03d3ec29a7e7
# ╟─67a62549-2b20-4819-80ca-7ae0aa0f2dc1
# ╟─ad7ca85f-43b7-471a-ac7f-fefb7ec9b4c0
