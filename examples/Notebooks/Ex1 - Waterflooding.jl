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

# ╔═╡ 2858e2e6-07dd-11ed-1828-1d526ab9abb9
begin
	import Pkg
	Pkg.activate("/home/rick/.julia/dev/AnalyticalEOR/Project.toml")
	using AnalyticalEOR, PlutoUI, Plots
end

# ╔═╡ 06e65237-45da-458e-811a-764920a4526e
md"""# Input"""

# ╔═╡ 03be79da-6828-4986-8a4f-1c92254be9ba
@bind nw Slider(1.:0.1:3.5, default=1.5)

# ╔═╡ 0798bba6-bafa-41ce-ad22-2fca894541eb
md"""nw = $nw"""

# ╔═╡ e9f5326f-98b3-4d51-8a89-07e487d0de5a
@bind no Slider(1.0:0.1:3.5, default=1.5)

# ╔═╡ 343d45e0-8762-43c5-abb4-84f0223521c1
begin
    kr = RelPerms(swr=0.2,
                    sor=0.2,
                    krw0=0.2,
                    kro0=0.5,
                    nw= nw,
                    no=no)

    μw = 1.0
    μo = 5.0

	wf = WaterFlooding(
					kr=kr,
					μw=μw,
					μo=μo
	)
end

# ╔═╡ cd331ff4-8945-4867-9d97-132008bdc793
md"""no = $no"""

# ╔═╡ 1f6ee5d0-bf26-47f2-a5f9-677b226a927b
@bind si Slider(0.2:0.01:0.8, default=0.2)

# ╔═╡ 8521cf42-a493-4c90-8ea1-fbd346337e62
md"""si = $si"""

# ╔═╡ bb13b791-46b4-4a38-b851-b30e598667d7
@bind sj Slider(0.2:0.01:0.8, default = 0.8)

# ╔═╡ e5583896-6876-41f3-826c-0e1918ac7371
md"""sj = $sj"""

# ╔═╡ aedf9d8f-5aa8-47a5-a4af-b0d191898424
# ╠═╡ show_logs = false
begin
   	sol = solve_wf(wf, si, sj)
    plot_fractional_flow(wf, sol)
	ylims!(0,1)
end

# ╔═╡ a7157322-a5d9-4374-b38d-7c2baabd73e5
@bind t Slider(0.01:0.01:1.5, default=0.2)

# ╔═╡ e05a97f3-a3ef-473f-aaaf-4aeac2a80ce2
md"""t = $t"""

# ╔═╡ a657c6d2-3ec4-4043-a811-8537b05fb2af
begin
    plot_sat_profile(sol, t)
end

# ╔═╡ Cell order:
# ╠═2858e2e6-07dd-11ed-1828-1d526ab9abb9
# ╟─06e65237-45da-458e-811a-764920a4526e
# ╠═343d45e0-8762-43c5-abb4-84f0223521c1
# ╟─0798bba6-bafa-41ce-ad22-2fca894541eb
# ╟─03be79da-6828-4986-8a4f-1c92254be9ba
# ╟─cd331ff4-8945-4867-9d97-132008bdc793
# ╟─e9f5326f-98b3-4d51-8a89-07e487d0de5a
# ╟─8521cf42-a493-4c90-8ea1-fbd346337e62
# ╟─1f6ee5d0-bf26-47f2-a5f9-677b226a927b
# ╟─e5583896-6876-41f3-826c-0e1918ac7371
# ╟─bb13b791-46b4-4a38-b851-b30e598667d7
# ╟─aedf9d8f-5aa8-47a5-a4af-b0d191898424
# ╟─e05a97f3-a3ef-473f-aaaf-4aeac2a80ce2
# ╟─a7157322-a5d9-4374-b38d-7c2baabd73e5
# ╟─a657c6d2-3ec4-4043-a811-8537b05fb2af
