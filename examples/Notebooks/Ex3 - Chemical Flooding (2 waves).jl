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

# ╔═╡ ce2623b0-0bba-11ed-0bb4-37c4e85dc528
# ╠═╡ show_logs = false
begin
	import Pkg
	Pkg.activate("/home/rick/.julia/dev/AnalyticalEOR/Project.toml")
	using Revise, AnalyticalEOR, PlutoUI, Plots, Roots
end

# ╔═╡ f48e4e8f-3b3c-402d-ad54-084f1d61010c
md""" ### Input"""

# ╔═╡ 6dd2d969-64c4-4ba7-a90f-a4e7e0a17bae
begin
    kr = RelPerms(swr=0.2,
                    sor=0.2,
                    krw0=0.2,
                    kro0=0.5,
                    nw= 2.0,
                    no=3.5)

    μw = 1.0
    μo = 5.0

    wf = WaterFlooding( kr=kr,
                        μw=μw,
                        μo=μo
                        )
end

# ╔═╡ 5f7ef8a1-37f0-4996-9018-3640514cd44d
@bind D1 Slider(0.:0.005:4, default=0.5)

# ╔═╡ fdfc1b64-6821-4f98-b1b7-e9d277e93009
md"""D1 = $D1"""

# ╔═╡ 44ff85ee-6968-4fd2-8e29-677d13074437
@bind D2 Slider(0.:0.005:D1, default=D1*0.5)

# ╔═╡ 7d515eed-e29c-4df6-a329-066b4315e19c
begin
    kr2 = RelPerms(swr=0.2,
                    sor=0.2,    
                    krw0=0.5,
                    kro0=1.0,
                    nw=3.0,
                    no=2.0)

    kr3 = RelPerms(swr=0.2,
                    sor=0.2,    
                    krw0=0.2,
                    kro0=1.0,
                    nw=3.0,
                    no=2.0)

    cf = ChemicalFlooding(
                    wf = wf,
                    krs = [kr3, kr2],
                    D = [D1, D2]
                    )
end

# ╔═╡ 4e11b8ef-89f7-4269-a2d3-0a05252a4c6f
md"""D2 = $D2"""

# ╔═╡ 534b1239-55da-45d9-9724-bc87179b322f
@bind si Slider(0.2:0.01:0.8, default=0.2)

# ╔═╡ 30b033cb-4272-4d21-b02d-25fcef74c569
md"""si = $si"""

# ╔═╡ 7b31a9f7-d6c8-4c23-a1be-a365d1035126
@bind sj Slider(0.2:0.01:0.8, default=0.8)

# ╔═╡ 7d7c06aa-02ec-4216-81d1-5e1a6296c17d
md"""sj = $sj"""

# ╔═╡ 854b12ea-383e-4db3-8d6f-160752bca530
# ╠═╡ show_logs = false
begin
    sol = solve_cf(cf, si, sj)
    plot_fractional_flow(cf, sol)
end

# ╔═╡ 0505ae15-e94c-49eb-b077-e50383ef0996
@bind t Slider(0.01:0.01:1.5, default=0.2)

# ╔═╡ 51d4444e-6837-401c-afb9-c157e9749e8f
md"""t = $t"""

# ╔═╡ 65daaf31-fdff-4b56-87f0-d70e7858d662
begin
    plot_sat_profile(sol, t)
end

# ╔═╡ Cell order:
# ╟─ce2623b0-0bba-11ed-0bb4-37c4e85dc528
# ╟─f48e4e8f-3b3c-402d-ad54-084f1d61010c
# ╠═6dd2d969-64c4-4ba7-a90f-a4e7e0a17bae
# ╠═7d515eed-e29c-4df6-a329-066b4315e19c
# ╟─fdfc1b64-6821-4f98-b1b7-e9d277e93009
# ╟─5f7ef8a1-37f0-4996-9018-3640514cd44d
# ╟─4e11b8ef-89f7-4269-a2d3-0a05252a4c6f
# ╟─44ff85ee-6968-4fd2-8e29-677d13074437
# ╟─30b033cb-4272-4d21-b02d-25fcef74c569
# ╟─534b1239-55da-45d9-9724-bc87179b322f
# ╟─7d7c06aa-02ec-4216-81d1-5e1a6296c17d
# ╟─7b31a9f7-d6c8-4c23-a1be-a365d1035126
# ╟─854b12ea-383e-4db3-8d6f-160752bca530
# ╟─51d4444e-6837-401c-afb9-c157e9749e8f
# ╟─0505ae15-e94c-49eb-b077-e50383ef0996
# ╟─65daaf31-fdff-4b56-87f0-d70e7858d662
