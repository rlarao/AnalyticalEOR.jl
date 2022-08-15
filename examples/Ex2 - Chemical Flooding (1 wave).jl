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

# ╔═╡ 35ac77e6-0ba7-11ed-1ad8-d9356b568460
begin
	import Pkg
	Pkg.activate("/home/rick/.julia/dev/AnalyticalEOR/Project.toml")
	using Revise, AnalyticalEOR, PlutoUI, Plots, Roots
end

# ╔═╡ db885e1c-569c-451a-a32b-21f50df62ca9
@bind si Slider(0.2:0.01:0.8, default=0.2)

# ╔═╡ 54cfd1de-b5cd-455c-9228-cd3e8ee26858
md"""si = $si"""

# ╔═╡ d8e0dc30-55cb-4961-bc10-a56433387af7
@bind sj Slider(0.2:0.01:0.8, default=0.8)

# ╔═╡ 945875a7-696c-4a5d-9252-1da652176f4f
md"""sj = $sj"""

# ╔═╡ 1e0de9ea-eeae-4879-9626-1c8020ff3436
@bind D Slider(0.:0.05:4.0, default=0.9)

# ╔═╡ 0a6e7335-1a57-4ae2-bbde-2c9de30c3f6e
begin
    kr = RelPerms(swr=0.2,
                    sor=0.2,
                    krw0=0.2,
                    kro0=0.5,
                    nw= 2.0,
                    no=3.5)

    μw = 1.0
    μo = 5.0

    wf = WaterFlooding(
                        kr=kr,
                        μw=μw,
                        μo=μo
                        )

    kr_2 = RelPerms(swr=0.2,
                    sor=0.2,    
                    krw0=0.5,
                    kro0=1.0,
                    nw=3.0,
                    no=2.0)
    
    cf = ChemicalFlooding(
                    wf = wf,
                    krs = [kr_2],
                    D = [D]
                    )

end

# ╔═╡ 29d661f4-3df7-4806-955d-b56965a46c60
md"""D = $D"""

# ╔═╡ f36f3a96-c2c8-4d14-8a10-efef0b4caedc
# ╠═╡ show_logs = false
begin
    sol = solve_cf(cf, si, sj)
    plot_fractional_flow(cf, sol)
end

# ╔═╡ 667d9bea-9495-4c1b-9634-9cf397d3b7cf
@bind t Slider(0.01:0.01:1.5, default=0.2)

# ╔═╡ 610b8346-b0cd-4116-a06c-b8e0c4325858
begin
    plot_sat_profile(sol, t)
end

# ╔═╡ Cell order:
# ╠═35ac77e6-0ba7-11ed-1ad8-d9356b568460
# ╠═0a6e7335-1a57-4ae2-bbde-2c9de30c3f6e
# ╟─54cfd1de-b5cd-455c-9228-cd3e8ee26858
# ╟─db885e1c-569c-451a-a32b-21f50df62ca9
# ╟─945875a7-696c-4a5d-9252-1da652176f4f
# ╟─d8e0dc30-55cb-4961-bc10-a56433387af7
# ╟─29d661f4-3df7-4806-955d-b56965a46c60
# ╠═1e0de9ea-eeae-4879-9626-1c8020ff3436
# ╟─f36f3a96-c2c8-4d14-8a10-efef0b4caedc
# ╠═610b8346-b0cd-4116-a06c-b8e0c4325858
# ╟─667d9bea-9495-4c1b-9634-9cf397d3b7cf
