using Base:Float64
using Revise
using AnalyticalEOR
using Plots

# Initial Concentrations in mol/L
# Injection concentrations in mol/L
# CEC in mol/L
# Equilibrium constants




# * Initial concentrations
ζᵢ = [	0.47, # Na
		0.049, # Mg
		0.0115, # Ca
		]
push!(ζᵢ, ζᵢ[1] + 2ζᵢ[2] + 2ζᵢ[3]) # Anion Cl conc. by charge balance

# * Injected concentrations
ζⱼ = [	0.00166, # Na
	0.00062, # Mg
	0.00310, # Ca
	]
	
push!(ζⱼ, ζⱼ[1] + 2ζⱼ[2] + 2ζⱼ[3] ) # Anion Cl conc. by charge balance

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

ec = ExchangeConstants(K₂₁, K₃₁, K₂₃, Z)


ie = solve_Ion_Transport(ζᵢ, ζⱼ, ν, ec)


sol2, sol3 = M2_ODE_solutions(ie.cₘ₂[3], ie.cⱼ, ie.cₘ₁, ec)

plot(sol2, color="red", label="W3")
plot!(sol3, color="blue", label="W2")
plot!([ie.cᵢ[3] ie.cₘ₁[3] ie.cₘ₂[3] ie.cⱼ[3] c3],
	  [ie.cᵢ[2] ie.cₘ₁[2] ie.cₘ₂[2] ie.cⱼ[2] sol2(c3)],
	  seriestype=:scatter, labels=["I" "M1" "M2" "J"]
	 )
plot!(ylim=(1e-5, 1), xlim=(1e-6, 1),
	 scale=:log, xlabel="Ca", ylabel="Mg",
)


t = 0.9
plot(ie.λ * t, [ie.c₁ ie.c₂ ie.c₃ ie.c₄],
layout=4, label=false, lw=2.5)
plot!(xlim=(0, 1), label=false, lw=2, ylabel="Concentration, M")

plot(ie.c₃, ie.λ  )

ie
ie.c₂

tspan = collect(range(0, 20, length=shifts))

plot(1 ./ ie.λ, ie.c₂*2, yscale=:log10)


c3r = collect(range(ie.cₘ₂[3], ie.cₘ₁[3], length=50))
c3l = collect(range(ie.cⱼ[3], ie.cₘ₂[3], length=50))

ie.λ₂

plot(c3r, 1 ./ ie.λ₂)
plot!(c3l, 1 ./ ie.λ₃)
plot!([c3r[end] c3r[1] c3l[1]], [1 ./ie.λ₂[end] 1 ./ie.λ₂[1] 1 ./ie.λ₃[1]], labels=["M1" "M2" "J"], seriestype=:scatter)
plot!(xscale=:log10)


plot(c3r,  ie.λ₂)
plot!(c3l,  ie.λ₃)
plot!([c3r[end] c3r[1] c3l[1]], [ie.λ₂[end] ie.λ₂[1] ie.λ₃[1]], labels=["M1" "M2" "J"], seriestype=:scatter)
plot!(xscale=:log10)


isotherm



1 ./ ie.λ₂
1 ./ ie.λ₃
# plot!(λ * t, ĉ', xlim=(0,1), layout=4,
# 		label=false, lw=2,)
# plot!(ylabel="Concentration, mol/l", xlabel="x", xguidefontsize=8,
# 	yguidefontsize=8, titlefontsize=10, legendtitlefontsize=8
# 	)
# plot!(λ * t, c+ĉ', xlim=(0,1), layout=4,
# 		title=["Na" "Mg" "Ca" "Cl"], label=false, lw=2, ls=:dot,
# 		ylabel="Concentration, mol/l")

# plot!(λ * (t2) / 100, c, layout=4, labels=false, lw=2)