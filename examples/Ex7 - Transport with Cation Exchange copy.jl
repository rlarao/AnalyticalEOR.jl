using Base:Float64
using Revise
using AnalyticalEOR

# Initial Concentrations in mol/L
# Injection concentrations in mol/L
# CEC in mol/L
# Equilibrium constants

# * Initial concentrations
ζᵢ = [	1e-9, # Na
		1e-9, # Mg
		5e-3, # Ca
		]
push!(ζᵢ, ζᵢ[1] + 2ζᵢ[2] + 2ζᵢ[3]) # Anion Cl conc. by charge balance

# * Injected concentrations
ζⱼ = [	95e-3, 	# Na
		2.6e-3, # Mg
		1e-9, 	# Ca
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

ec = ExchangeConstants(K₂₁, K₃₁, K₂₃, Z, ν)


ie = solve_Ion_Transport(ζᵢ, ζⱼ, ec)

using Plots
t = 0.3
plot(ie.λ * t, [ie.c₁ ie.c₂ ie.c₃ ie.c₄],
layout=4, label=false, lw=2.5)
plot!(xlim=(0, 1), label=false, lw=2, ylabel="Concentration, M")



# plot!(λ * t, ĉ', xlim=(0,1), layout=4,
# 		label=false, lw=2,)
# plot!(ylabel="Concentration, mol/l", xlabel="x", xguidefontsize=8,
# 	yguidefontsize=8, titlefontsize=10, legendtitlefontsize=8
# 	)
# plot!(λ * t, c+ĉ', xlim=(0,1), layout=4,
# 		title=["Na" "Mg" "Ca" "Cl"], label=false, lw=2, ls=:dot,
# 		ylabel="Concentration, mol/l")

# plot!(λ * (t2) / 100, c, layout=4, labels=false, lw=2)