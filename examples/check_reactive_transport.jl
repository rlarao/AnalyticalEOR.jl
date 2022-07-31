using Revise, AnalyticalEOR

# Initial Concentrations in mol/L
# Injection concentrations in mol/L
# CEC in mol/L
# Equilibrium constants

# * Initial concentrations
begin
	ν = [1, 1, 1, 1] # charges

	ζᵢ = [	0.0017, # Na
			0.00124, # Mg
			0.0062, # Ca
		]
    push!(ζᵢ, ν[1] * ζᵢ[1] + ν[2] * ζᵢ[2] + ν[3] * ζᵢ[3])
	
	# * Injected concentrations
	ζⱼ = [	0.47, 	# Na
			0.098, # Mg
			0.023, 	# Ca
		]
		
    push!(ζⱼ, ν[1] * ζⱼ[1] + ν[2] * ζⱼ[2] + ν[3] * ζⱼ[3])

    Z = 0.1
end
# * Equilibrium constants and 
begin
    # K₂₁ = 10^2.8 / (2Z)
    # K₃₁ = 10^1.78 / (2Z)
    K₂₁ = 10^2.8 
    K₃₁ = 10^1.78 
    K₂₃ = K₂₁ / K₃₁ 

# * Cation Exchange Capacity
# Z = 0.1
# ec = IonExchangeProblem(K₂₁, K₃₁, K₂₃, Z, ν)
end


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