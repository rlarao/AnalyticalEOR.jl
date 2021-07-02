using AnalyticalEOR
using Plots

# * Initial concentrations
ζᵢ = [	3.5553, # Na mol/L
		0.1213, # Mg
		0.5481, # Ca
		]


push!(ζᵢ, ζᵢ[1] + 2ζᵢ[2] + 2ζᵢ[3]) # Anion Cl conc. by charge balance

# * Injected concentrations
ζⱼ = [	1.00831, # Na
	0.0223, # Mg
	0.0726488, # Ca
	]
push!(ζⱼ, ζⱼ[1] + 2ζⱼ[2] + 2ζⱼ[3] ) # Anion Cl conc. by charge balance

ν =  [1, 2, 2, 1] # charges

# * Equilibrium constants and 
K₂₁ = 10^1.67148
K₃₁ = 10^1.83148     
K₂₃ =  10^-0.16

# * Cation Exchange Capacity
ρ = 2.8082271 # g/cm3
ϕ = 0.59
cec = 0.06 # 
Z = cec * ((1 - ϕ) / ϕ) * ρ # Conversion of cation exchange capacity into moles/liter

ec = ExchangeConstants(K₂₁, K₃₁, K₂₃, Z)


it = solve_Ion_Transport(ζᵢ, ζⱼ, ν, ec)

c = it.c
λ = it.λ
σ = it.σ

plot(σ.+1, c[:,2], yscale=:log10, ylim=(1e-5, 1), xlim=(0, 20),
		title=["Mg"], label=false, lw=2, ylabel="Concentration, M", marker=:circle)

# plot!(σ, c[:,2], yscale=:log10, ylim=(1e-5, 1), xlim=(0, 20),
# 		title=["Mg"], label=false, lw=2, ylabel="Concentration, M", marker=:circle)


it.W2
it.W3

t = 0.5
plot(λ * t, c, xlim=(0, 1), layout=4,
		title=["Na" "Mg" "Ca" "Cl"], label=false, lw=2, ylabel="Concentration, M", marker=:circle)

plot(λ * t, c, xlim=(0, 1), layout=4,
		title=["Na" "Mg" "Ca" "Cl"], label=false, lw=2, ylabel="Concentration, M", marker="circle")


plot(λ * t, c₃, xlim=(0, 1),
		title=["Ca"], label=false, lw=2, ylabel="Concentration, M", marker=:circle)


plot(λ * t, c₂, xlim=(0, 1),
		title=["Mg"], label=false, lw=2, ylabel="Concentration, M", marker=:circle)








σ₃ⱼ =  eigenvectors([cⱼ[2] cⱼ[3] cⱼ[4]], ec)[2] 
σ₃ₘ₂ = eigenvectors([cₘ₂[2] cₘ₂[3] cⱼ[4]], ec)[2]

c₃ₗ = collect(range(cⱼ[3], cₘ₂[3],  length=100))[2:end - 1]
σ₃ = [eigenvectors([c₂ c₃ cⱼ[4]], ec)[2] for (c₂, c₃) in zip(c₂ₗ, c₃ₗ)]

