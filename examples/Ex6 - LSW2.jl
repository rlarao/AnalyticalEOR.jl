using AnalyticalEOR
using Plots

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

plot!(σ, c[:,2], yscale=:log10, ylim=(1e-5, 1), xlim=(0, 20),
		title=["Mg"], label=false, lw=2, ylabel="Concentration, M", marker=:circle)



t = 0.5
plot(λ * t, c, xlim=(0, 1), layout=4,
		title=["Na" "Mg" "Ca" "Cl"], label=false, lw=2, ylabel="Concentration, M", marker=:circle)

c[:,2]
λ

cᵢ = ζᵢ .* ν
cⱼ = ζⱼ .* ν

# * First intermediate point
ĉᵢ = isotherm(cᵢ, ec)
ĉₘ₁ = ĉᵢ
cₘ₁ = flowingConcentrations(ĉₘ₁, cⱼ[4], ec)

# * Second intermediate point
cₘ₂ = solve_IntegralCurve(cₘ₁, cⱼ, ec)

# * Get intermediate points
sol2, sol3 = M2_ODE_solutions(cₘ₂[3], cⱼ, cₘ₁, ec)


σ₁ = 1
σ₂ₘ₂ = eigenvectors([cₘ₂[2] cₘ₂[3] cⱼ[4]], ec)[1]
σ₂ₘ₁ = eigenvectors([cₘ₁[2] cₘ₁[3] cⱼ[4]], ec)[1] 
σ₃ⱼ =  eigenvectors([cⱼ[2] cⱼ[3] cⱼ[4]], ec)[2] 
σ₃ₘ₂ = eigenvectors([cₘ₂[2] cₘ₂[3] cⱼ[4]], ec)[2]

σ₃[1]
σ₃[end]

if σ₃ₘ₂ >= σ₃ⱼ
    𝒲₃ = :shock
else
    𝒲₃ = :rarefication
end

if σ₂ₘ₁ >= σ₂ₘ₂
    𝒲₂ = :shock
else
    𝒲₂ = :rarefication
end

c₃ᵣ = collect(range(cₘ₂[3], cₘ₁[3],  length=100))[2:end - 1]
c₃ₗ = collect(range(cⱼ[3], cₘ₂[3],  length=100))[2:end - 1]

function rankine_hugoniot(c₂, c₃, cₘ₁, ec)
    ĉ₁, ĉ₂, ĉ₃, ĉ₄ = isotherm([c₂, c₃, cₘ₁[4]], ec)
    ĉ₁ₘ₁, ĉ₂ₘ₁, ĉ₃ₘ₁, ĉ₄ₘ₁ = isotherm(cₘ₁, ec)

    LHS = ((ĉ₂ + c₂) - (ĉ₂ₘ₁ + cₘ₁[2])) / (c₂ - cₘ₁[2])
    RHS = ((ĉ₃ + c₃) - (ĉ₃ₘ₁ + cₘ₁[3])) / (c₃ - cₘ₁[3])

    return LHS - RHS
end

function RH_eigenvalues(c₂, c₃, cₗ, ec)
    ĉ₁, ĉ₂, ĉ₃, ĉ₄ = isotherm([c₂, c₃, cₗ[4]], ec)
    ĉ₁ₗ, ĉ₂ₗ, ĉ₃ₗ, ĉ₄ₗ = isotherm(cₗ, ec)

    return ((ĉ₃ + c₃) - (ĉ₃ₗ + cₗ[3])) / (c₃ - cₗ[3])
end


loss(c₂, c₃) = rankine_hugoniot(c₂, c₃, cₘ₁, ec)
loss2(c₂, c₃) = rankine_hugoniot(c₂, c₃, cⱼ, ec)

c₂ᵣ = [find_zero(c -> loss(c, c₃), [cₘ₁[2], cₘ₂[2] + 1e-3]) for c₃ in c₃ᵣ]
c₂ₗ = [find_zero(c -> loss2(c, c₃), [cₘ₂[2], cⱼ[2] - 1e-3]) for c₃ in c₃ₗ]



plot(sol3, color="blue", label="W3")
plot!(sol2, color="red", label="W2")
plot!([cᵢ[3] cₘ₁[3] cₘ₂[3] cⱼ[3]],
[cᵢ[2] cₘ₁[2] cₘ₂[2] cⱼ[2]],
seriestype=:scatter, labels=["I" "M1" "M2" "J"]
)
plot!(ylim=(1e-5, 1), xlim=(1e-6, 1),
scale=:log, xlabel="Ca", ylabel="Mg",

)


plot(c₃ᵣ, c₂ᵣ, lw=3, label="Shock2", scale=:log10)
plot!(c₃ₗ, c₂ₗ, lw=3, label="Shock3", scale=:log10)
plot!(sol2, color="red", label="W2")
plot!(sol3, color="red", label="W3")
plot!([cᵢ[3] cₘ₁[3] cₘ₂[3] cⱼ[3]],
[cᵢ[2] cₘ₁[2] cₘ₂[2] cⱼ[2]],
seriestype=:scatter, labels=["I" "M1" "M2" "J"])
plot!(ylim=(1e-5, 1), xlim=(1e-6, 1),
scale=:log, xlabel="Ca", ylabel="Mg",)


σ̃₂ = [RH_eigenvalues(c₂, c₃, cₘ₁, ec) for (c₂, c₃) in zip(c₂ᵣ, c₃ᵣ)]
σ̃₃ = [RH_eigenvalues(c₂, c₃, cⱼ, ec) for (c₂, c₃) in zip(c₂ₗ, c₃ₗ)]

λ̃₂ = 1 ./ σ̃₂
λ̃₃ = 1 ./ σ̃₃


σ₂ = [eigenvectors([c₂ c₃ cⱼ[4]], ec)[1] for (c₂, c₃) in zip(c₂ᵣ, c₃ᵣ)]
σ₃ = [eigenvectors([c₂ c₃ cⱼ[4]], ec)[2] for (c₂, c₃) in zip(c₂ₗ, c₃ₗ)]
λ₂ = 1 ./ σ₂
λ₃ = 1 ./ σ₃

plot(c₃ᵣ,  λ̃₂, label="S2")
plot!(c₃ₗ,  λ̃₃, label="S3")
plot!(c₃ᵣ,  λ₂, label="R2")
plot!(c₃ₗ,  λ₃, label="R3")
plot!([cₘ₁[3] cₘ₂[3] cⱼ[3]],
[cₘ₁[2] cₘ₂[2] cⱼ[2]],
seriestype=:scatter, labels=["M1" "M2" "J"]
)



plot(c₃ᵣ,  σ̃₂, label="S2")
plot!(c₃ₗ,  σ̃₃, label="S3")
plot!(c₃ᵣ,  σ₂, label="R2")
plot!(c₃ₗ,  σ₃, label="R3")
plot!(yscale=:log10)


σ̃₂ₘ₂ =  RH_eigenvalues(cₘ₂[2], cₘ₂[3], cₘ₁, ec)
σ̃₃ⱼ =  RH_eigenvalues(cⱼ[2], cⱼ[3], cₘ₂, ec)
λ̃₂ₘ₂ = 1 ./ σ̃₂ₘ₂ 
λ̃₃ⱼ = 1 ./ σ̃₃ⱼ 


c₂ = [cⱼ[2];  cⱼ[2]; cₘ₂[2]  ; cₘ₂[2] ; cₘ₁[2]; cₘ₁[2];  cᵢ[2]; cᵢ[2]]
c₃ = [cⱼ[3];  cⱼ[3]; cₘ₂[2]  ; cₘ₂[2] ; cₘ₁[3]; cₘ₁[3];  cᵢ[3]; cᵢ[3]]
c₄ = [cⱼ[4];  cⱼ[4]; cⱼ[4]   ; cⱼ[4]  ;  cₘ₁[4]; cₘ₁[4]; cᵢ[4]; cᵢ[4]]
λ =  [ 0.0001  ;  λ̃₃ⱼ   ; λ̃₃ⱼ  ;  λ̃₂ₘ₂  ;  λ̃₂ₘ₂  ;  1    ;   1 ;   100 ] 
c₁ = c₄ .- c₃ .- c₂

c = [c₁ c₂ c₃ c₄]

t = 0.5
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

