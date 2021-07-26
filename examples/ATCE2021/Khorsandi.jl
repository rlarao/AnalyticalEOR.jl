using AnalyticalEOR
using Plots
using Roots

ν =  [1, 1, 2, 1] # charges
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

it = solve_Ion_Transport(ζᵢ, ζⱼ, ec)

plot_ODEs(it)
plot_velocities(it)

it.W2
it.W3

c = it.c
ĉ = it.ĉ
λ = it.λ
σ = it.σ

t = 0.9
plot(λ * t, c, xlim=(0, 1), layout=4,
		title=["COO" "S" "Gly" "Na"], label=false, lw=2, ylabel="Concentration, M", marker=:circle)

plot(λ * t, ĉ[:,1:3], xlim=(0, 1),
		labels=["COO" "S" "Gly"], label=false, lw=2, ylabel="Concentration, M", marker=:circle,
        colors=[:brown :darkorchid :orangered])
        plot!(ylims=(1e-4, 1), yscale=:log10)

        plot(λ * t, ĉ, xlim=(0, 1), layout=4,
		title=["COO" "S" "Gly" "Na"], label=false, lw=2, ylabel="Concentration, M", marker=:circle)


cᵢ = ζᵢ .* ν
cⱼ = ζⱼ .* ν
ĉᵢ = isotherm(cᵢ, ec)
ĉₘ₁ = ĉᵢ
cₘ₁ = flowingConcentrations(ĉₘ₁, cⱼ[4], ec)

c₃ₘ₂ = collect(range(cⱼ[3], cₘ₁[3], length=100000))
c₃₂ = c₃ₘ₂[1]
sol2 = M2_ODE2(c₃₂, cⱼ, cₘ₁, ec)

c₃₁ = c₃ₘ₂[end]
sol3 = M2_ODE3(c₃₁, cⱼ, cₘ₁, ec)


plot(sol2)
plot!(sol3, )
plot!([ cₘ₁[3] cⱼ[3]],
[ cₘ₁[2]  cⱼ[2]], seriestype=:scatter, labels=["M1" "J"],  legend=:outerright)

c₃ₘ₂ = fzero(c -> sol2(c) - sol3(c), 0.002)
c₂ₘ₂ = sol3(c₃ₘ₂)
cₘ₂ = [c₂ₘ₂, c₃ₘ₂, cⱼ[4]]
prepend!(cₘ₂, cₘ₂[3] - cₘ₂[2] - cₘ₂[1])


plot(sol2,xlim=(1e-12, 0.06))
plot!(sol3, )
plot!([ cₘ₁[3] cₘ₂[3] cⱼ[3]],
[ cₘ₁[2] cₘ₂[2]  cⱼ[2]],
seriestype=:scatter, labels=["M1" "M2" "J"]
,legend=:outerright)



σ₂ₘ₁ = eigenvectors([cₘ₁[2] cₘ₁[3] cⱼ[4]], ec)[1] 
σ₂ₘ₂ = eigenvectors([cₘ₂[2] cₘ₂[3] cⱼ[4]], ec)[1]
σ₃ₘ₂ = eigenvectors([cₘ₂[2] cₘ₂[3] cⱼ[4]], ec)[2]
σ₃ⱼ =  eigenvectors([cⱼ[2] cⱼ[3] cⱼ[4]], ec)[2] 

if σ₃ₘ₂ >= σ₃ⱼ
    𝒲₃ = :shock
    c₂ₗ, c₃ₗ, σₗ, λₗ  = RH_eigenvalues(cⱼ, cₘ₂, ec)
else
    𝒲₃ = :rarefaction
    c₂ₗ, c₃ₗ, σₗ, λₗ = integral_eigenvalues(cⱼ, cₘ₂, 3, sol3, ec)
end

𝒲₃


if σ₂ₘ₁ >= σ₂ₘ₂
    𝒲₂ = :shock
    c₂ᵣ, c₃ᵣ, σᵣ, λᵣ  = RH_eigenvalues(cₘ₂, cₘ₁, ec)
else
    𝒲₂ = :rarefication
    c₂ᵣ, c₃ᵣ, σᵣ, λᵣ = integral_eigenvalues(cₘ₂, cₘ₁, 2, sol2, ec)
end

𝒲₂

c₄ₗ = cⱼ[4] * ones(length(c₃ₗ))
c₄ᵣ = cⱼ[4] * ones(length(c₃ᵣ))

# * Get ions composition and their wave velocities
c₂ = [cⱼ[2];  c₂ₗ;cₘ₂[2]; c₂ᵣ ;cₘ₁[2]; cₘ₁[2]; cᵢ[2]; cᵢ[2]]
c₃ = [cⱼ[3];  c₃ₗ;cₘ₂[3]; c₃ᵣ ;cₘ₁[3]; cₘ₁[3]; cᵢ[3]; cᵢ[3]]
c₄ = [cⱼ[4];  c₄ₗ;cₘ₂[4]; c₄ᵣ ;cₘ₁[4]; cₘ₁[4]; cᵢ[4]; cᵢ[4]]
λ =  [1e-3;   λₗ; λₗ[end] ;λᵣ  ; λᵣ[end]; 1   ;   1  ;   10 ] 

σ = 1 ./ λ

c₁ = c₄ .- c₃ .- c₂

c = [c₁ c₂ c₃ c₄]

ĉ = zeros(size(c))
for i in 1:size(c)[1]
    ĉ[i,:] = isotherm(c[i,:], ec)
end


t = 0.9
plot(λ * t, c, xlim=(0, 1), layout=4,
		title=["Na" "Mg" "Ca" "Cl"], label=false, lw=2, ylabel="Concentration, M", marker=:circle)

plot(λ * t, ĉ,  xlim=(0, 1),layout=4,
		labels=["Na" "Mg" "Ca" "Cl"], label=false, lw=2, ylabel="Concentration, M", marker=:circle,
        colors=[:brown :darkorchid :orangered])
        plot!(ylims=(1e-4, 1), yscale=:log10)

        plot(λ * t, ĉ, xlim=(0, 1), layout=4,
		title=["COO" "S" "Gly" "Na"], label=false, lw=2, ylabel="Concentration, M", marker=:circle)

λ₁ =  σ[2] - 1
λ₂ = σ[4] - 1 
λ₃ = σ[6] - 1

Λ = 