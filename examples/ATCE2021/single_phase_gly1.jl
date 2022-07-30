using AnalyticalEOR
using Plots
using Roots
using LaTeXStrings
theme(:dao; grid_lw=0.2)

ν =  [1, 2, 1, 1] # charges

ζᵢ = [	25e-3,   # COO	
1e-3,  # S
1e-9,  # Gly
		]
push!(ζᵢ, ν[1]*ζᵢ[1] + ν[2]*ζᵢ[2] + ν[3]*ζᵢ[3]) # Anion Cl conc. by charge balance

# * Injected concentrations
ζⱼ =  [ 25e-3,   # COO
7e-3,  # S
120e-3,  # Gly
     	]
push!(ζⱼ, ν[1]*ζⱼ[1] + ν[2]*ζⱼ[2] + ν[3]*ζⱼ[3]) # Anion Cl conc. by charge balance

# * Equilibrium constants and 
K₂₁ = 10^1.46 #10^1.67148
K₃₁ = 10^1.14  
K₂₃ = K₂₁ / K₃₁

# * Cation Exchange Capacity
Z = 0.7

ec = ExchangeConstants(K₂₁, K₃₁, K₂₃, Z, ν)

cᵢ = ζᵢ .* ν
cⱼ = ζⱼ .* ν
ĉᵢ = isotherm(cᵢ, ec)
ĉₘ₁ = ĉᵢ
cₘ₁ = flowingConcentrations(ĉₘ₁, cⱼ[4], ec)

c₃ₘ₂ = collect(range(cⱼ[3], cₘ₁[3], length=10))

c₃₂ = c₃ₘ₂[1]
sol2 = M2_ODE2(c₃₂, cⱼ, cₘ₁, ec)

c₃₁ = c₃ₘ₂[end]
sol3 = M2_ODE3(c₃₁, cⱼ, cₘ₁, ec)


c₃ₘ₂ = fzero(c -> sol2(c) - sol3(c), 1e-8)
c₂ₘ₂ = sol2(c₃ₘ₂)
cₘ₂ = [c₂ₘ₂, c₃ₘ₂, cⱼ[4]]
prepend!(cₘ₂, cₘ₂[3] - cₘ₂[2] - cₘ₂[1])

plot([ cₘ₁[3] cₘ₂[3]  cⱼ[3]],
[ cₘ₁[2]  cₘ₂[2]  cⱼ[2]],
seriestype=:scatter, legend=false, ms=5, marker=:circle,
color=[:orangered3  :black :dodgerblue2 ])
plot!(sol2.t, sol2.u, lw=3, alpha=0.6, color=:orangered3)
plot!(sol3.t, sol3.u, lw=3, alpha=0.6, color=:dodgerblue2)
plot!(xlim=(1e-9,1), scale=:log10, size=(450, 420),)
plot!(xlabel=L"\mathrm{Gly \ c_3, eqmol/L}", ylabel=L"\mathrm{SO_4^{2-} \ c_2, eqmol/L}")
plot!(ann=(10^-8, 10^-1.36, "P"))
plot!(ann=(10^-0.5, 10^-1.85, "J"))
plot!(ann=(10^-8.4, 10^-0.83, "Q"))
plot!(ann=(10^-7.75, 10^-1.1, L"\mathcal{W_2}"))
plot!(ann=(10^-0.7, 10^-1, L"\mathcal{W_3}"))
plot!(title="a) Glycine 1 wt%")
savefig("gly_1wt_ODE_sols.svg")

cⱼ

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

ĉ 

# * Plot concentration profiles
t = 1.0

σ2 = σ .- 1 
λ2 = 1 ./ σ2 

λ2[6:7] .= 1.0
λ2[8] = 10

plot(λ * t, c, xlim=(0, 1.5), layout=4,
		title=[L"c_1, \ COO" L"c_2, \ S" L"c_3, \ Gly" L"c_4, \ Na"], label=false, alpha=0.8, lw=3,
        ls=:dash, color=:dodgerblue2)

plot!(λ * t, ĉ, xlim=(0, 1.5), layout=4,
		 label=false, lw=2, alpha=0.8, ylabel="Concentration, eqmol/L",
        size=(800,600), xlabel=L"{x} / {t}", color=:orangered3)

        savefig("conc_profile_gly1.svg")

plot(λ * t, ĉ[:,3], xlim=(0, 1), ylim=(1e-9,1.0), lw=3,
 alpha=0.6, yscale=:log10, legend=false)
plot!(λ * t, c[:,3], lw=2.5, alpha=0.6, ls=:dot)
plot!(xlabel="Dimensionless distance",
    ylabel="Concentration, eqmol/L",
    size=(460, 400), title=L"\mathrm{Glycine} \ c_3" )



    plot(λ * t, ĉ[:,2], xlim=(0, 1), ylim=(1e-3,1.0), lw=3,
    alpha=0.6, yscale=:log10, legend=false)
   plot!(λ * t, c[:,2], lw=2.5, alpha=0.6, ls=:dot)
   plot!(xlabel="Dimensionless distance",
       ylabel="Concentration, eqmol/L",
       size=(460, 400), title=L"SO_4^{2-} \ c_3" )

       plot(λ * t, ĉ[:,1], xlim=(0, 1), ylim=(1e-3,1.0), lw=3,
       alpha=0.6, yscale=:log10, legend=false)
      plot!(λ * t, c[:,1], lw=2.5, alpha=0.6, ls=:dot)
      plot!(xlabel="Dimensionless distance",
          ylabel="Concentration, eqmol/L",
          size=(460, 400), title=L"COO^- \ c_1" )
      
   

cᵢ
cⱼ
cₘ₁
cₘ₂
cⱼ
λ
ĉ[:,1]

ĉ[:,3]

hline!([ĉⱼ[1]], lw=1, ls=:dash)