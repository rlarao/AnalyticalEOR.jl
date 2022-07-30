using AnalyticalEOR
using Plots
using Roots
using LaTeXStrings
using DifferentialEquations
theme(:dao; grid_lw=0.2)

ζᵢ = [	0.0017, # Na
0.00124/2, # Mg
0.0062/2, # Ca
]
push!(ζᵢ, ζᵢ[1] + 2ζᵢ[2] + 2ζᵢ[3]) # Anion Cl conc. by charge balance

# * Injected concentrations
ζⱼ = [	0.47, 	# Na
0.098/2, # Mg
0.023/2, 	# Ca
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

ec = ExchangeConstants(K₂₁, K₃₁, K₂₃, Z, ν)

it = solve_Ion_Transport(ζᵢ, ζⱼ, ec)

    cᵢ = ζᵢ .* ec.ν
    cⱼ = ζⱼ .* ec.ν

    # * First intermediate point
    ĉᵢ = isotherm(cᵢ, ec)
    ĉₘ₁ = ĉᵢ
    cₘ₁ = flowingConcentrations(ĉₘ₁, cⱼ[4], ec)

	c₃ₘ₂ = collect(range(cⱼ[3], cₘ₁[3], length=10000))



	c₃ₘ₂ = collect(range(cⱼ[3], cₘ₁[3], length=100000))
	# c₃ₘ₂ = collect(10 .^ range(log10(cⱼ[3]), log10(cₘ₁[3]), length=10000))
	
    i = binary_search(try_M2_ODE2, c₃ₘ₂, cⱼ, cₘ₁, ec)
    sol2 = M2_ODE2(c₃ₘ₂[i], cⱼ, cₘ₁, ec)


    i = binary_search(try_M2_ODE3, c₃ₘ₂, cⱼ, cₘ₁, ec)
    sol3 = M2_ODE3(c₃ₘ₂[i-1], cⱼ, cₘ₁, ec)

    # * Second intermediate point
    cₘ₂, sol2, sol3 = solve_IntegralCurve(cₘ₁, cⱼ, ec)

    c₃ₘ₂ = fzero(c -> sol2(c) - sol3(c), c₃ₘ₂[end])

    cₘ₁
    cₘ₂

    σ₁ = 1
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

    if σ₂ₘ₁ >= σ₂ₘ₂
        𝒲₂ = :shock
        c₂ᵣ, c₃ᵣ, σᵣ, λᵣ  = RH_eigenvalues(cₘ₂, cₘ₁, ec)
    else
        𝒲₂ = :rarefication
        c₂ᵣ, c₃ᵣ, σᵣ, λᵣ = integral_eigenvalues(cₘ₂, cₘ₁, 2, sol2, ec)
    end

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
