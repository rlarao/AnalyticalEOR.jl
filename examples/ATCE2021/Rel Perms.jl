using Revise
using DifferentialEquations
using AnalyticalEOR
using Plots
using ForwardDiff:derivative
using Roots

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
K₂₁ = 10^0.6
K₃₁ = 10^0.8     
K₂₃ =  10^-0.16

# * Cation Exchange Capacity
ρ = 2.8082271 # g/cm3
ϕ = 0.30
# cec = 0.06 # 
cec = 1.1e-2
Z = cec * ((1 - ϕ) / ϕ) * ρ # Conversion of cation exchange capacity into moles/liter

ec = ExchangeConstants(K₂₁, K₃₁, K₂₃, Z)

cᵢ = ζᵢ .* ν
cⱼ = ζⱼ .* ν





it = solve_Ion_Transport(ζᵢ, ζⱼ, ν, ec)

c = it.c
ĉ = it.ĉ
λ = it.λ
σ = it.σ
it.W2
it.W3


ζ = c #./ ν'

plot(c[2:end-1,3], λ[2:end-1])



t = 0.5
plot(λ * t, c, xlim=(0, 1), layout=4,
		title=["Na" "Mg" "Ca" "Cl"], label=false, lw=2, ylabel="Concentration, M")

plot(λ * t, ĉ, xlim=(0, 1), layout=4,
		title=["Na" "Mg" "Ca" "Cl"], label=false, lw=2, ylabel="Concentration, M")





  # * First intermediate point
ĉᵢ = isotherm(cᵢ, ec)
ĉₘ₁ = ĉᵢ
cₘ₁ = flowingConcentrations(ĉₘ₁, cⱼ[4], ec)

c₃ₘ₂ = rand(range(cⱼ[3], cₘ₁[3], length=10000))


f2(u, p, t) = integralcurves(u, p, t)[1]

prob2 = ODEProblem(f2,
				cₘ₁[2], 			    # u0
				(cₘ₁[3], cⱼ[3]), 		# tspan
				(cⱼ[4], ec), 			# p
					) 
sol2 = DifferentialEquations.solve(prob2, reltol=1e-12)


f3(u, p, t) = integralcurves(u, p, t)[2]

prob3 = ODEProblem(f3, 
				cⱼ[2],				    # u0
				(cⱼ[3], cₘ₁[3]), 		# tspan
				(cⱼ[4], ec))			# p
sol3 = DifferentialEquations.solve(prob3, reltol=1e-12)



plot(sol3, ylim=(0, 0.07),xlim=(0.06, 0.16))
plot!(sol2, ls=:dash, ylim=(0, 0.07),xlim=(0.06, 0.16))
plot!([cₘ₁[3] cⱼ[3]],
[cₘ₁[2] cⱼ[2]], seriestype=:scatter, labels=["m1" "J"])
plot!(legend=:outerright)


c₃ₘ₂ = fzero(c -> sol2(c) - sol3(c), c₃ₘ₂)
c₂ₘ₂ = sol2(c₃ₘ₂)
cₘ₂ = [c₂ₘ₂, c₃ₘ₂, cⱼ[4]]
prepend!(cₘ₂, cₘ₂[3] - cₘ₂[2] - cₘ₂[1])

plot!([cₘ₂[3]],
[cₘ₂[2]], seriestype=:scatter, labels="m2")

σ₁ = 1
σ₂ₘ₁ = eigenvectors([cₘ₁[2] cₘ₁[3] cⱼ[4]], ec)[1] 
σ₂ₘ₂ = eigenvectors([cₘ₂[2] cₘ₂[3] cⱼ[4]], ec)[1]
σ₃ₘ₂ = eigenvectors([cₘ₂[2] cₘ₂[3] cⱼ[4]], ec)[2]
σ₃ⱼ =  eigenvectors([cⱼ[2] cⱼ[3] cⱼ[4]], ec)[2] 

σ₃ⱼ
σ₃ₘ₂

# if σ₃ₘ₂ > σ₃ⱼ
	𝒲₃ = :shock
	c₂ₗ, c₃ₗ, σₗ, λₗ  = RH_eigenvalues(cⱼ, cₘ₂, ec)
# else
# 	𝒲₃ = :rarefaction
# 	# c₂ₗ, c₃ₗ, σₗ, λₗ = integral_eigenvalues(cⱼ, cₘ₂, 3, sol3, ec)
# end

if σ₂ₘ₁ >= σ₂ₘ₂
	𝒲₂ = :shock
	c₂ᵣ, c₃ᵣ, σᵣ, λᵣ  = RH_eigenvalues(cₘ₂, cₘ₁, ec)
else
	𝒲₂ = :rarefication
	# c₂ᵣ, c₃ᵣ, σᵣ, λᵣ = integral_eigenvalues(cₘ₂, cₘ₂, 2, sol2, ec)
end


c₄ₗ = cⱼ[4] * ones(length(c₃ₗ))
c₄ᵣ = cⱼ[4] * ones(length(c₃ᵣ))

# * Get ions composition and their wave velocities
c₂ = [cⱼ[2]; cⱼ[2];  c₂ₗ; c₂ᵣ[end] ; c₂ᵣ ; cₘ₁[2]; cₘ₁[2]; cᵢ[2]; cᵢ[2]]
c₃ = [cⱼ[3]; cⱼ[3];  c₃ₗ; c₃ᵣ[end] ; c₃ᵣ ; cₘ₁[3]; cₘ₁[3]; cᵢ[3]; cᵢ[3]]
c₄ = [cⱼ[4]; cⱼ[4];  c₄ₗ; c₄ᵣ[end] ; c₄ᵣ ; cₘ₁[4]; cₘ₁[4]; cᵢ[4]; cᵢ[4]]
λ =  [1e-12; λₗ[1];  λₗ;  λₗ[end] ; λᵣ  ; λᵣ[end];  1   ;   1  ;   10 ] 
σ =  [1e12 ; σₗ[1];  σₗ;  σₗ[end] ; σᵣ ; σᵣ[end];  1   ;   1  ;   0.1] 
c₁ = c₄ .- c₃ .- c₂

c = [c₁ c₂ c₃ c₄]


cᵢ


t = 0.3
plot(λ * t, c, xlim=(0, 1), layout=4,
		title=["Na" "Mg" "Ca" "Cl"], label=false, lw=2, ylabel="Concentration, M", marker=:circle)

plot(c[2:end-1, 3], λ[2:end-1], marker=:circle)


c₂ = [cⱼ[2]; cⱼ[2];  c₂ₗ; c₂ᵣ[end] ; c₂ᵣ ; cₘ₁[2]; cₘ₁[2]; cᵢ[2]; cᵢ[2]]
c₃ = [cⱼ[3]; cⱼ[3];  c₃ₗ; c₃ᵣ[end] ; c₃ᵣ ; cₘ₁[3]; cₘ₁[3]; cᵢ[3]; cᵢ[3]]
c₄ = [cⱼ[4]; cⱼ[4];  c₄ₗ; c₄ᵣ[end] ; c₄ᵣ ; cₘ₁[4]; cₘ₁[4]; cᵢ[4]; cᵢ[4]]
λ =  [1e-12; λₗ[1];  λₗ;  λₗ[end] ; λᵣ  ; λᵣ[end];  1   ;   1  ;   10 ] 
σ =  [1e12 ; σₗ[1];  σₗ;  σₗ[end] ; σᵣ ; σᵣ[end];  1   ;   1  ;   0.1] 
c₁ = c₄ .- c₃ .- c₂

c
λ










cₘ₂
cₘ₁
cⱼ


c₂ₗ, c₃ₗ, σₗ, λₗ  = RH_eigenvalues(cⱼ, cₘ₁, ec)

λₗ


it = solve_Ion_Transport(ζᵢ, ζⱼ, ν, ec)

c = it.c
ĉ = it.ĉ
λ = it.λ
σ = it.σ


plot(λ * t, c[:,3], marker=:circle, xlim=(0,1))

plot(σ, c[:,2], yscale=:log10, ylim=(1e-5, 1), xlim=(0, 20),
		title=["Mg"], label=false, lw=2, ylabel="Concentration, M", marker=:circle)

# plot!(σ, c[:,2], yscale=:log10, ylim=(1e-5, 1), xlim=(0, 20),
# 		title=["Mg"], label=false, lw=2, ylabel="Concentration, M", marker=:circle)

c[:,2]
λ

plot(c[:,2], λ)

ĉ = zeros(size(c))




for i in 1:size(c)[1]
	ĉ[i,:] = isotherm(c[i,:], ec)
end

isotherm(c[4,:], ec)

it.W2
it.W3

		cₘ₁ = it.cₘ₁
		cₘ₂ = it.cₘ₂

	c₂ = [cⱼ[2]; cⱼ[2];  cₘ₂[2] ; cₘ₁[2]; cₘ₁[2]; cₘ₁[2]; cᵢ[2]; cᵢ[2]]
	c₃ = [cⱼ[3]; cⱼ[3];  cₘ₂[3] ; cₘ₁[3]; cₘ₁[3]; cₘ₁[3]; cᵢ[3]; cᵢ[3]]
	c₄ = [cⱼ[4]; cⱼ[4];  cₘ₂[4] ; cₘ₁[4]; cₘ₁[4]; cₘ₁[4]; cᵢ[4]; cᵢ[4]]
	λ2 =  [1e-3; 0.6236; 0.63236; 0.63236   ; 0.9862;  1   ;   1  ;   10 ] 
	c₁ = c₄ .- c₃ .- c₂

plot(λ2 * t, c₃, marker=:circle, xlim=(0,1))
# plot(λ2 * t, c₂)

cᵢ[3]
cₘ₁[3]
cₘ₂[3]
cⱼ[3]


function RH_eigenvalues(cₗ, cᵣ, ec)
	ĉ₁ᵣ, ĉ₂ᵣ, ĉ₃ᵣ, ĉ₄ᵣ = isotherm(cᵣ, ec)
	ĉ₁ₗ, ĉ₂ₗ, ĉ₃ₗ, ĉ₄ₗ = isotherm(cₗ, ec)
	σ̃ = ((ĉ₃ᵣ + cᵣ[3]) - (ĉ₃ₗ + cₗ[3])) / (cᵣ[3] - cₗ[3])
	λ̃ = 1 ./ σ̃ 
	return cₗ[2], cₗ[3], σ̃, λ̃
end


function derivative_functions(c, ec::ExchangeConstants)
    ∇(f, x) = derivative(f, x)
	
    ĉ₂(c₂, c₃, c₄) = isotherm([c₂, c₃, c₄], ec)[2]
    ĉ₃(c₂, c₃, c₄) = isotherm([c₂, c₃, c₄], ec)[3]

    ∂ĉ₂∂c₂(c₂, c₃, c₄) = ∇(c₂ -> ĉ₂(c₂, c₃, c₄), c₂)
    ∂ĉ₂∂c₃(c₂, c₃, c₄) = ∇(c₃ -> ĉ₂(c₂, c₃, c₄), c₃)

    ∂ĉ₃∂c₃(c₂, c₃, c₄) = ∇(c₃ -> ĉ₃(c₂, c₃, c₄), c₃)
    ∂ĉ₃∂c₂(c₂, c₃, c₄) = ∇(c₂ -> ĉ₃(c₂, c₃, c₄), c₂)

    return ∂ĉ₂∂c₂(c...), ∂ĉ₂∂c₃(c...), ∂ĉ₃∂c₂(c...), ∂ĉ₃∂c₃(c...)
end


function dc₂dc₃(c, ec::ExchangeConstants)
    ĉ₂₂, ĉ₂₃, ĉ₃₂, ĉ₃₃ = derivative_functions(c, ec::ExchangeConstants)

    σ₂ = 1 + (ĉ₂₂ + ĉ₃₃ - sqrt((ĉ₂₂ - ĉ₃₃)^2 + 4ĉ₂₃ * ĉ₃₂)) / 2
    σ₃ = 1 + (ĉ₂₂ + ĉ₃₃ + sqrt((ĉ₂₂ - ĉ₃₃)^2 + 4ĉ₂₃ * ĉ₃₂)) / 2
    
    return [ĉ₂₃ / (σ₂ - 1 - ĉ₂₂), ĉ₂₃ / (σ₃ - 1 - ĉ₂₂)]
end


function integralcurves(u, p, t)
    c₄, ec = p
    c₂, c₃ = u, t
    return dc₂dc₃([c₂, c₃, c₄], ec)
end


function eigenvectors(c, ec::ExchangeConstants)
    ĉ₂₂, ĉ₂₃, ĉ₃₂, ĉ₃₃ = derivative_functions(c, ec)

	σ₂ = 1 + (ĉ₂₂ + ĉ₃₃ - sqrt((ĉ₂₂ - ĉ₃₃)^2 + 4ĉ₂₃ * ĉ₃₂)) / 2
	σ₃ = 1 + (ĉ₂₂ + ĉ₃₃ + sqrt((ĉ₂₂ - ĉ₃₃)^2 + 4ĉ₂₃ * ĉ₃₂)) / 2

	return σ₂, σ₃
end


function RH_eigenvalues(cₗ, cᵣ, ec)
    ĉ₁ᵣ, ĉ₂ᵣ, ĉ₃ᵣ, ĉ₄ᵣ = isotherm(cᵣ, ec)
    ĉ₁ₗ, ĉ₂ₗ, ĉ₃ₗ, ĉ₄ₗ = isotherm(cₗ, ec)
    σ̃ = ((ĉ₃ᵣ + cᵣ[3]) - (ĉ₃ₗ + cₗ[3])) / (cᵣ[3] - cₗ[3])
    λ̃ = 1 ./ σ̃ 
    return cₗ[2], cₗ[3], σ̃, λ̃
end



function integral_eigenvalues(cₗ, cᵣ, p, sol, ec)
    c₄ = cₗ[4]
    c₃ = collect(range(cₗ[3], cᵣ[3], length=50))
    c₂ = [sol(c) for c in c₃]

    if p == 2
        σ = [eigenvectors([c₂ c₃ c₄], ec)[1] for (c₂, c₃) in zip(c₂, c₃)]
    elseif p == 3
        σ = [eigenvectors([c₂ c₃ c₄], ec)[2] for (c₂, c₃) in zip(c₂, c₃)]
    end

    λ = 1 ./ σ

    return c₂, c₃, σ, λ
end