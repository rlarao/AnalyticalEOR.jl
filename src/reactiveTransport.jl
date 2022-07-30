
function solve_Ion_Transport(ζᵢ, ζⱼ, ec::IonExchangeProblem)
    cᵢ = ζᵢ .* ec.ν
    cⱼ = ζⱼ .* ec.ν

    # * First intermediate point
    ĉᵢ = isotherm(cᵢ, ec)
    ĉₘ₁ = ĉᵢ
    cₘ₁ = flowingConcentrations(ĉₘ₁, cⱼ[4], ec)

    # * Second intermediate point
    cₘ₂, sol2, sol3 = solve_IntegralCurve(cₘ₁, cⱼ, ec)


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
        𝒲₂ = :rarefaction
        c₂ᵣ, c₃ᵣ, σᵣ, λᵣ = integral_eigenvalues(cₘ₂, cₘ₁, 2, sol2, ec)
    end

    c₄ₗ = cⱼ[4] * ones(length(c₃ₗ))
    c₄ᵣ = cⱼ[4] * ones(length(c₃ᵣ))

    # * Get ions composition and their wave velocities
    c₂ = [cⱼ[2];  c₂ₗ;cₘ₂[2]; c₂ᵣ ;cₘ₁[2]; cₘ₁[2]; cᵢ[2]; cᵢ[2]]
    c₃ = [cⱼ[3];  c₃ₗ;cₘ₂[3]; c₃ᵣ ;cₘ₁[3]; cₘ₁[3]; cᵢ[3]; cᵢ[3]]
    c₄ = [cⱼ[4];  c₄ₗ;cₘ₂[4]; c₄ᵣ ;cₘ₁[4]; cₘ₁[4]; cᵢ[4]; cᵢ[4]]
    λ =  [1e-3;   λₗ; λₗ[end] ;λᵣ  ; λᵣ[end]; 1   ;   1  ;   100 ] 
    
    σ = 1 ./ λ
    c₁ = c₄ .- c₃ .- c₂
    
    c = [c₁ c₂ c₃ c₄]
    
    ĉ = zeros(size(c))
    for i in 1:size(c)[1]
        ĉ[i,:] = isotherm(c[i,:], ec)
    end

    return IonExchangeSolution(
            ζᵢ,
            ζⱼ,
            ec.ν,
            ec.K₂₁,
            ec.K₃₁,
            ec.K₂₃,
            ec.Z,
            cᵢ, cₘ₁, cₘ₂, cⱼ,
            c, ĉ , λ, σ,
            𝒲₂, 𝒲₃,
            )
end



function isotherm(c::T, ec::IonExchangeProblem) where {T}
    K₂₁ = ec.K₂₁
    K₃₁ = ec.K₃₁
    Z = ec.Z
    ν = ec.ν

	if size(c)[1] == 4
		c₁, c₂, c₃ = c
	else
		c₂, c₃, c₄ = c
		c₁ = c₄ - c₃ - c₂
	end

    β = 1
    α = 0
    if ν[2] == 1
        β += K₂₁ * c₂ / c₁ 
    elseif ν[2] == 2
        α += K₂₁ * c₂ / c₁ ^ ν[2]
    end

    if ν[3] == 1
        β += K₃₁ * c₃ / c₁ 
    elseif ν[3] == 2
        α += K₃₁ * c₃ / c₁ ^ ν[3]
    end

    ĉ₁ = (-β + sqrt(β^2 + 4 * α * Z)) / (2α)
    
    ĉ₂ = K₂₁ * c₂ * ĉ₁^ν[2] / c₁^ν[2]
    ĉ₃ = K₃₁ * c₃ * ĉ₁^ν[3] / c₁^ν[3]
 
	return [ĉ₁, ĉ₂, ĉ₃, 0]
end


function flowingConcentrations(ĉ, cⱼ₄, ec::IonExchangeProblem)
    K₂₁ = ec.K₂₁
    K₃₁ = ec.K₃₁
    K₂₃ = ec.K₂₃
    ν = ec.ν

    ĉ₁, ĉ₂, ĉ₃ = ĉ
    
    α = 0
    β = 1

    η₂ =  ĉ₂ / K₂₁ / ĉ₁^ ν[2] 
    η₃ =  ĉ₃ / K₃₁ / ĉ₁^ ν[3] 

    if ν[2] == 1
        β += η₂
    elseif ν[2] == 2
        α += η₂
    end

    if ν[3] == 1
        β += η₃
    elseif ν[3] == 2
        α += η₃
    end

    c₁ = (-β + sqrt(β^2 + 4* α * cⱼ₄)) / (2α)

    c₂ = ĉ₂ * c₁^ ν[2] / K₂₁ / ĉ₁ ^ ν[2] 
    c₃ = ĉ₃ * c₁^ ν[3] / K₃₁/ ĉ₁ ^ ν[3] 

    return [c₁, c₂, c₃, cⱼ₄]
end



function eigenvectors(c, ec::IonExchangeProblem)
    ĉ₂₂, ĉ₂₃, ĉ₃₂, ĉ₃₃ = derivative_functions(c, ec)

	σ₂ = 1 + (ĉ₂₂ + ĉ₃₃ - sqrt((ĉ₂₂ - ĉ₃₃)^2 + 4ĉ₂₃ * ĉ₃₂)) / 2
	σ₃ = 1 + (ĉ₂₂ + ĉ₃₃ + sqrt((ĉ₂₂ - ĉ₃₃)^2 + 4ĉ₂₃ * ĉ₃₂)) / 2

	return σ₂, σ₃
end


function derivative_functions(c, ec::IonExchangeProblem)
    ∇(f, x) = derivative(f, x)
	
    ĉ₂(c₂, c₃, c₄) = isotherm([c₂, c₃, c₄], ec)[2]
    ĉ₃(c₂, c₃, c₄) = isotherm([c₂, c₃, c₄], ec)[3]

    ∂ĉ₂∂c₂(c₂, c₃, c₄) = ∇(c₂ -> ĉ₂(c₂, c₃, c₄), c₂)
    ∂ĉ₂∂c₃(c₂, c₃, c₄) = ∇(c₃ -> ĉ₂(c₂, c₃, c₄), c₃)

    ∂ĉ₃∂c₃(c₂, c₃, c₄) = ∇(c₃ -> ĉ₃(c₂, c₃, c₄), c₃)
    ∂ĉ₃∂c₂(c₂, c₃, c₄) = ∇(c₂ -> ĉ₃(c₂, c₃, c₄), c₂)

    return ∂ĉ₂∂c₂(c...), ∂ĉ₂∂c₃(c...), ∂ĉ₃∂c₂(c...), ∂ĉ₃∂c₃(c...)
end


function dc₂dc₃(c, ec::IonExchangeProblem)
    ĉ₂₂, ĉ₂₃, ĉ₃₂, ĉ₃₃ = derivative_functions(c, ec)

    σ₂ = 1 + (ĉ₂₂ + ĉ₃₃ - sqrt((ĉ₂₂ - ĉ₃₃)^2 + 4ĉ₂₃ * ĉ₃₂)) / 2
	σ₃ = 1 + (ĉ₂₂ + ĉ₃₃ + sqrt((ĉ₂₂ - ĉ₃₃)^2 + 4ĉ₂₃ * ĉ₃₂)) / 2

    return [ĉ₂₃ / (σ₂ - 1 - ĉ₂₂), ĉ₂₃ / (σ₃ - 1 - ĉ₂₂)]
end




function integralcurves(u, p, t)
    c₄, ec = p
    c₂, c₃ = u, t
    return dc₂dc₃([c₂, c₃, c₄], ec)
end


function M2_ODE2(c₃ₘ₂, cⱼ, cₘ₁, ec::IonExchangeProblem)
    

    c₃ = 10 .^ range(log10(cₘ₁[3]), log10(cⱼ[3]), length=10000)

    f2(u, p, t) = integralcurves(u, p, t)[1]   

    prob2 = ODEProblem(f2,
                    cₘ₁[2], 			    # u0
                    (cₘ₁[3], c₃ₘ₂), 		# tspan
                    (cⱼ[4], ec), 			# p
                        )

    sol2 = solve(prob2, RadauIIA5()  ,
                                            reltol=1e-12,
                                            abstol=1e-12,
                                            # alg_hints=[:interpolant],
                                            # maxiters=1e7,
                                            # alg_hints=[:stiff],
                                            saveat=c₃ )
    return sol2
end


function M2_ODE3(c₃ₘ₂, cⱼ, cₘ₁, ec::IonExchangeProblem)
    
    c₃ = 10 .^ range(log10(cⱼ[3]), log10(cₘ₁[3]),  length=10000)


    f3(u, p, t) = integralcurves(u, p, t)[2]

    prob3 = ODEProblem(f3, 
                    cⱼ[2],				    # u0
                    (cⱼ[3], c₃ₘ₂), 		    # tspan
                    (cⱼ[4], ec))			# p
                    
    sol3 = solve(prob3, RadauIIA5() ,
                                        reltol=1e-12,
                                        abstol=1e-12,
                                        maxiters=1e3,
                                        # alg_hints=[:interpolant],
                                        saveat=c₃,
                                        )

    return sol3
end


function try_M2_ODE3(c₃₁, cⱼ, cₘ₁, ec)
    try
        sol3 = M2_ODE3(c₃₁, cⱼ, cₘ₁, ec)
        return false
    catch err
        return true
    end
end


function try_M2_ODE2(c₃₂, cⱼ, cₘ₁, ec::IonExchangeProblem)
    try
        sol2 = M2_ODE2(c₃₂, cⱼ, cₘ₁, ec)
        return true
    catch err
        return false
    end
end

function binary_search(fun, c1, c2, c3, ec)
    left = 1
    right = length(c1)

    mid = left + (right - left) ÷ 2

    while left<right
        mid = left + (right - left) ÷ 2

        if fun(c1[mid], c2, c3, ec)
            right = mid
        else
            left = mid + 1
        end
    end

    return left
end


function solve_IntegralCurve(cₘ₁, cⱼ, ec::IonExchangeProblem)
	# c₃ₘ₂ = collect(range(cⱼ[3], cₘ₁[3], length=100000))
	c₃ₘ₂ = collect(10 .^ range(log10(cⱼ[3]), log10(cₘ₁[3]), length=10000))
	
    i = binary_search(try_M2_ODE2, c₃ₘ₂, cⱼ, cₘ₁, ec)
    sol2 = M2_ODE2(c₃ₘ₂[i], cⱼ, cₘ₁, ec)

    i = binary_search(try_M2_ODE3, c₃ₘ₂, cⱼ, cₘ₁, ec)
    sol3 = M2_ODE3(c₃ₘ₂[i-1], cⱼ, cₘ₁, ec)

    c₃ₘ₂ = fzero(c -> sol2(c) - sol3(c), c₃ₘ₂[end])

    c₂ₘ₂ = sol2(c₃ₘ₂)
	cₘ₂ = [c₂ₘ₂, c₃ₘ₂, cⱼ[4]]
	prepend!(cₘ₂, cₘ₂[3] - cₘ₂[2] - cₘ₂[1])

	return cₘ₂, sol2, sol3
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
    c₃ = collect(range(cₗ[3], cᵣ[3], length=20))
    c₂ = [sol(c) for c in c₃]

    if p == 2
        σ = [eigenvectors([c₂ c₃ c₄], ec)[1] for (c₂, c₃) in zip(c₂, c₃)]
    elseif p == 3
        σ = [eigenvectors([c₂ c₃ c₄], ec)[2] for (c₂, c₃) in zip(c₂, c₃)]
    end

    λ = 1 ./ σ

    return c₂, c₃, σ, λ
end