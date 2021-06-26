
function solve_Ion_Transport(ζᵢ, ζⱼ, ν, ec::ExchangeConstants)
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
    σ₂ₘ₁ = eigenvectors([cₘ₂[2] cₘ₂[3] cⱼ[4]], ec)[1]
    σ₂ₘ₂ = eigenvectors([cₘ₁[2] cₘ₁[3] cⱼ[4]], ec)[1] 
    σ₃ⱼ =  eigenvectors([cⱼ[2] cⱼ[3] cⱼ[4]], ec)[2] 
    σ₃ₘ₂ = eigenvectors([cₘ₂[2] cₘ₂[3] cⱼ[4]], ec)[2]


    if σ₃ₘ₂ >= σ₃ⱼ
        𝒲₃ = :shock
        c₂ᵣ, c₃ᵣ, σᵣ, λᵣ  = RH_eigenvalues(cⱼ, cₘ₂, ec)
        c₂ᵣ, c₃ᵣ, σᵣ, λᵣ = integral_eigenvalues(cⱼ, cₘ₂, 2, sol2, ec)
    else
        𝒲₃ = :rarefication
    end

    if σ₂ₘ₁ >= σ₂ₘ₂
        𝒲₂ = :shock
        c₂ₗ, c₃ₗ, σₗ, λₗ  = RH_eigenvalues(cₘ₂, cₘ₁, ec)
    else
        𝒲₂ = :rarefication
        c₂ₗ, c₃ₗ, σₗ, λₗ = integral_eigenvalues(cⱼ, cₘ₂, 3, sol3, ec)
    end


    # Build complete solution
    c₄ₗ = cⱼ[4] * ones(length(c₃ₗ))
    c₄ᵣ = cⱼ[4] * ones(length(c₃ᵣ))

    # * Get ions composition and their wave velocities
    c₂ = [cⱼ[2]; cⱼ[2]; c₂ₗ ; c₂ᵣ ; cₘ₁[2]; cₘ₁[2]; cᵢ[2]; cᵢ[2]]
    c₃ = [cⱼ[3]; cⱼ[3]; c₃ₗ ; c₃ᵣ ; cₘ₁[3]; cₘ₁[3]; cᵢ[3]; cᵢ[3]]
    c₄ = [cⱼ[4]; cⱼ[4]; c₄ₗ ; c₄ᵣ ; cₘ₁[4]; cₘ₁[4]; cᵢ[4]; cᵢ[4]]
    λ =  [ 0  ; λₗ[1];  λₗ ; λᵣ  ;λᵣ[end] ;  1   ;   1  ;   10 ] 
    c₁ = c₄ .- c₃ .- c₂

    c = [c₁ c₂ c₃ c₄]
    # IonExchangeTransport(ζᵢ,
    #         ζⱼ,
    #         ν,
    #         ec.K₂₁,
    #         ec.K₃₁,
    #         ec.K₂₃,
    #         ec.Z,
    #         cᵢ, cₘ₁, cₘ₂, cⱼ,
    #         λ₁, λ₂, λ₃,
    #         c₁, c₂, c₃, c₄,
    #         λ
    #                     )
    return c, λ
end



function isotherm(c::T, ec::ExchangeConstants) where {T}
    K₂₁ = ec.K₂₁
    K₃₁ = ec.K₃₁
    Z = ec.Z

	if size(c)[1] == 4
		c₁, c₂, c₃ = c

	else
		c₂, c₃, c₄ = c
		c₁ = c₄ - c₃ - c₂
	end

		ĉ₁ = (-1 + sqrt(1 + (4Z * (K₂₁ * c₂ + K₃₁ * c₃) / c₁^2))
				) / ( 2((K₂₁ * c₂ + K₃₁ * c₃) / c₁^2) )
		
		ĉ₂ = K₂₁ * c₂ * ĉ₁^2 / c₁^2
		ĉ₃ = K₃₁ * c₃ * ĉ₁^2 / c₁^2

	return [ĉ₁, ĉ₂, ĉ₃, 0]
end

function flowingConcentrations(ĉ, cⱼ₄, ec::ExchangeConstants)
    K₂₁ = ec.K₂₁
    K₂₃ = ec.K₂₃

    ĉ₁, ĉ₂, ĉ₃ = ĉ
    
    a = (1 + (ĉ₃ / ĉ₂ * K₂₃))^2
    b = -(2 * cⱼ₄ * (1 + (ĉ₃ / ĉ₂ * K₂₃)) + (K₂₁ * ĉ₁^2 / ĉ₂))
    c = cⱼ₄^2

    c₂ = (-b - sqrt(b^2 - (4a * c))) / (2a)
    c₁ = sqrt(K₂₁ * c₂ / ĉ₂) * ĉ₁
    c₃ = c₂ * K₂₃ * ĉ₃ / ĉ₂

    return [c₁, c₂, c₃, cⱼ₄]
end

function eigenvectors(c, ec::ExchangeConstants)
    ĉ₂₂, ĉ₂₃, ĉ₃₂, ĉ₃₃ = derivative_functions(c, ec)

	σ₂ = 1 + (ĉ₂₂ + ĉ₃₃ - sqrt((ĉ₂₂ - ĉ₃₃)^2 + 4ĉ₂₃ * ĉ₃₂)) / 2
	σ₃ = 1 + (ĉ₂₂ + ĉ₃₃ + sqrt((ĉ₂₂ - ĉ₃₃)^2 + 4ĉ₂₃ * ĉ₃₂)) / 2

	return σ₂, σ₃
end

function derivative_functions(c, ec::ExchangeConstants)
    ∇(f, x) = ForwardDiff.derivative(f, x)
	
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



function M2_ODE_solutions(c₃ₘ₂, cⱼ, cₘ₁, ec::ExchangeConstants)
    
    f2(u, p, t) = integralcurves(u, p, t)[1]   
    f3(u, p, t) = integralcurves(u, p, t)[2]

    prob2 = ODEProblem(f2,
                    cₘ₁[2], 				# u0
                    (cₘ₁[3], cⱼ[3]), 		# tspan
                    (cⱼ[4], ec), 					# p
                        ) 
    sol2 = DifferentialEquations.solve(prob2, BS3(), reltol=1e-12)

    prob3 = ODEProblem(f3, 
                    cⱼ[2],				# u0
                    (cⱼ[3], cₘ₁[3]), 		# tspan
                    (cⱼ[4], ec))				# p
    sol3 = DifferentialEquations.solve(prob3, BS3(), reltol=1e-12)

    return sol2, sol3
end

# function M2_minimize_loss(c₃ₘ₂, cⱼ, cₘ₁, sol2, sol3)
# 	loss_fun(x) = abs(sol2(x) - sol3(x))

# 	model = Model(NLopt.Optimizer)
#     # register(model, :loss_fun, 1, abs(sol2(x) - sol3(x)); autodiff=true)
# 	set_optimizer_attribute(model, "algorithm", :LD_MMA)
# 	@variable(model, x)
# 	@NLobjective(model, Min, loss_fun(x))
# 	@NLconstraint(model, x <= cₘ₁[3])
# 	@NLconstraint(model, x >= cⱼ[3])	

# 	set_start_value(x, c₃ₘ₂)
# 	JuMP.optimize!(model)

# 	return objective_value(model), value(x)
# end

# function solve_IntegralCurve(cₘ₁, cⱼ, ec::ExchangeConstants)
# 	i = 0
# 	loss = 1
# 	c₃ₘ₂ = rand(range(cⱼ[3], cₘ₁[3], length=10000))
	
	
# 	while loss >= 1e-18 && i <= 1000	
# 		try
# 			sol2, sol3 = M2_ODE_solutions(c₃ₘ₂, cⱼ, cₘ₁, ec)
# 			loss, c₃ₘ₂ = M2_minimize_loss(c₃ₘ₂, cⱼ, cₘ₁, sol2, sol3)
# 		catch err
# 			if isa(err, DomainError)
# 				c₃ₘ₂ = rand(range(cⱼ[3], cₘ₁[3], length=10000))
# 				continue
# 			end
# 		end
# 		i += 1
# 	end
	
# 	sol2, sol3 = M2_ODE_solutions(c₃ₘ₂, cⱼ, cₘ₁, ec)
# 	cₘ₂ = [sol3(c₃ₘ₂), c₃ₘ₂, cⱼ[4]]
# 	prepend!(cₘ₂, cₘ₂[3] - cₘ₂[2] - cₘ₂[1])
# 	return cₘ₂
# end


function solve_IntegralCurve(cₘ₁, cⱼ, ec::ExchangeConstants)
	i = 0
	loss = 1
	c₃ₘ₂ = rand(range(cⱼ[3], cₘ₁[3], length=10000))
	
	
    sol2, sol3 = M2_ODE_solutions(c₃ₘ₂, cⱼ, cₘ₁, ec)
    
    c₃ₘ₂ = fzero(c -> sol2(c) - sol3(c), (cⱼ[3] + cₘ₁[3]) / 2)
    c₂ₘ₂ = sol2(c₃ₘ₂)

	cₘ₂ = [c₂ₘ₂, c₃ₘ₂, cⱼ[4]]
	prepend!(cₘ₂, cₘ₂[3] - cₘ₂[2] - cₘ₂[1])
	return cₘ₂
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