using AnalyticalEOR, Roots, Revise, Plots

begin
	ν  =  [1, 2, 1, 1] 	# charges

	ζᵢ = [	25e-3,   	# COO	
			1e-3,  		# S
			1e-9,  		# Gly
			]
	push!(ζᵢ, ν[1]*ζᵢ[1] + ν[2]*ζᵢ[2] + ν[3]*ζᵢ[3]) # Anion Cl conc. by charge balance
	
	# * Injected concentrations
	ζⱼ =  [ 25e-3,   	# COO
			7e-3,  		# S
			120e-3,  	# Gly
	     	]
	push!(ζⱼ, ν[1]*ζⱼ[1] + ν[2]*ζⱼ[2] + ν[3]*ζⱼ[3]) # Anion Cl conc. by charge balance
	
	# * Equilibrium constants and 
	K₂₁ = 10^1.46 #10^1.67148
	K₃₁ = 10^1.14  
	K₂₃ = K₂₁ / K₃₁
	
	# * Cation Exchange Capacity
	Z = 0.7
	
	ec = ExchangeConstants(K₂₁, K₃₁, K₂₃, Z, ν)
end


it = solve_Ion_Transport(ζᵢ, ζⱼ, ec)

c_labels = ["COO" "SO4" "Gly" "Na"]
a_labels = ["XCOO" "X2SO4" "XGly"]
c_colors = [:brown :darkorchid :orangered :gold]
a_colors = [:brown :darkorchid :orangered]

plot(it.λ * t, it.c, xlim=(0, 1.5), label=c_labels, alpha=0.6, lw=3,
        color=c_colors, marker="o")

plot_flowing_conc(it,0.5, c_labels, c_colors)
plot_adsorbed_conc(it,0.5, a_labels, a_colors)

function plot_flowing_conc(it, t, labels, colors)
    plot(it.λ * t, it.c, xlim=(0, 1), label=labels, alpha=0.6, lw=3,
            color=colors, markershape=:circle)
    xlabel!("Dimensionless distance, x")
    ylabel!("Flowing concentrations eqmol/L")
end

function plot_c_velocities(it)
    plot(it.λ[begin+1:end-1], markershape=:circle)
end