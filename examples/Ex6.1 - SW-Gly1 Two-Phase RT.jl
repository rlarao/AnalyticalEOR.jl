using AnalyticalEOR, Plots, LaTeXStrings

theme(:dao; grid_lw=0.2)


begin
    c_colors = ["brown" "darkorchid" "orangered" "gold"]
    c_labels = [L"COO^-" L"SO_4^{2-}" L"Gly^-" L"Na^+"]
    a_colors = ["brown" "darkorchid" "orangered"]
    a_labels = [L"COOX" L"SO_4X_2" L"GlyX"]
end

begin
	ν = [1, 2, 1, 1] # charges

    ζᵢ = [9e-3,   # COO	
        7e-3, # S
        1e-9,  # Gly
    ]
	
	# Cation Na conc. by charge balance
    push!(ζᵢ, ν[1] * ζᵢ[1] + ν[2] * ζᵢ[2] + ν[3] * ζᵢ[3])

    # * Injected concentrations
    ζⱼ = [9e-3,   # COO
        1.961e-2,  # S
        0.1, # Gly
    ]
	
	# Cation Na conc. by charge balance
    push!(ζⱼ, ν[1] * ζⱼ[1] + ν[2] * ζⱼ[2] + ν[3] * ζⱼ[3])
    
    # * Cation Exchange Capacity
    Z = 0.7

    # * Equilibrium constants and 
    K₂₁ = 10^2.46 
    K₃₁ = 10^3.30
    K₂₃ = K₂₁ / K₃₁ 

	ec = IonExchangeProblem(K₂₁, K₃₁, K₂₃, Z, ν)

    it = solve_Ion_Transport(ζᵢ, ζⱼ, ec)
end

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

it.cᵢ * 1000
it.ĉ
it.cₘ₁ * 1000
it.cₘ₂ .* 1000
it.cⱼ .* 1000

cᵢ = ζᵢ .* ν
cⱼ = ζⱼ .* ν
ĉᵢ = isotherm(cᵢ, ec)
ĉₘ₁ = ĉᵢ
cₘ₁ = flowingConcentrations(ĉₘ₁, cⱼ[4], ec)


cₘ₂, sol2, sol3 = solve_IntegralCurve(cₘ₁, cⱼ, ec)

plot(sol2.t, sol2.u, lw=3, alpha=0.6, color=:orangered3)
plot!(sol3.t, sol3.u, lw=3, alpha=0.6, color=:dodgerblue2)
plot!([ cₘ₁[3] cₘ₂[3]  cⱼ[3]],
[ cₘ₁[2]  cₘ₂[2]  cⱼ[2]],
seriestype=:scatter, legend=false, ms=5, marker=:circle,
color=[:orangered3  :black :dodgerblue2 ])
plot!(xlim=(1e-10,1), xscale=:log10, size=(450, 420),)

plot!(xlabel=L"\mathrm{Gly \ c_2, eqmol/L}", ylabel=L"\mathrm{SO_4^{2-} \ c_1, eqmol/L}")
savefig("")
it.c

plot!(ann=(10^-8.2, 10^-0.91, "P"))
plot!(ann=(10^-0.7, 0.04, "J"))
plot!(ann=(10^-8.2, 10^-0.85, "Q"))
plot!(ann=(10^-9, 0.129, L"\mathcal{W_2}"))
plot!(ann=(10e-2, 0.1, L"\mathcal{W_3}"))

it.c .* 100
it.σ




begin
    kr_ww = RelPerms(swr=0.11,
                    sor=0.10,
                    krw0=0.2,
                    kro0=0.8,
                    nw=3.,
                    no=2.)

    kr_ow = RelPerms(swr=0.11,
                    sor=0.16,
                    krw0=0.4,
                    kro0=0.5,
                    nw=2.,
                    no=3.)

	swr(ω)  =    (1-ω) * kr_ww.swr   + ω * kr_ow.swr
	sor(ω)  =    (1-ω) * kr_ww.sor   + ω * kr_ow.sor
	krw0(ω) =    (1-ω) * kr_ww.krw0  + ω * kr_ow.krw0
	kro0(ω) =    (1-ω) * kr_ww.kro0  + ω * kr_ow.kro0
	no(ω)   =    (1-ω) * kr_ww.no    + ω * kr_ow.no
	nw(ω)   =    (1-ω) * kr_ww.nw    + ω * kr_ow.nw
end

# begin
	ĉⱼ = it.ĉ[2,:]
	ĉₘ₂ = it.ĉ[3,:]
	ĉᵢ = it.ĉ[7,:]

	σₘ₂ = it.σ[3]
	σₘ₁ = it.σ[5]
    it.σ
    σₘ₁
    σₘ₂
	σ₁ = σₘ₁  - 1
	σ₂ = σₘ₂  - 1

	coo_ww = 0
	coo_ow = ĉᵢ[1]
	coo_mw1 = ĉₘ₂[1]
	coo_mw2 = ĉⱼ[1]

	ω₁ = (coo_mw1 - coo_ww) / (coo_ow - coo_ww)
	ω₂ = (coo_mw2 - coo_ww) / (coo_ow - coo_ww)

    println((ω₁, ω₂))

    kr_mw1 = RelPerms(swr=swr(ω₁),
                    sor=sor(ω₁),
                    krw0=krw0(ω₁),
                    kro0=kro0(ω₁),
                    nw=nw(ω₁),
                    no=no(ω₁))

    kr_mw2 = RelPerms(swr=swr(ω₂),
                    sor=sor(ω₂),
                    krw0=krw0(ω₂),
                    kro0=kro0(ω₂),
                    nw=nw(ω₂),
                    no=no(ω₂))


   	μw = 1e-3
    μo = 5e-3

    wf = WaterFlooding( kr=kr_ow,
                        μw=μw,
                        μo=μo,
                        )

    cf = ChemicalFlooding(
                    wf = wf,
                    krs = [kr_mw2, kr_mw1],
                    D = [σ₂, σ₁]
                    )

    sol = solve_cf(cf, 0.11, 0.9)
end 

plot_fractional_flow(cf, sol)

kr_mw1

kr_mw2

kr_ow


wf = cf.wf

begin
s = collect(range(wf.kr.swr, 1-wf.kr.sor; length = 101))
plot(s, wf.f(s), xlims=(0,1), lw=2, alpha=0.7, label=false, color=:dodgerblue2)

s = collect(range(cf.kr[1].swr, 1-cf.kr[1].sor; length = 101))
display(plot!(s, cf.f[1](s), label=false, lw=2, alpha=0.7, color=:orangered))


s = collect(range(cf.kr[2].swr, 1-cf.kr[2].sor; length = 101))
display(plot!(s, cf.f[2](s), label=false, lw=2, alpha=0.7, color=:darkorchid))
end


xlabel!("Dimensionless Distance, x")
ylabel!("Water Fractional Flow, f")

D = copy(cf.D)

i = 1
for (s1, s2, wave, speed, f) in sol[1:2]
        if wave == :shock
            display(plot!([s1, s2], [f(s1), f(s2)], lw=2, color=:gray, alpha=0.6, label=false))
            if ~isempty(D)
                display(plot!([-popfirst!(D), s1], [0, f(s1)], lw=2, color=:black, alpha=0.5, ls=:dot, label=false))
            end
        elseif wave == :spreading
            s = collect(range(s1, s2; length=101))
            display(plot!(s, f(s), alpha=0.7, lw=2, color=:gray, label=false))
        end
        display(scatter!([s1, s2], [f(s1), f(s2)], color=:gray, alpha=0.5, label=false))
end


plot!(ann=(0.2, 0.96, L"\mathcal{W_3}"), annotationfontsize=12)


plot!([0.45, 0.5], [cf.f[1](0.45), cf.f[1](0.45)], color=:black, lw=0.5, label=false)
plot!(ann=(0.52, 0.26, L"\textbf{c}^Q, \ \omega \equal 0.008"), annotationfontsize=10, annotationhalign=:left)

plot!([0.32, 0.5], [cf.f[2](0.32), cf.f[2](0.32)], color=:black, lw=0.5, label=false)
plot!(ann=(0.52, 0.17, L"\textbf{c}^P, \ \omega \equal 0.447"), annotationfontsize=10, annotationhalign=:left)


plot!([0.2, 0.5], [wf.f(0.2), wf.f(0.2)], color=:black, lw=0.5, label=false)
plot!(ann=(0.52, 0.1, L"\textbf{c}^I, \ \omega \equal 1.0"), annotationfontsize=10, annotationhalign=:left)
end


plot!(ann=(0.895, 0.963, "J"))
plot!(ann=(0.855, 0.96, "1"))
plot!(ann=(0.7, 0.94, "2"))


v_list

for (s1, s2, wave, speed, f) in sol[3:4]
        if wave == :shock
            display(plot!([s1, s2], [f(s1), f(s2)], lw=2, color=:gray, alpha=0.6, label=false))
            if ~isempty(D)
                display(plot!([-popfirst!(D), s1], [0, f(s1)], lw=2, color=:black, alpha=0.5, ls=:dot, label=false))
            end
        elseif wave == :spreading
            s = collect(range(s1, s2; length=101))
            display(plot!(s, f(s), alpha=0.7, lw=2, color=:gray, label=false))
        end
        display(scatter!([s1, s2], [f(s1), f(s2)], color=:gray, alpha=0.5, label=false))
end



plot!(ann=(0.17, 0.48, L"\mathcal{W_2}"), annotationfontsize=12)
plot!(ann=(0.6, 0.86, "3"))
plot!(ann=(0.35, 0.56, "4"))



for (s1, s2, wave, speed, f) in sol[4:5]
        if wave == :shock
            display(plot!([s1, s2], [f(s1), f(s2)], lw=2, color=:gray, alpha=0.6, label=false))
            if ~isempty(D)
                display(plot!([-popfirst!(D), s1], [0, f(s1)], lw=2, color=:black, alpha=0.5, ls=:dot, label=false))
            end
        elseif wave == :spreading
            s = collect(range(s1, s2; length=101))
            display(plot!(s, f(s), alpha=0.7, lw=2, color=:gray, label=false))
        end
        display(scatter!([s1, s2], [f(s1), f(s2)], color=:gray, alpha=0.5, label=false))
end

s1, s2, wave, speed, f = sol[5]

plot!(ann=(0., 0.12, L"\mathcal{W_1}"), annotationfontsize=12)
plot!([0,s1], [0, wf.f(s1)], lw=2, color=:black, alpha=0.5, ls=:dot, label=false)

tracer = wf.f(s1) / s1

plot!(ann=(0.09, 0.03, "I"))

plot_chemical_waves(sol, 0.5)



s, v = get_saturation_speeds(sol)
plot(v .* t, s, lims=(0,1), lw=3, fill=(0, 0.5), label=false, color=color)
xlabel!("Dimensionless Distance, x")
ylabel!("Saturation, s")




    t  = 0.3
    s, v = get_saturation_speeds(sol)
    begin
    plot(v .* t, s, lims=(0,1), lw=2, label=false, color=:dodgerblue2, fill=(0,0.1), grid=false)
    hline!(s_list, ls=:dot, lw=1, color=:gray, label=false)
    vline!(v_list .* t, ls=:dot, lw=1, color=:gray, label=false)
    xlabel!("Dimensionless distance, x")
    ylabel!("Water saturation, s")

    plot!(ann=(0.08, 0.22, L"\mathcal{W_3}"), annotationfontsize=12)
    plot!(ann=(0.4, 0.22, L"\mathcal{W_2}"), annotationfontsize=12)
    plot!(ann=(0.55, 0.22, L"\mathcal{W_1}"), annotationfontsize=12)


    plot!(ann=(0.85, 0.14, L"{s_I}"), annotationfontsize=12)
    plot!(ann=(0.85, 0.38, L"{s_4}"), annotationfontsize=12)
    plot!(ann=(0.85, 0.63, L"{s_3}"), annotationfontsize=12)
    plot!(ann=(0.85, 0.73, L"{s_2}"), annotationfontsize=12)
    plot!(ann=(0.85, 0.88, L"{s_1}"), annotationfontsize=12)
    plot!(ann=(0.85, 0.93, L"{s_J}"), annotationfontsize=12)
    
    plot!(ann=(0.18, 0.95, L"t=0.35  PVI"), annotationfontsize=12)

    end


    pop!(v_list)
    push!(v_list, tracer)

    v_list = []
    s_list = []
    for (s1, s2, wave, speed, f) in sol
        if s1 ∉ s_list
            push!(s_list, s1)
        end

        if s2 ∉ s_list
            push!(s_list, s2)
        end

        if wave ==:shock 
            push!(v_list, speed)
        end

    end
    s_list

    

    colors = [:orangered, :darkorchid, :gold]

    s, v = get_saturation_speeds(sol)

    x(c_pos) = min.(v .*t, c_pos)

    shock_sols = [s for s in sol if s[3] == :shock][1:end-1]

    s2 = shock_sols[end-1][2] 
    f2 = shock_sols[end-1][5](s2)

    tracer_v = f2 / s2
    shock_sols = vcat(shock_sols, (s2, s2, :shock, tracer_v, f2))
    


    s_old = 1
    i = 1
    for (s1, s2, wave, speed, f) in shock_sols

        c_x = x(speed .* t)
        c_pos = c_x[end]
        s_i = findfirst(c_x .== c_pos) + 1
        display(plot!(c_x[s_old:s_i], s[s_old:s_i], fill=(0, 0.2, colors[i]), lw=0., label=false)) 

        # display(plot!([c_x[s_i], c_x[s_i]], [0, s[s_i]], lw=2, ls=:dash, alpha=0.5, color=colors[i], label=false))

        s_old = s_i
        i += 1
    end
end






 
 begin
	ĉᵢ = it.ĉ[7,:]
	ĉⱼ = it.ĉ[2,:]
	ĉₘ₁ = it.ĉ[6,:]
	ĉₘ₂ = it.ĉ[3,:]


	cᵢ = it.c[7,:]
	cⱼ = it.c[2,:]
	cₘ₁ = it.c[6,:]
	cₘ₂ = it.c[3,:]
 end

gr(size=(500,400), legend=false)
t = 0.3
begin
    p1 =  plot(it.λ * t, it.ĉ[:,1], xlim=(0,1), alpha=0.8, lw=3, color=a_colors[1], title=a_labels[1])
    hline!([ĉᵢ[1], ĉⱼ[1], ĉₘ₁[1], ĉₘ₂[1]], alpha=0.8, ls=:dot, color=:gray, lw=3)
    ylabel!("Concentration eqmol/L")

    p2 =  plot(it.λ * t, it.ĉ[:,2], xlim=(0,1), alpha=0.8, lw=3, color=a_colors[2], title=a_labels[2])
    hline!([ĉᵢ[2], ĉⱼ[2], ĉₘ₁[2], ĉₘ₂[2]], alpha=0.8, ls=:dot, color=:gray, lw=3)
    
    p3 =  plot(it.λ * t, it.ĉ[:,3], xlim=(0,1), alpha=0.8, lw=3, color=a_colors[3], title=a_labels[3])
    hline!([ĉᵢ[3], ĉⱼ[3], ĉₘ₁[3], ĉₘ₂[3]], alpha=0.8, ls=:dot, color=:gray, lw=3)

    p4 =  plot(it.λ * t, it.c[:,1], xlim=(0,1), alpha=0.8, lw=3, color=a_colors[1], title=c_labels[1])
    hline!([cᵢ[1], cⱼ[1], cₘ₁[1], cₘ₂[1]], alpha=0.8, ls=:dot, color=:gray, lw=3)
    xlabel!("Dimensionless distance")
    ylabel!("Concentration eqmol/L")

    p5 =  plot(it.λ * t, it.c[:,2], xlim=(0,1), alpha=0.8, lw=3, color=a_colors[2], title=c_labels[2])
    hline!([cᵢ[2], cⱼ[2], cₘ₁[2], cₘ₂[2]], alpha=0.8, ls=:dot, color=:gray, lw=3)
    xlabel!("Dimensionless distance")

    p6 =  plot(it.λ * t, it.c[:,3], xlim=(0,1), alpha=0.8, lw=3, color=a_colors[3], title=c_labels[3])
    hline!([cᵢ[3], cⱼ[3], cₘ₁[3], cₘ₂[3]], alpha=0.8, ls=:dot, color=:gray, lw=3)

    xlabel!("Dimensionless distance")

    plot(p1, p2, p3, p4, p5, p6,  layout=(2,3))
end