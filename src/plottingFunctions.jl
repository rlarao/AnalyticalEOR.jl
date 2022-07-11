
function plot_fw(wf::WaterFlooding)
    fw(sw) = fractional_flow(sw, wf.kr, wf.μw, wf.μo)
    swi = wf.Si
    swj = wf.Sj
    s̃w = wf.S̃
    swr = wf.kr.swr
    sor = wf.kr.sor
    sw = collect(swr:0.01:1 - sor)

    plot(sw, fw(sw), lw=3, label=false)
    plot!(lims=(0, 1), title="Water Fractional Flow")
    plot!(xlabel="Water Saturation, Sw", ylabel="Water Fractional Flow, fw")


    if !isnan(s̃w)
        plot!([swi s̃w swj], [fw(swi) fw(s̃w) fw(swj)], labels=["Swi" "S̃w" "Swj"], seriestype=:scatter)
        plot!([swi, s̃w], [fw(swi), fw(s̃w)], ls=:dash, label="Shock")
    else

        plot!([swi swj], [fw(swi) fw(swj)], labels=["Swi" "Swj"], seriestype=:scatter)

    end
end


function plot_fw(μw::Float64, μo::Float64, krs::RelPerms...)
    plot()

    for kr in krs
        fw(sw) = fractional_flow(sw, kr, μw, μo)
        sw = collect(kr.swr:0.01:1 - kr.sor)
        display(plot!(sw, fw(sw), lw=3, label=false))
    end
    plot!(lims=(0, 1))
    plot!(xlabel="Water Saturation, Sw", ylabel="Water Fractional Flow, fw")
end

function plot_sw_profile(wf::WaterFlooding, time::Float64...)
    plot()
    plot!(lims=(0, 1))
    plot!(xlabel="Position, x", ylabel="Water Saturation, Sw")
    for t in time
        x, sw = saturation_profile(wf, t)    
        t = round(t; digits= 2)
        display(plot!(x, sw, lw=3.0, fill=(0, 0.1), label="Time $t"))
    end
end



function animate_sw_profile(wf::WaterFlooding)
    dfw(sw) = fw_derivative(sw, wf.kr, wf.μw, wf.μo)
    tmax = 1 / dfw(wf.S̃)
    times = collect(range(0, 3 * tmax, length=300))
    
    anim = @animate for t in times
        plot_sw_profile(wf::WaterFlooding, t)
    end

    gif(anim, fps=10)
end


# Tracer

function plot_tracer_fw(wf, tracer)
    fw(s) = fractional_flow(s, wf.kr, wf.μw, wf.μo)
    plot_fw(wf)
    plot!([0, tracer.sc], [0, fw(tracer.sc)], marker=:circle, ls=:dot, label="Tracer")
end

function plot_sw_profile(wf::WaterFlooding, c::Tracer, time::Float64)
    plot_sw_profile(wf, time)
    s, x = tracer_profile(wf, c, time)
    xc = c.v * time

    plot!(x, s, lw=0, fill=(0, 0.1, :yellowgreen), label=false)
    plot!([xc, xc], [0, c.sc], lw=2, ls=:dash,  label="Tracer", color=:yellowgreen)

end

function animate_sw_profile(wf::WaterFlooding, c::Tracer)
    dfw(sw) = fw_derivative(sw, wf.kr, wf.μw, wf.μo)
    tmax = 1 / dfw(wf.S̃)
    times = collect(range(0, 3 * tmax, length=300))
    
    anim = @animate for t in times
        plot_sw_profile(wf, c, t)
    end

    gif(anim, fps=10)
end

# * Polymer Flood

function plot_sw_profile(pf::PolymerFlooding, time::Float64...)
    plot()
    for t in time
        x, s = saturation_profile(pf, t)    
        display(plot!(x, s, lw=3.0, fill=(0, 0.1), label="Time $t"))
    end
    plot!(lims=(0, 1))
    plot!(xlabel="Position, x", ylabel="Water Saturation, Sw")
end


function plot_sw_profile(pf::PolymerFlooding)
    plot()
    plot!(lims=(0, 1))
    plot!(xlabel="Position, x", ylabel="Water Saturation, Sw")

    x, s = saturation_profile(pf, t)    
    t = round(t; digits=2)
    plot!(x, s, lw=3.0, fill=(0, 0.1), label="Time $t", title="PVI = $t")
end




function plot_sw_profile(pf::PolymerFlooding, c::Tracer, time::Float64)
    x, s = saturation_profile(pf,time)
    xt, st = tracer_profile(pf, c, time)
    xp, sp = polymer_profile(pf, time)
    xc = c.v * time
    
    plot(x, s, lw=3.0, fill=(0,0.1), label="Water Saturation", lims=(0,1))
    
    plot!(xt, st, lw=0, fill=(0, 0.1, :yellowgreen), label=false)
    plot!([xc, xc], [0, c.sc], lw=3, ls=:dot,  label="Tracer", color=:yellowgreen)
    
    plot!(xp, sp, lw=0, fill=(0, 0.2, :deeppink2), label=false)
    plot!([xp[1], xp[1]], [0, pf.Sb], lw=3, ls=:dot,  label="Polymer", color=:deeppink2)

    plot!(xlabel="Position, x", ylabel="Water Saturation, Sw")
end


function animate_sw_profile(pf::PolymerFlooding, c::Tracer)
    dfw(sw) = fw_derivative(sw, pf.kr, pf.μw, pf.μo)
    tmax = 1 / pf.Vb
    times = collect(range(0, 2 * tmax, length=200))
    
    anim = @animate for t in times
        plot_sw_profile(pf, c, t)
    end

    gif(anim, fps=10)
end


function plot_fw(pf::PolymerFlooding)
    fww(s) = fractional_flow(s, pf.kr2, pf.μw, pf.μo)
    fow(s) = fractional_flow(s, pf.kr1, pf.μw, pf.μo)

    si = pf.Si
    sj = pf.Sj
    s̃ = pf.S̃
    sb = pf.Sb
    swr = pf.kr2.swr
    sor = pf.kr2.sor
    
    sw = collect(swr:0.01:1 - sor)

    plot_fw(pf.μw, pf.μo, pf.kr1, pf.kr2)
    plot!(lims=(0, 1))
    plot!([si sb s̃ sj], [fow(si) fow(sb) fww(s̃)  fww(sj) ], ms=5, labels=["si" "sb" "s*" "sj"], seriestype=:scatter)
    plot!([-pf.D, s̃], [0, fww(s̃)], ls=:dash, color=:limegreen, label="Wettability Shock", lw=2)
    plot!([sb, si], [fow(sb), fow(si)],  lw=2.5, ls=:dash, color=:maroon, label="Oil bank")
    plot!(xlabel="Water Saturation, s", ylabel="Water Fractional Flow, f", legend=:bottomright)
end

# Reactive Transport

function plot_ODEs(it::IonExchangeTransport)
    cₘ₁ = it.cₘ₁
    cⱼ = it.cⱼ
    cₘ₂ = it.cₘ₂

    plot(it.sol2, label=false)
    plot!(it.sol3, label=false)
    plot!([ cₘ₁[3] cₘ₂[3] cⱼ[3]],
    [ cₘ₁[2] cₘ₂[2] cⱼ[2]],
    seriestype=:scatter, labels=["M1" "M2" "J"]
    )
    plot!(xlabel="C3", ylabel="C2")
end


function plot_velocities(it::IonExchangeTransport)
    c = it.c
    λ = it.λ

    plot(λ[2:end-2], it.c[2:end-2,3],
	    label=false, lw=2,
        xlabel="Velocities",
        ylabel="C3, M", marker=:circle)

end