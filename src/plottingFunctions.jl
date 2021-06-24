
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
    plot!([swi s̃w swj], [fw(swi) fw(s̃w) fw(swj)], labels=["Swi" "S̃w" "Swj"], seriestype=:scatter)
    plot!([swi, s̃w], [fw(swi), fw(s̃w)], ls=:dash, label="Shock")
    plot!(xlabel="Water Saturation, Sw", ylabel="Water Fractional Flow, fw")
end



function plot_sw_profile(wf::WaterFlooding, time::Float64...)
    plot()
    for t in time
        sw, x = saturation_profile(wf, t)    
        display(plot!(x, sw, lw=3.0, fill=(0, 0.1), label="Time $t"))
    end
    plot!(lims=(0, 1))
    plot!(xlabel="Position, x", ylabel="Water Saturation, Sw")
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