function plot_sat_profile(sol::Vector, t::Real)
    s, v = get_saturation_speeds(sol)
    plot(v .* t, s, lims=(0,1), lw=3, fill=(0, 0.1), label=false)
    xlabel!("Dimensionless Distance, x")
    ylabel!("Saturation, s")
end


function animate_sat_profile(sol::Vector, t::Real; dt=0.01)
    s, v = get_saturation_speeds(sol)
    
    anim = @animate for i in dt:dt:t
        plot(v .* i, s, lims=(0,1), lw=3, fill=(0, 0.1), label=false)
        xlabel!("Dimensionless Distance, x")
        ylabel!("Saturation, s")
    end
    gif(anim)
end


function plot_fractional_flow(wf::WaterFlooding, sol::Vector)
    f = wf.f

    s = collect(range(wf.kr.swr, 1-wf.kr.sor; length = 101))
    plot(s, f(s), xlims=(0,1), lw=3, alpha=0.4, label=false)

    for (s1, s2, wave, speed, f) in sol
            if wave == :shock
                display(plot!([s1, s2], [f(s1), f(s2)], lw=3, color=:green, alpha=0.5, label=false))
            elseif wave == :spreading
                s = collect(range(s1, s2; length=101))
                display(plot!(s, f(s), alpha=0.5, lw=3, color=:green, label=false))
            end
            display(scatter!([s1, s2], [f(s1), f(s2)], color=:green, alpha=0.5, label=false))
    end

    xlabel!("Dimensionless Distance, x")
    ylabel!("Water Fractional Flow, f")
end


function plot_fractional_flow(cf::ChemicalFlooding, sol::Vector)
    wf = cf.wf

    s = collect(range(wf.kr.swr, 1-wf.kr.sor; length = 101))
    plot(s, wf.f(s), xlims=(0,1), lw=3, alpha=0.4, label=false)


    for (kr, fc) in zip(cf.kr, cf.f)
        s = collect(range(kr.swr, 1-kr.sor; length = 101))
        display(plot!(s, fc(s), label=false, lw=3, alpha=0.4))
    end
    D = copy(cf.D)

    i = 1
    for (s1, s2, wave, speed, f) in sol
            if wave == :shock
                display(plot!([s1, s2], [f(s1), f(s2)], lw=3, color=:green, alpha=0.5, label=false))
                if ~isempty(D)
                    display(plot!([-popfirst!(D), s1], [0, f(s1)], lw=2, color=:black, alpha=0.5, ls=:dot, label=false))
                end
            elseif wave == :spreading
                s = collect(range(s1, s2; length=101))
                display(plot!(s, f(s), alpha=0.5, lw=3, color=:green, label=false))
            end
            display(scatter!([s1, s2], [f(s1), f(s2)], color=:green, alpha=0.5, label=false))
    end
    xlabel!("Dimensionless Distance, x")
    ylabel!("Water Fractional Flow, f")
end



# Reactive Transport Single Phase

function plot_flowing_conc(it, t, labels, colors)
    plot(it.λ * t, it.c, xlim=(0, 1), label=labels, alpha=0.6, lw=3,
            color=colors, markershape=:circle, legend=:outertopright)
    xlabel!("Dimensionless distance, x")
    ylabel!("Flowing concentration eqmol/L")
end


function plot_adsorbed_conc(it, t, labels, colors)
    plot(it.λ * t, it.ĉ[:, begin:end-1], xlim=(0, 1), label=labels, alpha=0.6, lw=3,
            color=colors, markershape=:circle, legend=:outertopright)
    xlabel!("Dimensionless distance, x")
    ylabel!("Adsorbed concentration eqmol/L")
end
