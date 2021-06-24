function solve_tracer(wf::WaterFlooding; phase::Symbol=:water)
    si = wf.Si
    sj = wf.Sj
    s̃ = wf.S̃
    
    fw(s) = fractional_flow(s, wf.kr, wf.μw, wf.μo)
    dfw(s) = fw_derivative(s, wf.kr, wf.μw, wf.μo)
    Δfw(s) = (fw(s) - fw(si)) / (s - si) 

    # vw(s) = s < s̃ ? Δfw(s̃) : dfw(s)
    vc(s) = fw(s) ./ s

    sc = fzero(s -> vc.(s) - dfw(s), (si + sj) / 2)

    return Tracer(:water, sc, vc(sc)) 
end



function tracer_profile(wf::WaterFlooding, c::Tracer, time::Real)
    si = wf.Si
    sj = wf.Sj
    s̃ = wf.S̃
    sc = c.sc

    fw(sw) = fractional_flow(sw, wf.kr, wf.μw, wf.μo)
    dfw(sw) = fw_derivative(sw, wf.kr, wf.μw, wf.μo)
    Δfw(sw) = (fw(sw) - fw(si)) / (sw - si) 

    s = vcat(collect(sc:0.01:sj), [sj])
    x = Vector{Float64}(undef, length(s))
    x[1] = c.v * time
    x[2:end - 1] = dfw(s[2:end - 1]) * time
    x[end] = 0.0

    return s, x
end