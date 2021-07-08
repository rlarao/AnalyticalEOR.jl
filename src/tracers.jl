function solve_tracer(wf::WaterFlooding; phase::Symbol=:water)
    si = wf.Si
    sj = wf.Sj
    s̃ = wf.S̃
    
    fw(s) = fractional_flow(s, wf.kr, wf.μw, wf.μo)
    dfw(s) = fw_derivative(s, wf.kr, wf.μw, wf.μo)
    Δfw(s) = (fw(s) - fw(si)) / (s - si) 

    # vw(s) = s < s̃ ? Δfw(s̃) : dfw(s)
    vc(s) = fw(s) ./ s

    sc = fzero(s -> vc.(s) - dfw(s), 0.9*sj)

    return Tracer(:water, sc, vc(sc)) 
end

function solve_tracer(pf::PolymerFlooding; phase::Symbol=:water)
    sb = pf.Sb
    f(s) = fractional_flow(s, pf.kr1, pf.μw, pf.μo)
    vc = f(sb) / sb

    return Tracer(:water, sb, vc) 
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

    return x, s
end


function tracer_profile(pf::PolymerFlooding, c::Tracer, time::Real)
    si = pf.Si
    sj = pf.Sj
    sb = pf.Sb
    ŝ = pf.S̃
    sc = c.sc
    vc = c.v

    df(s) = fw_derivative(s, pf.kr2, pf.μw, pf.μo)

    s = vcat([sc, sb ,sb, ŝ], collect(ŝ:0.01:sj))
    x = zeros(length(s))
    x[1] = vc * time
    x[2:3] .= pf.VΔc * time
    x[4:end] = df(s[4:end]) * time

    return x, s
end