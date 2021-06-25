function solve_polymerflood(Si::T, Sj::T, kr1::RelPerms, kr2::RelPerms, D::T, μw::T, μo::T) where T <: Float64
    
    f1(s) = fractional_flow(s, kr1, μw, μo)
    df1(s) = fw_derivative(s, k1, μw, μo)
    f2(s) = fractional_flow(s, kr2, μw, μo)
    df2(s) = fw_derivative(s, kr2, μw, μo)

    # Velocity of polymer shock
    vΔc(s) = f2(s) ./ (s + D)

    # * Step 1 - Calculate Ŝ1 and vΔc
    Ŝ = fzero(s -> vΔc.(s) - df2(s), 0.6)
    VΔc = vΔc(Ŝ)

    # * Step2 - Calculate S1B and v2b
    Sb = fzero(s -> (s + D) * vΔc(Ŝ) - f1(s), 0.4)
    Vb = (f1(Sb) - f1(Si)) / (Sb - Si)

    PolymerFlooding(Ŝ, Si, Sj, kr1, kr2, μw, μo, Sb, D, Vb, VΔc)

end



function saturation_profile(pf::PolymerFlooding, time::Real)
    si = pf.Si
    sj = pf.Sj
    sb = pf.Sb
    ŝ = pf.S̃

    df(s) = fw_derivative(s, pf.kr2, pf.μw, pf.μo)

    s = vcat([si, si, sb ,sb, ŝ], collect(ŝ:0.01:sj))
    x = zeros(length(s))
    x[1] = 1.0
    x[2:3] .= pf.Vb * time
    x[4:5] .= pf.VΔc * time
    x[6:end] = df(s[6:end]) * time

    return x, s
end


function polymer_profile(pf::PolymerFlooding, time::Real)
    sj = pf.Sj
    ŝ = pf.S̃

    df(s) = fw_derivative(s, pf.kr2, pf.μw, pf.μo)

    s = collect(ŝ:0.01:sj)
    x = zeros(length(s))
    x = df(s) * time

    return x, s
end