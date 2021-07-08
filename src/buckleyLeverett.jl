function solve_waterflooding(si, sj, kr, μw, μo)
    fw(sw) = fractional_flow(sw, kr, μw, μo)
    dfw(sw) = fw_derivative(sw, kr, μw, μo)
    Δfw(sw) = (fw(sw) - fw(si)) / (sw - si) 

    S̃  = fzero(sw -> dfw(sw) - Δfw.(sw), sj*0.9)

    WaterFlooding(S̃, si, sj, kr, μw, μo)
end




function saturation_profile(wf::WaterFlooding, time::Real)
    si = wf.Si
    sj = wf.Sj
    s̃ = wf.S̃

    fw(sw) = fractional_flow(sw, wf.kr, wf.μw, wf.μo)
    dfw(sw) = fw_derivative(sw, wf.kr, wf.μw, wf.μo)
    Δfw(sw) = (fw(sw) - fw(si)) / (sw - si) 

    s = vcat([si,si], collect(s̃:0.01:sj), [sj])
    x = Vector{Float64}(undef, length(s))
    x[1] = 1.0
    x[2:3] .= Δfw(s̃) * time
    x[findall(s .> s̃)] = dfw(s[s .> s̃]) * time
    x[end] = 0.0

    return x, s
end




# struct WaterFlooding{T <: Real}
#     S̃::T
#     Si::T
#     Sj::T
#     RelPerms::RelPerms
#     μw::T
#     μo::T

#     function WaterFlooding(Si::T, Sj::T, krs::RelPerms,  μw::T, μo::T) where {T <: Real}
#         μw, μo = μ
#         df(s) = fw_derivative(sw, krs, μw, μo)    
#         Δfw(sw) = (fw(sw) - fw(si)) / (sw - si) 
#         S̃  = find_zero(sw -> dfw(sw) - Δfw.(sw), (si + 0.0001, sj))
#         new(S̃, Si, Sj, krs, μw, μo)
#     end
# end


