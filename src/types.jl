struct RelPerms{T <: Real}
    swr::T
    sor::T
    kro0::T
    krw0::T
    nw::T
    no::T

    # function RelPerms(;swr::T, sor::T, krw0::T, kro0::T, nw::T, no::T) where T <: Real
    #     new(swr, sor, krw0, kro0, nw, no)
    # end
end



struct WaterFlooding{T <: Real}
    S̃::T
    Si::T
    Sj::T
    kr::RelPerms
    μw::T
    μo::T

    # function WaterFlooding(Si::T, Sj::T, krs::RelPerms,  μw::T, μo::T) where {T <: Real}
    #     μw, μo = μ
    #     df(s) = fw_derivative(sw, krs, μw, μo)    
    #     Δfw(sw) = (fw(sw) - fw(si)) / (sw - si) 
    #     S̃  = find_zero(sw -> dfw(sw) - Δfw.(sw), (si + 0.0001, sj))
    #     new(S̃, Si, Sj, krs, μw, μo)
    # end
end