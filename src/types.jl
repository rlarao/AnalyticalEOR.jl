struct RelPerms{T <: Float64}
    swr::T
    sor::T
    krw0::T
    kro0::T
    nw::T
    no::T

    function RelPerms(;swr::T, sor::T, krw0::T, kro0::T, nw::T, no::T) where T <: Float64
        new{Float64}(swr, sor, krw0, kro0, nw, no)
    end
end


struct WaterFlooding{T <: Real}
    S̃::T
    Si::T
    Sj::T
    kr::RelPerms
    μw::T
    μo::T
end 

struct PolymerFlooding{T <: Real}
    S̃::T
    Si::T
    Sj::T
    kr1::RelPerms
    kr2::RelPerms
    μw::T
    μo::T
    Sb::T
    D::T
    Vb::T
    VΔc::T
end




struct Tracer
    phase::Symbol
    sc::Float64
    v::Float64
end