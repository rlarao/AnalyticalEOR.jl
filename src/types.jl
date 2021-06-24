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
end 

struct Tracer
    phase::Symbol
    sc::Float64
    v::Float64
end