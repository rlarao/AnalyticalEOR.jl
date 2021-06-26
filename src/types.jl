
struct ExchangeConstants
    K₂₁::Float64
    K₃₁::Float64
    K₂₃::Float64
    Z::Float64
end

struct IonExchangeTransport{T <: Float64}
    ζᵢ::Vector{T}
    ζⱼ::Vector{T}
    ν::Vector{Int64}
    K₂₁::T
    K₃₁::T
    K₂₃::T
    Z::T
    cᵢ::Vector{T}
    cₘ₁::Vector{T}
    cₘ₂::Vector{T}
    cⱼ::Vector{T}
    λ₁::T
	λ₂::Vector{T} 
	λ₃::Vector{T}
    c₁::Vector{T}
    c₂::Vector{T}
    c₃::Vector{T}
    c₄::Vector{T}
    λ::Vector{T}
    
end

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

