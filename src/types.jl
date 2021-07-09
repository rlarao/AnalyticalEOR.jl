
mutable struct ExchangeConstants
    K₂₁::Float64
    K₃₁::Float64
    K₂₃::Float64
    Z::Float64
    ν::Vector{Int64}
end

mutable struct IonExchangeTransport
    ζᵢ::Vector{Float64}
    ζⱼ::Vector{Float64}
    ν::Vector{Int64}
    K₂₁::Float64
    K₃₁::Float64
    K₂₃::Float64
    Z::Float64
    cᵢ::Vector{Float64}
    cₘ₁::Vector{Float64}
    cₘ₂::Vector{Float64}
    cⱼ::Vector{Float64}
	c::Matrix{Float64} 
	ĉ::Matrix{Float64} 
    λ::Vector{Float64}
    σ::Vector{Float64}
    W2::Symbol
    W3::Symbol
    sol2
    sol3
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
    s::Vector{Float64}
    λ::Vector{Float64}
end




struct Tracer
    phase::Symbol
    sc::Float64
    v::Float64
end

