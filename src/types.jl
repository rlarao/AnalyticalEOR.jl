using Roots, ForwardDiff

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
    si::T
    sj::T
    kr::RelPerms
    f::Function
    df::Function
    d²f::Function

    function WaterFlooding(;si::T, sj::T, kr::RelPerms, μw::T, μo::T) where T <: Float64

        f(s) = fractional_flow(s, kr, μw, μo)
        df(s) = ForwardDiff.derivative.(s -> f(s), s)
        d²f(s) = ForwardDiff.derivative.(s -> df(s), s)
        new{Float64}(si, sj, kr, f, df, d²f)
    end
end 


struct ChemicalFlooding{T <: Real}
    sj::T
    kr::Vector{RelPerms}
    f::Vector{Function}
    df::Vector{Function}
    d²f::Vector{Function}
    D::Vector{T}
    wf::WaterFlooding

    function ChemicalFlooding(;si::T, sj::T, kr::RelPerms, μw::T, μo::T) where T <: Float64
        f(s) = fractional_flow(s, kr, μw, μo)
        df(s) = ForwardDiff.derivative.(s -> f(s), s)
        d²f(s) = ForwardDiff.derivative.(s -> df(s), s)
        new{Float64}(si, sj, kr, μw, μo, f, df, d²f)
    end
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

