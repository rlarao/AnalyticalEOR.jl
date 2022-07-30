using Roots, ForwardDiff

mutable struct IonExchangeProblem
    K₂₁::Float64
    K₃₁::Float64
    K₂₃::Float64
    Z::Float64
    ν::Vector{Int64}
end

mutable struct IonExchangeSolution
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
    kr::RelPerms
    μw::Float64
    μo::Float64
    f::Function
    df::Function
    d²f::Function

    function WaterFlooding(;kr::RelPerms, μw::T, μo::T) where T <: Float64
        f(s) = fractional_flow(s, kr, μw, μo)
        df(s) = ForwardDiff.derivative.(s -> f(s), s)
        d²f(s) = ForwardDiff.derivative.(s -> df(s), s)
        new{Float64}(kr, μw, μo, f, df, d²f)
    end
end 


struct ChemicalFlooding{T <: Real}
    kr::Vector{RelPerms{Float64}}
    f::Vector{Function}
    df::Vector{Function}
    d²f::Vector{Function}
    D::Vector{T}
    wf::WaterFlooding

    function ChemicalFlooding(;wf::WaterFlooding, krs::Vector{RelPerms{Float64}}, D::Vector{T}) where T <: Float64
        fs = [s -> fractional_flow(s, kr, wf.μw, wf.μo) for kr in krs]
        dfs = [s -> ForwardDiff.derivative.(s -> f(s), s) for f in fs]
        d²fs = [s -> ForwardDiff.derivative.(s -> df(s), s) for df in dfs]
        new{Float64}(krs, fs, dfs, d²fs, D, wf)
    end
end 
