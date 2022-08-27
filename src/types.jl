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

"""
    RelPerms(swr::T, sor::T, krw0::T, kro0, nw::T, no::T) where T <: Float64

Type that defines the water and oil relative permeabilities using Brooks-Corey model.


### Example

```julia
kr = RelPerms(swr = 0.2,
              sor = 0.2,
              krw0 = 0.2,
              kro0 = 0.5,
              nw = 1.5, 
              no = 1.5, 
              )

```
"""
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

"""
    WaterFlooding(;kr::RelPerms, μw::T, μo::T) where T <: Float64

Type that defines the waterflooding problem in terms of the oil and water 
relative permeabilities and viscosities.

### Example
```julia
 kr = RelPerms(swr=0.2,
                    sor=0.2,
                    krw0=0.2,
                    kro0=0.5,
                    nw= nw,
                    no=no)

    μw = 1.0
    μo = 5.0

	wf = WaterFlooding(
					kr=kr,
					μw=μw,
					μo=μo
	)
```

The WaterFlooding struct holds the function to calculate the water
fractional flow:

### Example
```
s = collect(range(kr.swr, 1 - kr.sor); length=100)
f = wf.f(s) 
```

"""
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

"""
ChemicalFlooding(;wf::WaterFlooding, krs::Vector{RelPerms{Float64}}, D::Vector{Float64})

Type that defines the chemical EOR problem where a chemical is injected
to change the fractional flow curve to conditions favorable to oil recovery.

The waterflooding problem is first defined based on the fluids viscosities
and relative permeabilities. Then, a second set of relative permeabilities
is defined in correspondance to the impact of the injected chemical on the 
oil and water mobilities.

More than one chemical wave can be defined. In the following example there
are two chemical waves with their corresponding retardation factors `D`.

### Example
```julia
# Define the waterflooding problem
 kr = RelPerms(swr=0.2,
                    sor=0.2,
                    krw0=0.2,
                    kro0=0.5,
                    nw= nw,
                    no=no
                    )

μw = 1.0
μo = 5.0

wf = WaterFlooding(
                kr=kr,
                μw=μw,
                μo=μo
                )

# Define second relative permeability
kr_2 = RelPerms(swr=0.2,
                sor=0.2,    
                krw0=0.5,
                kro0=1.0,
                nw=3.0,
                no=2.0
                )


kr3 = RelPerms(swr=0.2,
                sor=0.2,    
                krw0=0.2,
                kro0=1.0,
                nw=3.0,
                no=2.0
                )

cf = ChemicalFlooding(
                wf = wf,
                krs = [kr3, kr2],
                D = [0.5, 0.25]
                )

```

"""
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
