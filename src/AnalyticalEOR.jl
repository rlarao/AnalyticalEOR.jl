module AnalyticalEOR

using Base:Real
using ForwardDiff
using Roots
using Plots

include("types.jl")
include("flowFunctions.jl")
include("buckleyLeverett.jl")
include("polymerFlooding.jl")
include("tracers.jl")
include("plottingFunctions.jl")

export RelPerms, WaterFlooding
export water_rel_perm, oil_rel_perm, krw_derivative, kro_derivative
export fractional_flow, fw_derivative
export solve_waterflooding, saturation_profile
export solve_tracer, tracer_profile
export plot_fw, plot_sw_profile, animate_sw_profile
export plot_tracer_fw
export solve_polymerflood

end
