module AnalyticalEOR

using Base:Real
using ForwardDiff:derivative
using Roots
using Plots
using DifferentialEquations

include("types.jl")
include("flowFunctions.jl")
include("buckleyLeverett.jl")
include("polymerFlooding.jl")
include("tracers.jl")
include("plottingFunctions.jl")
include("reactiveTransport.jl")

export RelPerms, WaterFlooding
export water_rel_perm, oil_rel_perm, krw_derivative, kro_derivative
export fractional_flow, fw_derivative
export solve_waterflooding, saturation_profile
export solve_tracer, tracer_profile
export plot_fw, plot_sw_profile, animate_sw_profile
export plot_tracer_fw
export solve_polymerflood
export ExchangeConstants, isotherm, flowingConcentrations
export solve_IntegralCurve, solve_IntegralCurve2, M2_ODE_solutions, IonExchangeTransport, solve_Ion_Transport

export M2_ODE2, M2_ODE3, M2_ODE_solutions, integralcurves, dc₂dc₃, derivative_functions
export eigenvectors, integral_eigenvalues, RH_eigenvalues
export plot_ODEs, plot_velocities
export isotherm2, flowingConcentrations2

end
