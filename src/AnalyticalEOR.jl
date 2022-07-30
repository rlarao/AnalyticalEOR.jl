module AnalyticalEOR

using Base:Real
using ForwardDiff:derivative
using Roots
using Plots
using OrdinaryDiffEq

include("types.jl")
include("flowFunctions.jl")
include("buckleyLeverett.jl")
include("plottingFunctions.jl")
include("reactiveTransport.jl")

export RelPerms, WaterFlooding, ChemicalFlooding, solve_cf
export water_rel_perm, oil_rel_perm, krw_derivative, kro_derivative
export plot_sat_profile, animate_sat_profile, plot_fractional_flow, solve_chemical_flooding
export fractional_flow, fw_derivative, fw_derivative2, fw_derivative3
export solve_waterflooding, get_saturation_speeds, solve_wf, shock_is_admissible, find_sat_shock_spreading
export plot_fw, plot_sw_profile, animate_sw_profile
export plot_tracer_fw
export IonExchangeProblem, isotherm, flowingConcentrations
export solve_IntegralCurve, solve_IntegralCurve2, M2_ODE_solutions, IonExchangedSolution, solve_Ion_Transport

export M2_ODE2, M2_ODE3, M2_ODE_solutions, integralcurves, dc₂dc₃, derivative_functions
export eigenvectors, integral_eigenvalues, RH_eigenvalues
export plot_ODEs, plot_velocities
export isotherm2, flowingConcentrations2
export plot_flowing_conc, plot_adsorbed_conc

end
