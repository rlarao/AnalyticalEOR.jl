module AnalyticalEOR

using Base:Real
using ForwardDiff:derivative
using Roots
using Plots
using OrdinaryDiffEq

include("types.jl")
include("flowFunctions.jl")
include("solvers.jl")
include("plottingFunctions.jl")
include("reactiveTransport.jl")

export RelPerms, WaterFlooding, ChemicalFlooding, IonExchangeProblem
export get_saturation_speeds, solve_wf, solve_cf
export plot_sat_profile, animate_sat_profile, plot_fractional_flow
export plot_flowing_conc, plot_adsorbed_conc


end
