module SoilDifferentialEquations
  
# using DifferentialEquations
using HydroTools
using Parameters

include("Soil.jl")
include("Soil_depth.jl")
include("soil_moisture_Q0.jl")
include("soil_moisture.jl")
include("Equations.jl")

export Soil, ParamVanGenuchten, van_genuchten_K, van_genuchten_ψ
export soil_moisture!, soil_moisture_Q0!
export TsoilEquation, RichardsEquation

# HydroTools
export soil_temperature
export soil_depth_init, soil_thermal_properties, θ_S, ρ_wat, K0, matric_potential

end
