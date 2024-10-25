module SoilDifferentialEquations
  
# using DifferentialEquations
using HydroTools
using Parameters

import HydroTools: soil_temperature_delta, soil_temperature, soil_thermal_properties

include("Soil.jl")
include("Soil_depth.jl")
include("soil_moisture_Q0.jl")
include("soil_moisture.jl")
include("Equations/TsoilEquation.jl")
include("Equations/RichardsEquation.jl")


export dir_soil
dir_soil = "$(@__DIR__)/.." |> abspath

export Soil, ParamVanGenuchten, van_genuchten_K, van_genuchten_ψ
export soil_moisture!, soil_moisture_Q0!

export TsoilEquation, TsoilEquation_partial, RichardsEquation

# HydroTools
export soil_temperature
export soil_depth_init, soil_thermal_properties, θ_S, ρ_wat, K0, matric_potential

end
