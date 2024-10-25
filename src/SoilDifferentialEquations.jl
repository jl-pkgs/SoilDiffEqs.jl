module SoilDifferentialEquations
  
# using DifferentialEquations
import HydroTools: θ_S, ρ_wat, ρ_ice, 
  SAND, SILT, CLAY, 
  λ_fus, tfrz, 
  tridiagonal_solver,
  van_Genuchten, 
  K0, matric_potential, soil_depth_init
using Parameters

include("Soil.jl")
include("Soil_depth.jl")

include("soil_moisture_Q0.jl")
include("soil_moisture.jl")

include("soil_thermal_properties.jl")
include("soil_temperature.jl")
include("soil_temperature_F0.jl")

include("Equations/TsoilEquation.jl")
include("Equations/RichardsEquation.jl")


dir_soil = "$(@__DIR__)/.." |> abspath
export dir_soil

export Soil, ParamVanGenuchten, van_genuchten_K, van_genuchten_ψ
export soil_moisture!, soil_moisture_Q0!
export TsoilEquation, TsoilEquation_partial, RichardsEquation

# HydroTools
export soil_temperature_F0, soil_temperature
export soil_depth_init, soil_thermal_properties, θ_S, ρ_wat, K0, matric_potential

end
