module SoilDifferentialEquations
  
# using OrdinaryDiffEq
import HydroTools: θ_S, ρ_wat, ρ_ice, 
  SAND, SILT, CLAY, 
  λ_fus, tfrz, 
  tridiagonal_solver,
  K0, matric_potential
import HydroTools: sceua, GOF, of_KGE, of_NSE

using Parameters
using DiffEqBase
of_MSE(yobs, ysim) = mean((yobs .- ysim) .^ 2)

include("ψ_Cambell.jl")
include("ψ_van_Genuchten.jl")
include("Soil.jl")
include("Soil_depth.jl")

include("soil_moisture.jl")
include("soil_moisture_Q0.jl")

include("soil_thermal_properties.jl")
include("soil_temperature.jl")
include("soil_temperature_F0.jl")

include("Equations/TsoilEquation.jl")
include("Equations/RichardsEquation.jl")

include("Solve_Tsoil.jl")
include("Solve_SM.jl")
include("soil_texture.jl")

export solve_Tsoil_Bonan

dir_soil = "$(@__DIR__)/.." |> abspath
export dir_soil

export Soil
export ParamVanGenuchten, van_Genuchten, van_Genuchten_K, van_Genuchten_ψ
export soil_moisture!, soil_moisture_Q0!
export soil_temperature!, soil_temperature_F0!
export TsoilEquation, TsoilEquation_partial
export RichardsEquation, RichardsEquation_partial
export soil_thermal_properties

# HydroTools
export sceua, GOF, of_KGE, of_NSE
export soil_depth_init, θ_S, ρ_wat, K0

end
