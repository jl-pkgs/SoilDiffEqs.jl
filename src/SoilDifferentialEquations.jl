module SoilDifferentialEquations
  
# using OrdinaryDiffEq
import HydroTools: sceua, GOF, of_KGE, of_NSE
using Parameters
using DiffEqBase

const ρ_wat = 1000.0                       # Density of water, [kg/m3]
const ρ_ice = 917.0                        # Density of ice, [kg/m3]
const λ_fus = 0.3337 * 1e6                 # Heat of fusion for water at 0℃, [J/kg]
const K0 = 273.15                          # zero degree Celsius, [K]
const tfrz = K0                            # freezing Temperature (K)

of_MSE(yobs, ysim) = mean((yobs .- ysim) .^ 2)

# 2.5x faster power method
"Faster method for exponentiation"
pow(x, y) = x^y
# @fastmath pow(x::Real, y::Real) = exp(y * log(x))

include("case_Bonan2019.jl")
include("tridiagonal_solver.jl")

include("Soil.jl")
include("Soil_depth.jl")
include("soil_texture.jl")


include("SoilMoisture/ψ_van_Genuchten.jl")
include("SoilMoisture/ψ_Campbell.jl")
include("SoilMoisture/soil_properties.jl")
include("SoilMoisture/soil_moisture.jl")
include("SoilMoisture/soil_moisture_Q0.jl")
include("SoilMoisture/soil_moisture_BEPS.jl")
include("SoilMoisture/Solve_SM.jl")
include("SoilMoisture/EquationRichards.jl")

include("SoilTemperature/soil_properties_thermal.jl")
include("SoilTemperature/soil_temperature.jl")
include("SoilTemperature/soil_temperature_F0.jl")
include("SoilTemperature/EquationTsoil.jl")
include("SoilTemperature/Solve_Tsoil.jl")

include("GroundWater/GroundWater.jl")

include("GlobalOptions.jl")

export solve_Tsoil_Bonan

dir_soil = "$(@__DIR__)/.." |> abspath
export dir_soil

export f_SM_Batesville
dir_root = @__DIR__
f_SM_Batesville = "$(dir_root)/../data/SM_AR_Batesville_8_WNW_2024.csv" |> abspath

export Soil
export ParamVanGenuchten, van_Genuchten, van_Genuchten_K, van_Genuchten_ψ
export soil_moisture!, soil_moisture_Q0!
export soil_temperature!, soil_temperature_F0!
export TsoilEquation, TsoilEquation_partial
export RichardsEquation, RichardsEquation_partial
export soil_properties_thermal

# HydroTools
export sceua, GOF, of_KGE, of_NSE
export soil_depth_init, θ_S, ρ_wat, K0

end
