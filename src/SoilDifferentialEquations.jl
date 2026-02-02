module SoilDifferentialEquations

# using OrdinaryDiffEqTsit5
import ModelParams: sceua, GOF, of_KGE, of_NSE
using Parameters
using Reexport
# using DiffEqBase
using Printf
using OffsetArrays
using StructArrays
using ProgressMeter

include("ultilize.jl")
include("GlobalOptions.jl")
@reexport using .GlobalOptions

include("soil_texture.jl")
@reexport import .USDA

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

include("SoilParam.jl")
include("Soil.jl")

include("SoilMoisture/SoilMoisture.jl")
include("SPAC.jl")

include("SoilTemperature/soil_properties_thermal.jl")
include("SoilTemperature/soil_temperature.jl")
include("SoilTemperature/soil_temperature_F0.jl")
include("SoilTemperature/EquationTsoil.jl")
include("SoilTemperature/Solve_Tsoil.jl")

include("GroundWater/GroundWater.jl")

include("Vegetation/root_fraction.jl")
include("Vegetation/soil_moisture_constraint.jl")

export solve_Tsoil_Bonan

dir_soil = "$(@__DIR__)/.." |> abspath
export dir_soil

export f_SM_Batesville
dir_root = @__DIR__
f_SM_Batesville = "$(dir_root)/../data/SM_AR_Batesville_8_WNW_2024.csv" |> abspath

export Soil
export VanGenuchten, van_Genuchten, van_Genuchten_K, van_Genuchten_ψ
export soil_moisture!, soil_moisture_Q0!
export soil_temperature!, soil_temperature_F0!
export TsoilEquation, TsoilEquation_partial
export RichardsEquation, RichardsEquation_partial
export soil_properties_thermal

# ModelParams
export sceua, GOF, of_KGE, of_NSE
export soil_depth_init, θ_S, ρ_wat, K0

function _ODEProblem end
function _solve end


end
