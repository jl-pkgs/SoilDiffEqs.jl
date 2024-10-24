using SoilDifferentialEquations, Plots, Test, RTableTools, Dates
using DifferentialEquations
import Ipaper: drop_missing


include("test-soil_moisture.jl")
include("test-soil_moisture_Q.jl")
include("test-Tsoil.jl")
