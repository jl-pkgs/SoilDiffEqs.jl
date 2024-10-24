using SoilDifferentialEquations, Plots, Test, RTableTools, Dates
using DifferentialEquations
import Ipaper: drop_missing


include("test-soil_moisture.jl")
include("test-soil_moisture_Q.jl")
include("test-Tsoil.jl")


dir_root = "$(@__DIR__)/.."
d = fread("$dir_root/data/CUG_TS_202306.csv")

t = d.time
A = Matrix(d[:, 2:end]) |> drop_missing 
TS0 = A[:, 1]
# yobs = A[end, :]
