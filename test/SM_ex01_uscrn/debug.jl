using SoilDifferentialEquations, Plots, Test, RTableTools, Dates
using OrdinaryDiffEq, Ipaper
import HydroTools: sceua, GOF, of_KGE, of_NSE
using Artifacts
includet("main_Tsoil.jl")


begin
  # :SOLARAD, RH_HR_AVG
  vars_SM = [:P_CALC, :SOIL_MOISTURE_5, :SOIL_MOISTURE_10, :SOIL_MOISTURE_20, :SOIL_MOISTURE_50, :SOIL_MOISTURE_100]

  f_uscrn2024 = artifact"USCRN2024" * "/USCRN_hourly_2024_sp54_Apr-Jun.csv"
  df = fread(f_uscrn2024)
  sites = unique_sort(df.site)

  df.time = DateTime.(df.time, "yyyy-mm-ddTHH:MM:SSZ")

  SITE = sites[i]
  d = df[df.site.==SITE, :][1:24*7*8, [:time; vars_SM]]
end
 
yobs_full = d[:, 3:end] |> Matrix |> drop_missing
