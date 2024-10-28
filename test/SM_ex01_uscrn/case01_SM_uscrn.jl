using SoilDifferentialEquations, Test, Dates
using Plots, Printf
import RTableTools: fread
using OrdinaryDiffEq, Ipaper
using Artifacts
# includet("main_Tsoil.jl")
gr(framestyle=:box)


function plot_SM(i; ibeg, ysim=nothing)
  band = vars_SM[i]
  t = d[:, :time]
  x = d[:, band]
  k = i - ibeg + 1

  time_min, time_max = minimum(t), maximum(t)
  ticks = time_min:Dates.Day(7):time_max
  xticks = ticks, format.(ticks, "mm-dd")

  p = plot(title=string(band); xticks,
    xrot=30, tickfonthalign=:center, tickfontvalign=:bottom)
  plot!(p, t, x, label="OBS")
  if k >= 1 && ysim !== nothing
    # i2 = i - 1
    @show i, k
    # title = @sprintf("layer %d: depth = %d cm", i2, -z[i2] * 100)
    plot!(p, t, ysim[:, k], label="SIM")
  end
  return p
end

z = -[1.25, 5, 10, 20, 50, 100.0] ./ 100# 第一层是虚拟的

function init_soil(; θ0, dt=3600.0, ibeg=2, soil_type=7)
  z = -[1.25, 5, 10, 20, 50, 100.0] ./ 100 # 第一层是虚拟的
  Δz = cal_Δz(z)
  z, z₊ₕ, dz₊ₕ = soil_depth_init(Δz)
  # dz = [2.5, 5, 5, 15, 45, 55]
  
  n = length(Δz)
  z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  # m_sat = θ_S[soil_type] * ρ_wat * Δz # kg/m2
  θ = fill(0.2, n)
  θ[ibeg:end] .= θ0
  param_water = get_soilpar(soil_type)
  Soil{Float64}(; n, ibeg, dt, z, z₊ₕ, Δz, Δz₊ₕ, θ, param_water)
end


function model_sim(theta)
  soil = init_soil(; θ0, soil_type=8, ibeg)
  soil.param_water = get_soilpar(theta)
  ysim = solve_SM_Bonan(soil, θ_surf)
  return ysim
end

function goal(theta; kw...)
  ysim = model_sim(theta)
  obs = yobs[:, 2:end]
  sim = ysim[:, 2:end]
  -GOF(obs[:], sim[:]).R2
end


begin
  # :SOLARAD, RH_HR_AVG
  vars_SM = [:P_CALC, :SOIL_MOISTURE_5, :SOIL_MOISTURE_10, :SOIL_MOISTURE_20, :SOIL_MOISTURE_50, :SOIL_MOISTURE_100]

  f_uscrn2024 = artifact"USCRN2024" * "/USCRN_hourly_2024_sp54_Apr-Jun.csv"
  df = fread(f_uscrn2024)
  sites = unique_sort(df.site)

  df.time = DateTime.(df.time, "yyyy-mm-ddTHH:MM:SSZ")
end


begin
  i = 6
  SITE = sites[i]
  d = df[df.site.==SITE, :][1:24*7*12, [:time; vars_SM]]

  ibeg = 2
  # [5, 10, 20, 50, 100]
  yobs_full = d[:, 3:end] |> Matrix |> drop_missing
  yobs = yobs_full[:, max(ibeg-1, 1):end]
  θ0 = yobs_full[1, max(ibeg-1, 1):end]
  θ_surf = yobs_full[:, ibeg-1]

  # θ_sat, θ_res, Ksat, α, n, m
  lower = [0.25, 0.03, 0.002 / 3600, 0.002, 1.05, 0.1]
  upper = [0.50, 0.20, 60.0 / 3600, 0.300, 4.00, 10.0]

  soil = init_soil(; θ0, soil_type=7, ibeg)
  theta0 = soil.param_water |> Vector
  theta = theta0
  goal(theta0)

  f(theta) = goal(theta)
  @time theta, feval, exitflag = sceua(f, theta0, lower, upper; maxn=Int(5e4))

  ysim = model_sim(theta)
  # goal(theta)  
end

get_soilpar(theta)

plot([plot_SM(i; ibeg, ysim) for i in 1:length(vars_SM)]..., size=(1200, 600))
