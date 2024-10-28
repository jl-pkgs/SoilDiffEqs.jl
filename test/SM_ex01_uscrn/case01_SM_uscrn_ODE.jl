using SoilDifferentialEquations, Test, Dates
using Plots, Printf
import RTableTools: fread
using OrdinaryDiffEq, Ipaper
using LazyArtifacts
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

begin
  # :SOLARAD, RH_HR_AVG
  vars_SM = [:P_CALC, :SOIL_MOISTURE_5, :SOIL_MOISTURE_10, :SOIL_MOISTURE_20, :SOIL_MOISTURE_50, :SOIL_MOISTURE_100]

  f_uscrn2024 = artifact"USCRN2024" * "/USCRN_hourly_2024_sp54_Apr-Jun.csv"
  df = fread(f_uscrn2024)
  sites = unique_sort(df.site)

  df.time = DateTime.(df.time, "yyyy-mm-ddTHH:MM:SSZ")
end



function model_sim(theta; method="Bonan")
  soil = init_soil(; θ0, soil_type=8, ibeg)
  soil.param_water = get_soilpar(theta)
  if method == "Bonan"
    ysim = solve_SM_Bonan(soil, θ_surf)
  else
    # solver = Rosenbrock23()
    # solver = Rodas5(autodiff=false)
    solver = Tsit5()
    ysim = solve_SM_ODE(soil, θ_surf; solver)
  end
  return ysim
end

function goal(theta; method="Bonan", kw...)
  ysim = model_sim(theta; method)
  obs = yobs[:, 2:end]
  sim = ysim[:, 2:end]
  -GOF(obs[:], sim[:]).NSE
end


begin
  i = 6
  SITE = sites[i]
  d = df[df.site.==SITE, [:time; vars_SM]]
  # d = d[1:24*7*4, ]

  ibeg = 2
  # [5, 10, 20, 50, 100]
  yobs_full = d[:, 3:end] |> Matrix |> drop_missing
  yobs = yobs_full[:, max(ibeg - 1, 1):end]
  θ0 = yobs_full[1, max(ibeg - 1, 1):end]
  θ_surf = yobs_full[:, ibeg-1]

  # θ_sat, θ_res, Ksat, α, n, m
  lower = [0.25, 0.03, 0.002 / 3600, 0.002, 1.05, 0.1]
  upper = [0.50, 0.20, 60.0 / 3600, 0.300, 4.00, 10.0]

  soil = init_soil(; θ0, soil_type=7, ibeg)
  theta0 = soil.param_water |> Vector
  theta = theta0
  goal(theta0)

  method = "ODE"
  # method = "Bonan"
  f(theta) = goal(theta; method)
  @time theta, feval, exitflag = sceua(f, theta0, lower, upper; maxn=Int(5e4))

  ysim = model_sim(theta; method)
  plot([plot_SM(i; ibeg, ysim) for i in 1:length(vars_SM)]..., size=(1200, 600))
end

soil = init_soil(; θ0, soil_type=7, ibeg)
param = get_soilpar(theta)
ysim = model_sim(theta; method)
plot([plot_SM(i; ibeg, ysim) for i in 1:length(vars_SM)]..., size=(1200, 600))
# goal(theta)

param = get_soilpar(theta)
ψ = van_Genuchten_ψ.(θ_surf; param)
θ2 = [van_Genuchten(x; param)[1] for x in ψ]
ψ = van_Genuchten_ψ.(θ_surf; param)
p1 = plot(θ_surf, label="θ_surf")
plot!(p1, θ2, label="reconstructed θ")
plot(p1, plot(ψ))
