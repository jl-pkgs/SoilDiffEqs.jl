using SoilDifferentialEquations, Test, Dates, Ipaper
using LazyArtifacts
import RTableTools: fread
# using OrdinaryDiffEq
includet("main_optim.jl")


begin
  # :SOLARAD, RH_HR_AVG
  vars_SM = [:P_CALC, :SOIL_MOISTURE_5, :SOIL_MOISTURE_10, :SOIL_MOISTURE_20, :SOIL_MOISTURE_50, :SOIL_MOISTURE_100]

  f_uscrn2024 = artifact"USCRN2024" * "/USCRN_hourly_2024_sp54_Apr-Jun.csv"
  df = fread(f_uscrn2024)
  sites = unique_sort(df.site)

  df.time = DateTime.(df.time, "yyyy-mm-ddTHH:MM:SSZ")
end

function plot_result(theta; same_layer=true)
  ysim = model_sim(theta; same_layer)
  plot([plot_sim(i; ibeg, ysim) for i in 1:length(vars_SM)]..., size=(1200, 600))
end


# 每次优化一个参数，可能结果更稳定
begin
  # method = "ODE"
  method = "Bonan"
  same_layer = false

  i = 3
  SITE = sites[i]
  d = df[df.site.==SITE, [:time; vars_SM]]
  # d = d[1:24*7*4, ]

  ibeg = 2
  # [5, 10, 20, 50, 100]
  yobs_full = d[:, 3:end] |> Matrix |> drop_missing
  yobs = yobs_full[:, max(ibeg - 1, 1):end]
  θ0 = yobs_full[1, max(ibeg - 1, 1):end]
  θ_surf = yobs_full[:, ibeg-1]

  soil = init_soil(; θ0, soil_type=7, ibeg, same_layer)
  lower, upper = SM_paramBound(soil)
  # theta0 = soil.param_water |> Vector
  theta0 = SM_param2theta(soil)
  ysim = model_sim(theta0; same_layer)
  goal(theta0; same_layer)

  plot_result(theta0; same_layer)
  f(theta) = goal(theta; method, same_layer, ibeg=1)
  @time theta, feval, exitflag = sceua(f, theta0, lower, upper; maxn=Int(2e4))
end

plot_result(theta; same_layer)
UpdateSoilParam!(soil, theta);
soil
goal(theta; same_layer)

# soil = init_soil(; θ0, soil_type=7, ibeg)
# param = get_soilpar(theta)
# ysim = model_sim(theta; method)
# plot([plot_SM(i; ibeg, ysim) for i in 1:length(vars_SM)]..., size=(1200, 600))
# # goal(theta)

# param = get_soilpar(theta)
# ψ = van_Genuchten_ψ.(θ_surf; param)
# θ2 = [van_Genuchten(x; param)[1] for x in ψ]
# ψ = van_Genuchten_ψ.(θ_surf; param)
# p1 = plot(θ_surf, label="θ_surf")
# plot!(p1, θ2, label="reconstructed θ")
# plot(p1, plot(ψ))
