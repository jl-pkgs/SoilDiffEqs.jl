using Plots, Printf
gr(framestyle=:box)

function plot_sim(i; ibeg, ysim=nothing)
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

function plot_result(theta; same_layer=true)
  ysim = model_sim(theta; same_layer)
  plot([plot_sim(i; ibeg, ysim) for i in 1:length(vars_SM)]..., size=(1200, 600))
end


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


GlobalOptions.options = Options()
options = GlobalOptions.options


begin
  d = fread(f_SM_Batesville)
  ibeg = 2
  yobs_full = d[:, 3:end] |> Matrix #|> drop_missing
  yobs = yobs_full[:, max(ibeg - 1, 1):end]
  θ0 = yobs_full[1, max(ibeg - 1, 1):end]
  θ_surf = yobs_full[:, ibeg-1]

  set_option!(; yobs, θ_surf, ibeg, same_layer=true)
  options
end

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
  goal(theta0; same_layer, yobs)

  plot_result(theta0; same_layer)
  f(theta) = goal(theta; method, same_layer, ibeg=1)
  @time theta, feval, exitflag = sceua(f, theta0, lower, upper; maxn=Int(2e4))
end

# plot_result(theta; same_layer)
# UpdateSoilParam!(soil, theta);
# soil
# goal(theta; same_layer)
