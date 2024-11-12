using SoilDifferentialEquations, Test, Dates, Ipaper
using LazyArtifacts
import RTableTools: fread
# using OrdinaryDiffEq
includet("main_optim.jl")
include("main_plot.jl")

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
options.method_retention = "van_Genuchten"


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
  set_option!(; yobs, θ0, θ_surf, ibeg, same_layer=false)
  options

  soil = init_soil(; θ0, soil_type=7)
  soil.param.θ_sat .= 0.30
  soil.param.θ_res .= 0.03

  lower, upper = SM_paramBound(soil)
  theta0 = SM_param2theta(soil)
  @time ysim = model_sim(theta0;)
end



begin
  lower, upper = SM_paramBound(soil)
  theta0 = SM_param2theta(soil)
  ysim = model_sim(theta0;)
  goal(theta0;)
  # plot_result(theta0)
  f(theta) = goal(theta; ibeg=1)
  @time theta, feval, exitflag = sceua(f, theta0, lower, upper; maxn=Int(2e4))
end

ysim = model_sim(theta0;)
plot_result(ysim)
# UpdateSoilParam!(soil, theta);
# soil
# goal(theta; same_layer)
