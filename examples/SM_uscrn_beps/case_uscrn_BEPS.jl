using SoilDifferentialEquations, Test, Dates, Ipaper
using LazyArtifacts
import RTableTools: fread
# using OrdinaryDiffEqTsit5

function model_sim(theta)
  (; θ_surf, θ0) = options
  soil = init_soil(; θ0, soil_type=8)
  SM_UpdateParam!(soil, theta) # update param
  ysim = soil_moisture_BEPS(soil, θ0, θ_surf)
  return ysim
end

function goal(theta)
  yobs = options.yobs
  ysim = model_sim(theta)
  n = size(ysim, 2)
  ∑ = 0.0
  map(i -> begin
      obs = @view yobs[:, i]
      sim = @view ysim[:, i]
      ∑ += -of_KGE(obs, sim)
    end, 1:n)
  return ∑ / n # mean of NSE
end

begin
  # :SOLARAD, RH_HR_AVG
  vars_SM = [:P_CALC, :SOIL_MOISTURE_5, :SOIL_MOISTURE_10, :SOIL_MOISTURE_20, :SOIL_MOISTURE_50, :SOIL_MOISTURE_100]
  f_uscrn2024 = artifact"USCRN2024" * "/USCRN_hourly_2024_sp54_Apr-Jun.csv"
  df = fread(f_uscrn2024)
  sites = unique_sort(df.site)

  df.time = DateTime.(df.time, "yyyy-mm-ddTHH:MM:SSZ")
end

i = 6
SITE = sites[i]
d = df[df.site.==SITE, [:time; vars_SM]]


z = -[1.25, 5, 10, 20, 50, 100.0] ./ 100# 第一层是虚拟的
function init_soil(; θ0=0.3, dt=3600.0, soil_type=7)
  (; method_retention, same_layer, ibeg) = options
  # dz = [2.5, 5, 5, 15, 45, 55]
  z = -[1.25, 5, 10, 20, 50, 100.0] ./ 100 # 第一层是虚拟的
  Δz = cal_Δz(z)
  N = length(Δz)
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)

  i0 = max(ibeg - 1, 1)
  θ = fill(0.2, N)
  θ[i0:end] .= θ0
  θ[1:i0-1] .= θ0[1]

  par = get_soilpar(soil_type)
  par = VanGenuchten(;
    θ_sat=0.20,
    θ_res=0.03,
    Ksat=1.04,
    α=0.036,
    n=1.56
  )
  method_retention = "van_Genuchten"
  param = SoilParam(N, par;
    use_m=false, method_retention, same_layer)
  soil = Soil{Float64}(; N, ibeg, dt, z, z₊ₕ, Δz, Δz₊ₕ, θ, param, method_retention)
  soil.param.θ_sat .= 0.30
  soil.param.θ_res .= 0.03
  soil
end

set_option!(; same_layer=true, method_retention="van_Genuchten")
options

begin
  d = d[1:24*7, :] # smoke-test
  ibeg = 3
  N_fake = 1
  i0 = max(ibeg - 1, 1) - N_fake
  # [5, 10, 20, 50, 100]
  yobs_full = d[:, 3:end] |> Matrix |> drop_missing # [time, depth], `depth`: 2:N (5cm~100cm)
  yobs = yobs_full[:, i0+1:end]  # 从3th开始模拟(10cm)，2th (5cm)作为输入, 1th虚假层, 
  θ0 = yobs_full[1, i0:end]      # t0初始SM值, 
  θ_surf = yobs_full[:, i0]      # i0所有时刻

  set_option!(; yobs, θ0, θ_surf, ibeg, same_layer=false)
  options

  soil = init_soil(; θ0, soil_type=7)
  soil.param.θ_sat .= 0.30
  soil.param.θ_res .= 0.03

  lower, upper = SM_paramBound(soil)
  theta0 = SM_param2theta(soil)
  @time ysim = model_sim(theta0)
  goal(theta0)
end

# @profview for i = 1:10
#   goal(theta0)
# end

@time theta, feval, exitflag = sceua(goal, theta0, lower, upper; maxn=10_000)
# SM_UpdateParam!(soil, theta)
# include("main_plot.jl")
# plot_result(theta)
# plot_result(theta0)
