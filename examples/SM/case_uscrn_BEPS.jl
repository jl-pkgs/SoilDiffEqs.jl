using SoilDifferentialEquations, Ipaper, RTableTools, Dates, YAML
using LazyArtifacts
import ModelParams: sceua
include("../main_plot.jl")

cfg_file = isempty(ARGS) ? joinpath(@__DIR__, "config_BEPS.yaml") : ARGS[1]
config = load_config(cfg_file)

(; zs_obs_orgin, zs_obs, scale_factor,
  dt, maxn, plot_file, soil_params) = config

begin
  # 加载 USCRN 数据
  site_index = 3  # sites[3] = AR_Batesville_8_WNW
  vars_SM = [:P_CALC, :SOIL_MOISTURE_5, :SOIL_MOISTURE_10, :SOIL_MOISTURE_20, :SOIL_MOISTURE_50, :SOIL_MOISTURE_100]
  f_uscrn2024 = artifact"USCRN2024" * "/USCRN_hourly_2024_sp54_Apr-Jun.csv"
  df = fread(f_uscrn2024)
  sites = unique_sort(df.site)

  SITE = sites[site_index]
  d = df[df.site.==SITE, [:time; vars_SM]]
  d.time = DateTime.(d.time, "yyyy-mm-ddTHH:MM:SSZ")

  d = d[1:24*15, :]
end

# 准备观测数据
data_org = d[:, 2:end] |> Matrix |> drop_missing
data_obs = interp_data_depths(data_org .* scale_factor, zs_obs_orgin, zs_obs)

# 初始化土壤（使用配置中的自定义参数）
soil, θ_surf, yobs = InitSoil(config, data_obs)

begin
  lower, upper = SM_paramBound(soil)
  theta0 = SM_param2theta(soil)

  println("\nOptimizing site '$SITE' (SCE-UA, $(config.objective), maxn=$maxn)...")
  @time theta_opt, feval, _ = sceua(
    theta -> SM_goal(config, data_obs, theta),
    theta0, lower, upper; maxn)

  f = joinpath(@__DIR__, "output", "theta_$(SITE)_BEPS")
  mkpath(dirname(f))
  serialize(f, theta_opt)
  SM_UpdateParam!(soil, theta_opt)

  final_cost = SM_goal(config, data_obs, theta_opt)
  println("\nSite: $SITE, Best $(config.objective): $(-final_cost)")  # 转为正数显示
end

# 绘图
if !isempty(plot_file)
  dates = d[:, 1]
  depths = round.(Int, -soil.z[soil.inds_obs] .* 100)

  ysim, _ = SM_simulate(config, data_obs, theta_opt)
  ysim0, _ = SM_simulate(config, data_obs, theta0)

  img_dir = joinpath(@__DIR__, "images")
  mkpath(img_dir)
  plot_result(; ysim, yobs, dates, depths,
    fout=joinpath(img_dir, "plot_BEPS_$(SITE).png"))

  # 初始参数对比
  plot_result(; ysim=ysim0, yobs, dates, depths,
    fout=joinpath(img_dir, "plot_BEPS_$(SITE)_initial.png"))
end
