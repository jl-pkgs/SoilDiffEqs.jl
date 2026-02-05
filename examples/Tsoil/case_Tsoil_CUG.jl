using SoilDifferentialEquations, Ipaper, RTableTools
include("../main_plot.jl")

fileConfig = isempty(ARGS) ? joinpath(@__DIR__, "case_Tsoil_CUG.yaml") : ARGS[1]
config = load_config(fileConfig)

# 加载数据
d = fread(config.file)
dates = d[:, config.col_time]

# CUG 数据深度与 zs_center 完全对应，无需插值
data_obs = Matrix(d[:, config.col_obs_start:end])

# 运行模拟
soil, theta_opt, best_cost = Soil_main(config, data_obs, "CUG", dates; 
  plot_fun=plot_result, maxn=config.maxn)

println("Optimization completed. Best cost: $best_cost")
if config.grid.N <= 5
  println("Optimized parameters:")
  println("  κ = $(soil.param.κ) W/m/K")
  println("  cv = $(soil.param.cv ./ 1e6) MJ/m³/K")
else
  println("Optimized parameters: κ = $(soil.param.κ[1]) W/m/K, cv = $(soil.param.cv[1]/1e6) MJ/m³/K")
end
