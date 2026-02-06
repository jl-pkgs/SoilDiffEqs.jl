using SoilDifferentialEquations, Ipaper, RTableTools, Dates
include("../main_plot.jl")

fileConfig = isempty(ARGS) ? joinpath(@__DIR__, "case_Tsoil_isusm.yaml") : ARGS[1]
config = load_config(fileConfig)

# 加载数据
d = fread(config.file)
dates = DateTime.(d.time, "yyyy-mm-ddTHH:MM:SSZ")

# 深度插值
(; zs_obs_orgin, zs_obs, scale_factor) = config
data_origin = Matrix(d[:, config.col_obs_start:end])
plot_obs(data_origin, dates, zs_obs_orgin; fout="outdir/data_origin.png")

data_obs = interp_data_depths(data_origin .* scale_factor, zs_obs_orgin, zs_obs)

# 运行优化
soil, theta_opt, best_cost = Soil_main(config, data_obs, "isusm", dates;
  plot_fun=plot_result, maxn=config.maxn)

println("Best cost: $best_cost, κ=$(soil.param.κ[1]) W/m/K, cv=$(soil.param.cv[1]/1e6) MJ/m³/K")
