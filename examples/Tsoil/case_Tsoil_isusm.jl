using SoilDifferentialEquations, Ipaper, RTableTools, Dates
include("../main_plot.jl")

fileConfig = isempty(ARGS) ? joinpath(@__DIR__, "case_Tsoil_isusm.yaml") : ARGS[1]
config = load_config(fileConfig)

# 加载数据
d = fread(config.file)
dates_raw = d[:, config.col_time]

# 转换日期类型 (处理 ISO8601 格式: 2022-07-01T00:00:00Z)
if dates_raw isa Vector{<:Dates.TimeType}
  dates = dates_raw
else
  date_strs = String.(dates_raw)
  # 移除 Z 后缀并解析
  dates = [DateTime(ds[1:end-1]) for ds in date_strs]
end

# 提取原始观测数据并进行深度插值
# ISUSM 数据: 从第4列开始是土壤温度数据 (sv_t04, sv_t12, ...)
data_origin = Matrix(d[:, config.col_obs_start:end])

# 深度插值：将观测深度映射到模拟深度
(; zs_obs_orgin, zs_obs, scale_factor) = config
data_obs = interp_data_depths(data_origin .* scale_factor, zs_obs_orgin, zs_obs)

# 运行模拟
soil, theta_opt, best_cost = Soil_main(config, data_obs, "isusm", dates; 
  plot_fun=plot_result, maxn=config.maxn)

println("Optimization completed. Best cost: $best_cost")
if config.grid.N <= 5
  println("Optimized parameters:")
  println("  κ = $(soil.param.κ) W/m/K")
  println("  cv = $(soil.param.cv ./ 1e6) MJ/m³/K")
else
  println("Optimized parameters: κ = $(soil.param.κ[1]) W/m/K, cv = $(soil.param.cv[1]/1e6) MJ/m³/K")
end
