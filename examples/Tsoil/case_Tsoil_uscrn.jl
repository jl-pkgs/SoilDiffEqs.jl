using SoilDifferentialEquations, Ipaper, RTableTools, Dates
using LazyArtifacts
include("../main_plot.jl")

fileConfig = isempty(ARGS) ? joinpath(@__DIR__, "case_Tsoil_uscrn_03.yaml") : ARGS[1]
config = load_config(fileConfig)

# 从 artifact 加载 USCRN 数据
f_uscrn = artifact"USCRN2024" * "/USCRN_hourly_2024_sp54_Apr-Jun.csv"
df = fread(f_uscrn)
sites = unique_sort(df.site)

# 站点选择 - 简化处理，只测试第一个站点
SITE = sites[1]
println("="^60)
println("Processing site: $SITE")
println("="^60)

# 提取站点数据 (默认14天)
d_site = df[df.site.==SITE, :][1:min(24*14, nrow(df[df.site.==SITE, :])), :]

# USCRN 变量列
vars_TS = [:SUR_TEMP, :SOIL_TEMP_5, :SOIL_TEMP_10, :SOIL_TEMP_20, :SOIL_TEMP_50, :SOIL_TEMP_100]
d = d_site[:, [:time; vars_TS]]
dates = d.time

# 转换为矩阵 (跳过 time 列)，并处理缺失值
data_raw = Matrix(d[:, 2:end])
# 将 Missing 替换为 NaN 并转换为 Float64
data_obs = Matrix{Float64}(replace(data_raw, missing => NaN))
# 移除包含 NaN 的行
valid_rows = .!any(isnan.(data_obs), dims=2)[:]
data_obs = data_obs[valid_rows, :]
dates = dates[valid_rows]
# 转换日期 (处理 ISO8601 格式: 2024-04-01T00:00:00Z)
dates = [DateTime(String(ds)[1:end-1]) for ds in dates]

# 运行模拟
soil, theta_opt, best_cost = Soil_main(config, data_obs, String(SITE), dates; 
  plot_fun=plot_result, maxn=config.maxn)

println("Site $SITE processed!")
println("Optimization completed. Best cost: $best_cost")
if config.grid.N <= 5
  println("Optimized parameters:")
  println("  κ = $(soil.param.κ) W/m/K")
  println("  cv = $(soil.param.cv ./ 1e6) MJ/m³/K")
else
  println("Optimized parameters: κ = $(soil.param.κ[1]) W/m/K, cv = $(soil.param.cv[1]/1e6) MJ/m³/K")
end
