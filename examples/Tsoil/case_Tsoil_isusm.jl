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
data_obs = Matrix(d[:, config.col_obs_start:end])

# 运行模拟
model_main(config, data_obs, "isusm", dates; plot_fun=plot_result)
