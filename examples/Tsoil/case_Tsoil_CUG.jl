using SoilDifferentialEquations, Ipaper, RTableTools
include("../main_plot.jl")

fileConfig = isempty(ARGS) ? joinpath(@__DIR__, "case_Tsoil_CUG.yaml") : ARGS[1]
config = load_config(fileConfig)

# 加载数据
d = fread(config.file)
dates = d[:, config.col_time]
data_obs = Matrix(d[:, config.col_obs_start:end])

# 运行模拟
model_main(config, data_obs, "CUG", dates; plot_fun=plot_result)
