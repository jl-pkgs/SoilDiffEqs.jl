using SoilDifferentialEquations, Ipaper, RTableTools
include("../main_plot.jl")
plot_fun = plot_result

fileConfig = isempty(ARGS) ? joinpath("examples/SM/case_SM_Bonan_China.yaml") : ARGS[1]
config = load_config(fileConfig)

(; file, zs_obs_orgin, zs_obs, scale_factor) = config

# Load data
SITE = splitext(basename(file))[1]
d = fread(file) # file是相对路径
dates = d[:, 1]

data_origin = d[:, 2:end] |> Matrix |> drop_missing
data_obs = interp_data_depths(data_origin .* scale_factor, zs_obs_orgin, zs_obs)


soil, state, θ_top, yobs = setup(config, data_obs)


maxn = 10000
Soil_main(config, data_obs, SITE, dates; plot_fun, method_retention="van_Genuchten")
Soil_main(config, data_obs, SITE, dates; plot_fun, method_retention="Campbell")
