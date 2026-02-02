using SoilDifferentialEquations, Ipaper, RTableTools
include("../main_plot.jl")

fileConfig = isempty(ARGS) ? joinpath(@__DIR__, "case_SM_Bonan_China.yaml") : ARGS[1]
config = load_config(fileConfig)

# Load data
(; file, zs_obs_orgin, zs_obs, scale_factor) = config
d = fread(joinpath(dirname(fileConfig), file))
data_origin = d[:, 2:end] |> Matrix |> drop_missing
data_obs = interp_data_depths(data_origin .* scale_factor, zs_obs_orgin, zs_obs)

# Site name from file (without extension)
site_name = splitext(basename(file))[1]
outdir = joinpath(dirname(fileConfig), "output")

SM_main(config, data_obs, site_name, d[:, 1];
  outdir, plot_fun=plot_result)
