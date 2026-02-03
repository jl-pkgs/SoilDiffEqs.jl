using SoilDifferentialEquations, Ipaper, RTableTools, Dates
using LazyArtifacts
include("../main_plot.jl")

fileConfig = isempty(ARGS) ? joinpath(@__DIR__, "case_SM_Bonan_uscrn.yaml") : ARGS[1]
config = load_config(fileConfig)

# Load USCRN data from artifact
function load_uscrn_data(site_index::Int)
  vars_SM = Symbol.(["SOIL_MOISTURE_5", "SOIL_MOISTURE_10", "SOIL_MOISTURE_20", "SOIL_MOISTURE_50", "SOIL_MOISTURE_100"])

  f_uscrn2024 = artifact"USCRN2024" * "/USCRN_hourly_2024_sp54_Apr-Jun.csv"
  df = fread(f_uscrn2024)
  sites = unique_sort(df.site)

  SITE = sites[site_index]
  d = df[df.site.==SITE, [:time; vars_SM]]
  d.time = DateTime.(d.time, "yyyy-mm-ddTHH:MM:SSZ")
  return d, SITE
end

# Main execution
(; zs_obs_orgin, zs_obs, scale_factor) = config

site_index = 3  # 站点索引
d, SITE = load_uscrn_data(site_index)
data_origin = d[:, 2:end] |> Matrix |> drop_missing
data_obs = interp_data_depths(data_origin .* scale_factor, zs_obs_orgin, zs_obs)


plot_fun = plot_result
SM_main(config, data_obs, SITE, d.time; plot_fun, method_retention="van_Genuchten")
SM_main(config, data_obs, SITE, d.time; plot_fun, method_retention="Campbell")
