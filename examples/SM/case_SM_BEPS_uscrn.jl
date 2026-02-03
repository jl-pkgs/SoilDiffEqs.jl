using SoilDifferentialEquations, Ipaper, RTableTools, Dates, LazyArtifacts
include("../main_plot.jl")

fileConfig = isempty(ARGS) ? joinpath(@__DIR__, "case_SM_BEPS_uscrn.yaml") : ARGS[1]
config = load_config(fileConfig)

(; zs_obs_orgin, zs_obs, scale_factor) = config

# Load USCRN data
function load_uscrn_data(site_index::Int)
  vars_SM = [:P_CALC, :SOIL_MOISTURE_5, :SOIL_MOISTURE_10, :SOIL_MOISTURE_20, :SOIL_MOISTURE_50, :SOIL_MOISTURE_100]
  f_uscrn2024 = artifact"USCRN2024" * "/USCRN_hourly_2024_sp54_Apr-Jun.csv"
  df = fread(f_uscrn2024)
  sites = unique_sort(df.site)

  SITE = sites[site_index]
  d = df[df.site.==SITE, [:time; vars_SM]]
  d.time = DateTime.(d.time, "yyyy-mm-ddTHH:MM:SSZ")
  return d, SITE
end

# Main execution
site_index = 3  # sites[3] = AR_Batesville_8_WNW
d, SITE = load_uscrn_data(site_index)
d = d[1:24*15, :]  # Limit to 15 days for faster optimization

# Prepare observation data
data_org = d[:, 2:end] |> Matrix |> drop_missing
data_obs = interp_data_depths(data_org .* scale_factor, zs_obs_orgin, zs_obs)
outdir = joinpath(@__DIR__, "output")

SM_main(config, data_obs, SITE, d.time;
  outdir, plot_fun=plot_result, plot_initial=true)
