using SoilDifferentialEquations, Ipaper, RTableTools, Dates, YAML
using LazyArtifacts
include("../main_plot.jl")

cfg_file = isempty(ARGS) ? joinpath(@__DIR__, "case_SM_uscrn.yaml") : ARGS[1]
config = load_config(cfg_file)

(; zs_obs_orgin, zs_obs, scale_factor,
  dt, maxn, objective, of_fun, plot_file) = config

# Load USCRN data from artifact (site_index is USCRN-specific)
begin
  site_index = 3          # 站点索引（sites[3]）
  vars_SM = Symbol.(["SOIL_MOISTURE_5", "SOIL_MOISTURE_10", "SOIL_MOISTURE_20", "SOIL_MOISTURE_50", "SOIL_MOISTURE_100"])

  f_uscrn2024 = artifact"USCRN2024" * "/USCRN_hourly_2024_sp54_Apr-Jun.csv"
  df = fread(f_uscrn2024)
  sites = unique_sort(df.site)

  ## Select site
  SITE = sites[site_index] # 
  d = df[df.site.==SITE, [:time; vars_SM]]
  d.time = DateTime.(d.time, "yyyy-mm-ddTHH:MM:SSZ")
end

# Prepare observation data
data_origin = d[:, 2:end] |> Matrix |> drop_missing
data_obs = interp_data_depths(data_origin .* scale_factor, zs_obs_orgin, zs_obs)


function main(; method_retention=nothing)
  !isnothing(method_retention) && (config.method_retention = method_retention)
  method_retention = config.method_retention

  soil, θ_surf, yobs = InitSoil(config, data_obs)

  ## Run optimization
  # println("Initial loss: $(SM_goal(config, data_obs, SM_param2theta(soil)))")

  lower, upper = SM_paramBound(soil)
  theta0 = SM_param2theta(soil)

  println("Optimizing site '$SITE' $method_retention (SCE-UA, $objective, maxn=$maxn)...")
  @time theta_opt, feval, _ = sceua(theta -> SM_goal(config, data_obs, theta),
    theta0, lower, upper; maxn)

  f = joinpath(dirname(cfg_file), "output/theta_$(SITE)_$(method_retention)")
  mkpath(dirname(f))
  serialize(f, theta_opt)
  SM_UpdateParam!(soil, theta_opt)

  ## Plot
  if !isempty(plot_file)
    dates = d[:, 1]
    # depths from observed layers
    depths = round.(Int, -soil.z[soil.inds_obs] .* 100)

    theta = SM_param2theta(soil)
    ysim, yobs = SM_simulate(config, data_obs, theta)

    fout = joinpath(dirname(cfg_file), "images", "$(method_retention)_$plot_file")
    plot_result(; ysim, yobs, dates, depths, fout)
  end

  println("Site: $SITE, Best Cost: $(SM_goal(config, data_obs, SM_param2theta(soil)))")
end

main(; method_retention="van_Genuchten")
main(; method_retention="Campbell")
