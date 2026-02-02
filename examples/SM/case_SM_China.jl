using SoilDifferentialEquations, Ipaper, RTableTools
include("../main_plot.jl")

cfg_file = isempty(ARGS) ? joinpath(@__DIR__, "case_SM_China.yaml") : ARGS[1]
config = load_config(cfg_file)

(; file, scale_factor,
  zs_obs_orgin, zs_obs,
  dt, maxn, objective, of_fun, plot_file) = config

# Load data
d = fread(joinpath(dirname(cfg_file), file))
data_origin = d[:, 2:end] |> Matrix |> drop_missing
data_obs = interp_data_depths(data_origin .* scale_factor, zs_obs_orgin, zs_obs)

soil, θ_surf, yobs = InitSoil(config, data_obs) # θ_surf: boundary layer

# Run
if config.optim
  lower, upper = SM_paramBound(soil)
  theta0 = SM_param2theta(soil)

  println("Optimizing (SCE-UA, $objective, maxn=$maxn)...")
  @time theta_opt, feval, _ = sceua(theta -> SM_goal(config, data_obs, theta),
    theta0, lower, upper; maxn)

  f = joinpath(dirname(cfg_file), "output/theta")
  serialize(f, theta_opt)
  SM_UpdateParam!(soil, theta_opt) # update soil with optimized param
else
  println("Initial loss: $(SM_goal(config, data_obs, SM_param2theta(soil)))")
end

# Plot
if !isempty(plot_file)
  dates = d[:, 1]
  depths = round.(Int, -soil.z[soil.inds_obs] .* 100) # [m] -> [cm]

  theta = SM_param2theta(soil)
  ysim, yobs = SM_simulate(config, data_obs, theta)

  fout = joinpath(dirname(cfg_file), plot_file)
  plot_result(; ysim, yobs, dates, depths, fout)
end

# Iteration =  12, nEvals = 10739, Best Cost = -0.80769
