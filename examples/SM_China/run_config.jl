using SoilDifferentialEquations, Ipaper, RTableTools
import ModelParams: sceua
using OrdinaryDiffEqTsit5
include("src/main_plot.jl")

cfg_file = isempty(ARGS) ? joinpath(@__DIR__, "config.yaml") : ARGS[1]
config = load_config(cfg_file)

(; file, scale_factor,
  zs_obs_orgin, zs_obs, z_bound_top,
  zs_center,
  method_retention, method_solve, same_layer, soil_type, dt,
  maxn, objective, of_fun, plot_file) = config

# Config
z₊ₕ = center_to_face(zs_center)
Δz = face_to_thickness(z₊ₕ) ./ 100.0 # [cm] -> [m]

# Load data
d = fread(joinpath(dirname(cfg_file), file))
A_origin = d[:, 2:end] |> Matrix |> drop_missing
A = interp_data_depths(A_origin .* scale_factor, zs_obs_orgin, zs_obs)

# use the global variable: ibeg
function init_soil(; dt=3600.0)
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  N = length(Δz)
  zs_sim = z[1:N] .* 100 # [m] -> [cm]
  ibeg = findfirst(==(z_bound_top), abs.(zs_sim)) + 1 # index of soil model start (simulation layer)
  itop = findfirst(==(z_bound_top), abs.(zs_obs))      # index of obs top boundary layer  

  θ_surf = A[:, itop]        # 边界层（ibeg-1层）的数据作为地表输入
  yobs = A[:, (itop+1):end]  # 从模拟起始层（ibeg层）开始的观测数据 

  θ0 = A[1, :]
  θ = fill(0.2, N)
  n = min(length(θ0), N)
  θ[1:n] .= θ0[1:n]
  θ[n+1:N] .= θ0[end]

  par = get_soilpar(soil_type; method_retention)
  param = SoilParam(N, par; use_m=false, method_retention, same_layer)
  soil = Soil{Float64}(; N, ibeg, dt, z, z₊ₕ, Δz, Δz₊ₕ, θ, method_retention, param)
  cal_ψ!(soil)
  cal_K!(soil)
  soil, θ_surf, yobs
end

soil, θ_surf, yobs = init_soil(; dt)

# Model functions
function model_sim(theta)
  soil, θ_surf, _ = init_soil(; dt)
  SM_UpdateParam!(soil, theta)
  method_solve == "Bonan" ?
  solve_SM_Bonan(soil, θ_surf) :
  solve_SM_ODE(soil, θ_surf; solver=Tsit5()) # [ibeg:N]
end

function goal(theta)
  ysim = model_sim(theta)
  ncol = size(yobs, 2)
  loss = 0.0
  for i in 1:ncol
    obs = view(yobs, :, i)
    sim = view(ysim, :, i)
    loss += -of_fun(obs, sim)
  end
  loss / ncol
end

# Run
if config.optim
  lower, upper = SM_paramBound(soil)
  theta0 = SM_param2theta(soil)

  println("Optimizing (SCE-UA, $objective, maxn=$maxn)...")
  @time theta_opt, feval, _ = sceua(goal, theta0, lower, upper; maxn)
  println("Done. feval=$feval")

  f = joinpath(dirname(cfg_file), "output/theta")
  serialize(f, theta_opt)

  SM_UpdateParam!(soil, theta_opt) # update soil with optimized param
else
  println("Initial loss: $(goal(SM_param2theta(soil)))")
end

# Plot
if !isempty(plot_file)
  theta_plot = SM_param2theta(soil)
  ysim = model_sim(theta_plot)
  dates = d[:, 1]
  depths = round.(Int, -soil.z[soil.inds_obs] .* 100) # [m] -> [cm]
  fout = joinpath(dirname(cfg_file), plot_file)

  plot_result(; ysim, yobs, dates, depths, fout)
end
