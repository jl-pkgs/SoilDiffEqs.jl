using SoilDifferentialEquations, Ipaper, RTableTools
import ModelParams: sceua
using OrdinaryDiffEqTsit5
include("src/main_config.jl")
include("src/main_plot.jl")

cfg_file = isempty(ARGS) ? joinpath(@__DIR__, "config.yaml") : ARGS[1]
config = load_config(cfg_file)

(; file, scale_factor,
  zs_obs_orgin, zs_obs,
  zs_center, z_bound_top,
  method_retention, method_solve, same_layer, soil_type, dt,
  enable, maxn, objective, of_fun, plot_file) = config

# Load data
d = fread(joinpath(dirname(cfg_file), file))
A_origin = d[:, 2:end] |> Matrix |> drop_missing
A = interp_data_depths(A_origin .* scale_factor, zs_obs_orgin, zs_obs)

# Config
ibeg = find_ibeg(z_bound_top, zs_obs) # 这是在模型中的位置
itop = ibeg - 1
θ_surf = A[:, itop] # 应该是第一层
yobs = A[:, itop:end]
θ0 = A[1, :]

set_option!(; method_retention, yobs, θ_surf, ibeg, same_layer, method_solve)
z₊ₕ = center_to_face(zs_center)
Δz = face_to_thickness(z₊ₕ) ./ 100.0 # [cm] -> [m]


# use the global variable: ibeg
function init_soil(; θ0, dt=3600.0)
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  N = length(Δz)
  θ = fill(0.2, N)
  n = min(length(θ0), N)
  θ[1:n] .= θ0[1:n]
  θ[n+1:N] .= θ0[end]

  par = get_soilpar(soil_type; method_retention)
  param = SoilParam(N, par; use_m=false, method_retention, same_layer)
  soil = Soil{Float64}(; N, ibeg, dt, z, z₊ₕ, Δz, Δz₊ₕ, θ, method_retention, param)
  cal_ψ!(soil)
  cal_K!(soil)
  soil
end

soil = init_soil(; θ0, dt)

# Model functions
function model_sim(theta)
  soil = init_soil(; θ0, dt)
  SM_UpdateParam!(soil, theta)
  method_solve == "Bonan" ?
  solve_SM_Bonan(soil, θ_surf) :
  solve_SM_ODE(soil, θ_surf; solver=Tsit5())
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
if enable
  lower, upper = SM_paramBound(soil)
  theta0 = SM_param2theta(soil)

  println("Optimizing (SCE-UA, $objective, maxn=$maxn)...")
  @time theta_opt, feval, _ = sceua(goal, theta0, lower, upper; maxn)
  println("Done. feval=$feval")

  f = joinpath(dirname(cfg_file), "output/theta")
  serialize(f, theta_opt)

  SM_UpdateParam!(soil, theta_opt)
else
  println("Initial loss: $(goal(SM_param2theta(soil)))")
end

# Plot
if !isempty(plot_file)
  theta_plot = SM_param2theta(soil)
  ysim_plot = model_sim(theta_plot)
  dates_plot = d[:, 1]
  plot_path = joinpath(dirname(cfg_file), plot_file)
  plot_result(; ysim=ysim_plot, yobs, dates=dates_plot, depths=zs_obs, ibeg, filename=plot_path)
end
