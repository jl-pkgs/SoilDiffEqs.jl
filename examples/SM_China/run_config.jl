using YAML, SoilDifferentialEquations, Ipaper, RTableTools
import ModelParams: sceua

# Load config
cfg_file = isempty(ARGS) ? joinpath(@__DIR__, "config.yaml") : ARGS[1]
cfg = YAML.load_file(cfg_file)

data_cfg = get(cfg, "data", Dict())
model_cfg = get(cfg, "model", Dict())
opt_cfg = get(cfg, "optimization", Dict())
output_cfg = get(cfg, "output", Dict())

# Load data
d = fread(joinpath(dirname(cfg_file), data_cfg["file"]))
_A = d[:, 2:end] |> Matrix |> drop_missing
_A .*= get(data_cfg, "scale_factor", 1.0)

_depths = Float64.(data_cfg["depths_cm"])
target_depths = Float64.(get(data_cfg, "target_depths_cm", _depths))

function interp_data(x, A, xout)
  ntime = size(A, 1)
  yout = zeros(ntime, length(xout))
  for i in 1:ntime
    yout[i, :] .= approx(x, view(A, i, :), xout)
  end
  yout
end

A = interp_data(_depths, _A, target_depths)
println("Data: $(size(A,1))×$(size(A,2))")

# Config
# 自动计算 ibeg：根据 surface_depth_cm 找到对应的 soil layer
surface_depth_cm = model_cfg["surface_depth_cm"]
# 找到 surface_depth_cm 在 target_depths 中的索引
surface_idx = findfirst(==(surface_depth_cm), target_depths)
surface_idx === nothing && error("surface_depth_cm=$surface_depth_cm not found in target_depths")

# 映射关系：target_depths[i] 对应 soil layer (i+1)
# 因为 layer 1 是虚拟层 (-2.5cm)，layer 2 才是第一个真实层 (-10cm)
# 所以 target_depths[1]=10cm 对应 soil layer 2
# surface layer = surface_idx + 1，模拟从下一层开始
ibeg = surface_idx + 2  # +1 因为 target_depths[1] 对应 layer 2，再 +1 从下一层开始

θ_surf = A[:, surface_idx]
θ0 = A[1, :]
# yobs 从 ibeg-1 开始，因为 target_depths[ibeg-1] 对应 soil layer ibeg
yobs = A[:, ibeg-1:end]

method_retention = model_cfg["method_retention"]
method_solve = get(model_cfg, "method_solve", "Bonan")
same_layer = get(model_cfg, "same_layer", false)
soil_type = model_cfg["soil_type"]
dt = model_cfg["dt"]

set_option!(; method_retention, yobs, θ_surf, ibeg, same_layer, method_solve)

# Init soil
Δz = Float64.(model_cfg["layer_thickness_cm"]) ./ 100.0

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
if method_solve != "Bonan"
  using OrdinaryDiffEqTsit5
end

model_sim(theta) = begin
  s = init_soil(; θ0, dt)
  SM_UpdateParam!(s, theta)
  method_solve == "Bonan" ? solve_SM_Bonan(s, θ_surf) : solve_SM_ODE(s, θ_surf; solver=Tsit5())
end

obj_type = get(opt_cfg, "objective", "NSE")

goal(theta) = begin
  ysim = model_sim(theta)
  ncol = size(yobs, 2)
  loss = 0.0
  for i in 1:ncol
    obs = view(yobs, :, i)
    sim = view(ysim, :, i)
    loss += obj_type == "NSE" ? -of_NSE(obs, sim) : -of_KGE(obs, sim)
  end
  loss / ncol
end

# Run
if get(opt_cfg, "enable", false)
  lower, upper = SM_paramBound(soil)
  theta0 = SM_param2theta(soil)
  maxn = opt_cfg["maxn"]

  println("Optimizing (SCE-UA, $obj_type, maxn=$maxn)...")
  @time theta_opt, feval, _ = sceua(goal, theta0, lower, upper; maxn)
  println("Done. feval=$feval")

  f = joinpath(dirname(cfg_file), output_cfg["theta_file"])
  serialize(f, theta_opt)
  println("Saved: $f")

  SM_UpdateParam!(soil, theta_opt)
else
  println("Initial loss: $(goal(SM_param2theta(soil)))")
end

# Plot
if get(output_cfg, "save_plot", false)
  include("main_plot.jl")
  theta_plot = SM_param2theta(soil)
  ysim_plot = model_sim(theta_plot)
  dates_plot = d[:, 1]
  plot_file = get(output_cfg, "plot_file", "plot.png")
  plot_path = joinpath(dirname(cfg_file), plot_file)
  plot_result(; ysim=ysim_plot, yobs, dates=dates_plot, depths=target_depths, ibeg, filename=plot_path)
end
