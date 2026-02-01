using YAML, SoilDifferentialEquations
using Ipaper, RTableTools
import ModelParams: sceua

"""
Config-driven: Soil moisture simulation for SM_China
Usage:
  julia --project=. examples/SM_China/run_config.jl
  julia --project=. examples/SM_China/run_config.jl examples/SM_China/config.yaml
"""

# 1) Load config
cfg_file::String = isempty(ARGS) ? joinpath(@__DIR__, "config.yaml") : ARGS[1]
cfg::Dict{String, Any} = YAML.load_file(cfg_file)
println("Config: $cfg_file")

data_cfg = get(cfg, "data", Dict{String, Any}())
model_cfg = get(cfg, "model", Dict{String, Any}())
opt_cfg = get(cfg, "optimization", Dict{String, Any}())
output_cfg = get(cfg, "output", Dict{String, Any}())

(data_cfg === nothing) && (data_cfg = Dict{String, Any}())
(model_cfg === nothing) && (model_cfg = Dict{String, Any}())
(opt_cfg === nothing) && (opt_cfg = Dict{String, Any}())
(output_cfg === nothing) && (output_cfg = Dict{String, Any}())

# 2) Load observation data
file::String = data_cfg["file"]
if !isfile(file)
  file = joinpath(dirname(cfg_file), file)
end

# 使用已有的 load_data 逻辑
d = fread(file)
time_col::Int = get(data_cfg, "time_col", 1)
obs_start_col::Int = get(data_cfg, "obs_start_col", 2)

# 提取日期和观测数据
dates = d[:, time_col]
_A = d[:, obs_start_col:end] |> Matrix

# 处理缺失值
_A = Ipaper.drop_missing(_A)

# 原始深度和缩放
_depths = Float64.(get(data_cfg, "depths_cm", [10, 20, 30, 40, 50, 60, 80, 100]))
scale_factor::Float64 = get(data_cfg, "scale_factor", 1.0)
_A = _A .* scale_factor

# 插值到目标深度
target_depths = Float64.(get(data_cfg, "target_depths_cm", _depths))

function interp_data(x::AbstractVector, A::AbstractMatrix, xout::AbstractVector)
  ntime = size(A, 1)
  yout = zeros(ntime, length(xout))
  for i = 1:ntime
    y = @view A[i, :]
    yout[i, :] .= approx(x, y, xout)
  end
  yout
end

A = interp_data(_depths, _A, target_depths)
println("Data loaded: $(size(A,1)) time steps, $(size(A,2)) depths")

# 3) Model configuration
ibeg::Int = get(model_cfg, "ibeg", 4)
surface_idx::Int = get(model_cfg, "surface_depth_index", 2)

# 提取输入数据
θ_surf = A[:, surface_idx]           # 表层SM时间序列
θ0 = A[1, :]                         # 初始SM
yobs = A[:, ibeg-1:end]              # 观测值（从ibeg-1开始）

method_retention::String = get(model_cfg, "method_retention", "van_Genuchten")
method_solve::String = get(model_cfg, "method_solve", "Bonan")
same_layer::Bool = get(model_cfg, "same_layer", false)
soil_type::Int = get(model_cfg, "soil_type", 7)
dt::Float64 = get(model_cfg, "dt", 3600.0)

# 设置全局选项
set_option!(; method_retention, yobs, θ_surf, ibeg, same_layer, method_solve)

# 4) Initialize soil
layer_thickness = Float64.(get(model_cfg, "layer_thickness_cm", [5.0, repeat([10], 10)...])) ./ 100.0
Δz = layer_thickness

function init_soil_sm_china(; θ0, dt=3600.0, soil_type=7, method_retention="van_Genuchten", same_layer=false)
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  N = length(Δz)
  
  θ = fill(0.2, N)
  i0 = max(ibeg - 1, 1)
  if length(θ0) >= N
    θ[1:N] .= θ0[1:N]
  else
    θ[1:length(θ0)] .= θ0
    θ[length(θ0)+1:end] .= θ0[end]
  end
  
  par = get_soilpar(soil_type; method_retention)
  param = SoilParam(N, par; use_m=false, method_retention, same_layer)
  soil = Soil{Float64}(; N, ibeg, dt, z, z₊ₕ, Δz, Δz₊ₕ, θ, method_retention, param)
  cal_ψ!(soil)
  cal_K!(soil)
  return soil
end

soil = init_soil_sm_china(; θ0, dt, soil_type, method_retention, same_layer)
println("Soil initialized: N=$(soil.N) layers")

# 5) Define simulation and objective functions
# Load ODE solver if needed
if method_solve != "Bonan"
  using OrdinaryDiffEqTsit5
end

function model_sim(theta)
  soil_sim = init_soil_sm_china(; θ0, dt, soil_type, method_retention, same_layer)
  SM_UpdateParam!(soil_sim, theta)
  
  if method_solve == "Bonan"
    ysim = solve_SM_Bonan(soil_sim, θ_surf)
  else
    solver = Tsit5()
    ysim = solve_SM_ODE(soil_sim, θ_surf; solver)
  end
  return ysim
end

objective_type::String = get(opt_cfg, "objective", "NSE")

function goal(theta)
  ysim = model_sim(theta)
  ncol = size(yobs, 2)
  ∑ = 0.0
  
  # solve_SM_Bonan 返回的 ysim 已经通过 inds_obs 映射到观测层
  # 所以 yobs[:, i] 直接对应 ysim[:, i]
  for i in 1:ncol
    obs = @view yobs[:, i]
    sim = @view ysim[:, i]
    if objective_type == "NSE"
      ∑ += -of_NSE(obs, sim)
    else
      ∑ += -of_KGE(obs, sim)
    end
  end
  return ∑ / ncol
end

# 6) Optimization
if get(opt_cfg, "enable", false)
  maxn::Int = get(opt_cfg, "maxn", 10000)
  
  lower, upper = SM_paramBound(soil)
  theta0 = SM_param2theta(soil)
  
  println("Starting optimization (SCE-UA)...")
  println("Objective: $objective_type, maxn: $maxn")
  @time theta_opt, feval, exitflag = sceua(goal, theta0, lower, upper; maxn)
  println("Optimization done. feval=$feval exit=$exitflag")
  
  # Save results
  theta_file::String = get(output_cfg, "theta_file", "theta")
  if !isabspath(theta_file)
    theta_file = joinpath(dirname(cfg_file), theta_file)
  end
  serialize(theta_file, theta_opt)
  println("Results saved to: $theta_file")
  
  # Update soil with optimized parameters
  SM_UpdateParam!(soil, theta_opt)
else
  # Run with initial parameters
  println("Running with initial parameters...")
  theta0 = SM_param2theta(soil)
  ysim = model_sim(theta0)
  loss = goal(theta0)
  println("Initial objective value: $loss")
end

# 7) Plot results if requested
if get(output_cfg, "save_plot", false)
  include(joinpath(@__DIR__, "main_plot.jl"))
  try
    theta_current = SM_param2theta(soil)
    plot_result(theta_current)
    println("Plot saved")
  catch e
    println("Plotting failed: $e")
  end
end

println("Done!")
