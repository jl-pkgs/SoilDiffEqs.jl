using YAML, Random, Printf
using SoilDifferentialEquations, Plots, RTableTools
using OrdinaryDiffEqTsit5, OffsetArrays
import ModelParams: sceua, GOF
import NetCDFTools: approx

# =============================================================================
# Helper Functions
# =============================================================================

"""
生成 CLM 指数型分层厚度

# Arguments
- `n::Int`: 层数 (默认 10)
- `fs::Float64`: 比例因子 (默认 0.025)

# Returns
- `Vector{Float64}`: 各层厚度 (m)
"""
function get_clm_layers(n::Int=10; fs::Float64=0.025)::Vector{Float64}
  z_h::Vector{Float64} = [fs * (exp(0.5 * i) - 1) for i in 0:n]
  return diff(z_h)
end

"""
线性插值：在 z_grid 定义的剖面上，求 z_tgt 处的值

# Arguments
- `z_grid::Vector{Float64}`: 网格深度 (m)
- `profile::Vector{Float64}`: 剖面值
- `z_tgt::Float64`: 目标深度 (m)

# Returns
- `Float64`: 插值结果
"""
function interp_profile(z_grid::Vector{Float64}, profile::Vector{Float64}, z_tgt::Float64)::Float64
  z_tgt <= z_grid[1] && return profile[1]
  z_tgt >= z_grid[end] && return profile[end]
  
  for i::Int in 1:length(z_grid)-1
    if z_grid[i] <= z_tgt <= z_grid[i+1]
      r::Float64 = (z_tgt - z_grid[i]) / (z_grid[i+1] - z_grid[i])
      return profile[i] + r * (profile[i+1] - profile[i])
    end
  end
  return NaN
end

"""
求解土壤温度（返回完整剖面）

# Arguments
- `soil::Soil{Float64}`: 土壤结构体
- `Tsurf::Vector{Float64}`: 上边界温度序列
- `solver`: ODE 求解器

# Returns
- `Matrix{Float64}`: 结果矩阵 (时间 × 层数)
"""
function solve_Tsoil(soil::Soil{Float64}, Tsurf::Vector{Float64}; solver)::Matrix{Float64}
  N::Int = soil.N
  ibeg::Int = soil.ibeg
  dt::Float64 = soil.dt
  ntime::Int = length(Tsurf)
  
  # 初始化 ODE 问题
  f(dT, T, p, t) = TsoilEquation_partial(dT, T, p, t; ibeg)
  prob = SoilDifferentialEquations._ODEProblem(f, soil.Tsoil[ibeg:N], (0.0, dt), soil)
  
  # 结果矩阵: 时间 × 层数
  R::Matrix{Float64} = zeros(Float64, ntime, N - ibeg + 1)
  R[1, :] .= soil.Tsoil[ibeg:N]
  
  for i::Int in 2:ntime
    soil.Tsurf = Tsurf[i]
    prob.u0 .= soil.Tsoil[ibeg:N]
    sol = SoilDifferentialEquations._solve(prob, solver; saveat=dt)
    soil.Tsoil[ibeg:N] .= sol.u[end]
    R[i, :] .= soil.Tsoil[ibeg:N]
  end
  return R
end

"""
初始化 Soil 结构体

# Arguments
- `Δz::Vector{Float64}`: 各层厚度 (m)
- `Tsoil0::Vector{Float64}`: 初始温度剖面 (°C)
- `ibeg::Int`: 起始层索引
- `inds_obs::Vector{Int}`: 观测层索引
- `soil_type::Int`: 土壤类型 (默认 1)
- `dt::Float64`: 时间步长 (s, 默认 3600.0)

# Returns
- `Soil{Float64}`: 土壤结构体
"""
function init_soil(Δz::Vector{Float64}, Tsoil0::Vector{Float64}, 
    ibeg::Int, inds_obs::Vector{Int}; 
    soil_type::Int=1, dt::Float64=3600.0)::Soil{Float64}
  
  N::Int = length(Δz)
  _z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  z::OffsetVector{Float64} = _z  # soil_depth_init 返回 OffsetArray
  
  m_sat::Vector{Float64} = θ_S[soil_type] * ρ_wat * Δz
  m_liq::Vector{Float64} = 0.8 * m_sat
  m_ice::Vector{Float64} = 0.0 * m_sat
  
  κ::Vector{Float64}, cv::Vector{Float64} = soil_properties_thermal(Δz, Tsoil0, m_liq, m_ice; soil_type)
  param::SoilParam{Float64} = SoilParam{Float64}(; N, κ, cv)
  
  return Soil{Float64}(; N, ibeg, inds_obs, dt, z, z₊ₕ, Δz, Δz₊ₕ, param, 
    Tsurf=0.0, Tsoil=deepcopy(Tsoil0))
end

# =============================================================================
# Main
# =============================================================================

# 1. 加载配置
fileConfig::String = isempty(ARGS) ? joinpath(@__DIR__, "config.yaml") : ARGS[1]
cfg::Dict{String, Any} = YAML.load_file(fileConfig)
println("Config: $fileConfig")

# 2. 加载数据
data_file::String = cfg["data"]["file"]
isfile(data_file) || (data_file = joinpath(dirname(fileConfig), cfg["data"]["file"]))

d = fread(data_file)
nsteps::Int = cfg["data"]["time_steps"]
t::Vector = d[1:nsteps, cfg["data"]["col_time"]]
yobs_full::Matrix{Float64} = Matrix{Float64}(d[1:nsteps, cfg["data"]["col_obs_start"]:end])
println("Data: $data_file ($nsteps steps)")

# 3. 设置网格
Δz::Vector{Float64} = if haskey(cfg["model"], "layers_thickness")
  Float64.(cfg["model"]["layers_thickness"])
elseif get(cfg["model"], "use_clm_formula", false)
  get_clm_layers(get(cfg["model"], "n_layers", 10))
else
  dz::Float64 = cfg["model"]["dz"]
  fill(dz, ceil(Int, 1.35/dz))
end

nlayer::Int = length(Δz)
z_nodes::Vector{Float64} = cumsum(Δz) .- 0.5Δz  # 各层中心深度

# 观测深度 (m)
z_obs::Vector{Float64} = haskey(cfg["model"], "depths_cm") ? 
  Float64.(cfg["model"]["depths_cm"]) ./ 100 :
  Float64.(cfg["model"]["depths_inch"]) .* 0.0254

# 观测点 -> 网格层索引 (最近邻)
inds_obs::Vector{Int} = [argmin(abs.(z_nodes .- z)) for z in z_obs]
println("Grid: $nlayer layers, Obs depths: $(round.(z_obs.*100, digits=1)) cm")

# 4. 边界条件与初始场
ibeg_obs::Int = cfg["model"]["Tsurf_layer_index"]
ibeg_grid::Int = inds_obs[ibeg_obs]
Tsurf::Vector{Float64} = Float64.(yobs_full[:, ibeg_obs])
Tsoil0::Vector{Float64} = approx(inds_obs, yobs_full[1, :], 1:nlayer)

# 5. 目标函数设置
yobs::Matrix{Float64} = yobs_full[:, ibeg_obs:end]
obs_flat::Vector{Float64} = yobs[:, 2:end][:]
z_tgt::Vector{Float64} = z_obs[ibeg_obs+1:end]
solver = Tsit5()

function goal(theta::Vector{Float64})::Float64
  soil::Soil{Float64} = init_soil(Δz, Tsoil0, ibeg_grid, inds_obs[ibeg_obs:end]; 
    soil_type=cfg["model"]["soil_type"], dt=Float64(cfg["simulation"]["dt"]))
  soil.param.κ .= theta[1:nlayer]
  soil.param.cv .= theta[nlayer+1:end]
  
  R::Matrix{Float64} = solve_Tsoil(soil, Tsurf; solver)
  z_grid::Vector{Float64} = soil.z[ibeg_grid:end]
  
  sim::Matrix{Float64} = [interp_profile(z_grid, R[i, :], z) for i in axes(R, 1), z in z_tgt]
  return -GOF(obs_flat, sim[:]).NSE
end

# 6. 参数优化
soil::Soil{Float64} = init_soil(Δz, Tsoil0, ibeg_grid, inds_obs[ibeg_obs:end]; 
  soil_type=cfg["model"]["soil_type"], dt=Float64(cfg["simulation"]["dt"]))

if cfg["optimization"]["enable"]
  println("Optimizing...")
  Random.seed!(1)
  
  θ0::Vector{Float64} = [soil.param.κ; soil.param.cv]
  lower::Vector{Float64} = [fill(cfg["optimization"]["bounds"]["kappa"][1], nlayer); 
                            fill(cfg["optimization"]["bounds"]["cv"][1], nlayer)]
  upper::Vector{Float64} = [fill(cfg["optimization"]["bounds"]["kappa"][2], nlayer); 
                            fill(cfg["optimization"]["bounds"]["cv"][2], nlayer)]
  
  @time θ_opt::Vector{Float64}, _::Int, best_nse::Float64 = sceua(goal, θ0, lower, upper; 
    maxn=cfg["optimization"]["max_iterations"])
  
  soil.param.κ .= θ_opt[1:nlayer]
  soil.param.cv .= θ_opt[nlayer+1:end]
  @printf("Optimization done. Best cost = %.4f (NSE ≈ %.4f)\n", best_nse, 1 - best_nse)
end

# 7. 最终模拟与绘图
R::Matrix{Float64} = solve_Tsoil(soil, Tsurf; solver)
z_grid::Vector{Float64} = soil.z[ibeg_grid:end]

if cfg["output"]["plot"]
  plts = map(1:size(yobs, 2)) do i::Int
    z::Float64 = z_obs[ibeg_obs + i - 1]
    sim::Vector{Float64} = [interp_profile(z_grid, R[j, :], z) for j in axes(R, 1)]
    plot(t, yobs[:, i], label="Obs", lw=1.5)
    plot!(t, sim, label="Sim", lw=1.5, ls=:dash, 
      title=@sprintf("Depth: %.0f cm", z*100), legend=:topright)
  end
  
  fig = plot(plts..., layout=(ceil(Int, length(plts)/3), 3), 
    size=Tuple(cfg["output"]["plot_size"]))
  savefig(fig, cfg["output"]["save_fig"])
  println("Saved: $(cfg["output"]["save_fig"])")
end
