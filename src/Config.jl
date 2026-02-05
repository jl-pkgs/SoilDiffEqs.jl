# using SoilDifferentialEquations
using Parameters, YAML

# z: 全部都采用cm
@with_kw mutable struct Config
  ## 模型类型: "SM" | "Tsoil"
  model_type::String = "SM"

  ## data
  file::String = ""
  fileConfig::String = ""

  col_time::Int = 1
  col_obs_start::Int = 2

  scale_factor::Float64 = 1.0

  zs_obs_orgin::Vector{Float64} = Float64[]  # 原始观测深度
  zs_obs::Vector{Float64} = Float64[]        # 插值后的观测深度
  z_bound_top::Float64 = 10.0                # top boundary layer depth [cm]
  itop::Int = 1                              # [derived] index of top boundary 

  ## model
  soil_type::Int = 7
  same_layer::Bool = false
  method_retention::String = "van_Genuchten"  # SM 特有
  method_solve::String = "Bonan"              # SM/Tsoil 通用: "Bonan" | "ODE"
  dt::Float64 = 3600.0
  zs_center::Vector{Float64} = Float64[]
  grid::NamedTuple = (;)

  # 自定义土壤参数（可选，用于覆盖标准参数）- SM 特有
  soil_params::Union{Dict{String,Float64},Nothing} = nothing

  ## optimization
  optim::Bool = false
  maxn::Int = 100
  objective::String = "NSE"
  of_fun::Function = of_NSE

  # output
  plot_file::String = ""

  # internal
  config_file::String = ""  # 配置文件路径（用于日志命名）
end


function load_config(fileConfig::String)
  cfg = YAML.load_file(fileConfig)
  data_cfg = get(cfg, "data", Dict())
  model_cfg = get(cfg, "model", Dict())
  opt_cfg = get(cfg, "optimization", Dict())

  ## data
  file = get(data_cfg, "file", "")
  !isabspath(file) && (file = joinpath(dirname(fileConfig), file)) # make sure: abspath

  col_time = Int(get(data_cfg, "col_time", 1))
  col_obs_start = Int(get(data_cfg, "col_obs_start", 2))
  scale_factor = Float64(get(data_cfg, "scale_factor", 1.0))

  # zs_obs_orgin 和 zs_obs 对 SM 是必需的，对 Tsoil 可选
  model_type = get(model_cfg, "model_type", "SM")
  if haskey(data_cfg, "zs_obs_orgin")
    zs_obs_orgin = Float64.(data_cfg["zs_obs_orgin"])
    zs_obs = Float64.(get(data_cfg, "zs_obs", zs_obs_orgin))
  else
    # Tsoil 模式：使用 zs_center 作为默认值
    zs_center_temp = Float64.(get(model_cfg, "zs_center", Float64[]))
    zs_obs_orgin = zs_center_temp
    zs_obs = zs_center_temp
  end

  z_bound_top = Float64(get(data_cfg, "z_bound_top",
    haskey(model_cfg, "z_bound_top") ? model_cfg["z_bound_top"] : 10.0))

  ## model
  # model_type 已在上面解析
  soil_type = Int(get(model_cfg, "soil_type", 7))
  same_layer = Bool(get(model_cfg, "same_layer", false))
  method_retention = get(model_cfg, "method_retention", "van_Genuchten")
  method_solve = get(model_cfg, "method_solve", "Bonan")
  dt = Float64(get(model_cfg, "dt", 3600.0))
  zs_center = Float64.(get(model_cfg, "zs_center", Float64[]))


  # 自定义土壤参数
  soil_params = get(model_cfg, "soil_params", nothing)
  if soil_params !== nothing
    soil_params = Dict{String,Float64}(k => Float64(v) for (k, v) in soil_params)
  end

  ## optimization
  optim = Bool(get(opt_cfg, "enable", false))
  maxn = Int(get(opt_cfg, "maxn", 100))
  objective = cfg["optimization"]["objective"]
  of_fun = objective == "NSE" ? of_NSE : of_KGE

  # output
  plot_file = cfg["output"]["plot_file"]

  ## grid info
  z₊ₕ = center_to_face(zs_center)
  Δz = face_to_thickness(z₊ₕ) ./ 100.0 # [cm] -> [m]
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  N = length(Δz)

  itop = findfirst(==(z_bound_top), abs.(zs_obs))        # index of obs top boundary layer  
  ibeg = findfirst(==(z_bound_top), abs.(zs_center)) + 1 # index of soil modelling start
  grid = (; z, z₊ₕ, Δz, Δz₊ₕ, N, ibeg, itop)

  Config(;
    model_type,  # 模型类型
    file, fileConfig, col_time, col_obs_start, scale_factor, zs_obs_orgin, zs_obs, z_bound_top, itop, # data
    soil_type, same_layer, method_retention, method_solve, dt, zs_center, soil_params, grid, # model
    optim, maxn, objective, of_fun, # optimization
    plot_file, # output
    config_file=fileConfig  # internal
  )
end



function guess_outdir(config::Config, outdir=nothing)
  isnothing(outdir) ? joinpath(dirname(config.file), "outdir") : outdir
end

function open_log(config::Config, log_file=nothing)
  isnothing(log_file) && (log_file = replace(config.fileConfig, r"\.yaml$" => ".log"))
  io = open(log_file, "a")
  return io
end

function log(io, msg)
  println(msg)
  if !isnothing(io)
    println(io, msg)
    flush(io)
  end
end


export Config, load_config
