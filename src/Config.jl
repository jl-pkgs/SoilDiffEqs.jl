# using SoilDifferentialEquations
using Parameters, YAML

# z: 全部都采用cm
@with_kw struct Config
  ## data
  file::String = ""
  col_time::Int = 1
  col_obs_start::Int = 2

  scale_factor::Float64 = 1.0

  zs_obs_orgin::Vector{Float64} = Float64[]  # 原始观测深度
  zs_obs::Vector{Float64} = Float64[]        # 插值后的观测深度
  z_bound_top::Float64 = 10.0                # top boundary layer depth [cm]

  ## model
  soil_type::Int = 7
  same_layer::Bool = false
  method_retention::String = "van_Genuchten"
  method_solve::String = "Bonan"
  dt::Float64 = 3600.0
  zs_center::Vector{Float64} = Float64[]

  ## optimization
  optim::Bool = false
  maxn::Int = 100
  objective::String = "NSE"
  of_fun::Function = of_NSE

  # output
  plot_file::String = ""
end

function load_config(cfg_file)
  cfg = YAML.load_file(cfg_file)
  data_cfg = get(cfg, "data", Dict())
  model_cfg = get(cfg, "model", Dict())
  opt_cfg = get(cfg, "optimization", Dict())

  ## data
  file = get(data_cfg, "file", "")
  col_time = Int(get(data_cfg, "col_time", 1))
  col_obs_start = Int(get(data_cfg, "col_obs_start", 2))
  scale_factor = Float64(get(data_cfg, "scale_factor", 1.0))
  zs_obs_orgin = Float64.(data_cfg["zs_obs_orgin"])
  zs_obs = Float64.(get(data_cfg, "zs_obs", zs_obs_orgin))
  z_bound_top = Float64(data_cfg["z_bound_top"])

  ## model
  soil_type = Int(get(model_cfg, "soil_type", 7))
  same_layer = Bool(get(model_cfg, "same_layer", false))
  method_retention = get(model_cfg, "method_retention", "van_Genuchten")
  method_solve = get(model_cfg, "method_solve", "Bonan")
  dt = Float64(get(model_cfg, "dt", 3600.0))
  zs_center = Float64.(get(model_cfg, "zs_center", Float64[]))

  ## optimization
  optim = Bool(get(opt_cfg, "enable", false))
  maxn = Int(get(opt_cfg, "maxn", 100))
  objective = cfg["optimization"]["objective"]
  of_fun = objective == "NSE" ? of_NSE : of_KGE

  # output
  plot_file = cfg["output"]["plot_file"]

  Config(;
    file, col_time, col_obs_start, scale_factor, zs_obs_orgin, zs_obs, z_bound_top, # data
    soil_type, same_layer, method_retention, method_solve, dt, zs_center,           # model
    optim, maxn, of_fun, # optimization
    plot_file # output
  )
end


export Config, load_config
