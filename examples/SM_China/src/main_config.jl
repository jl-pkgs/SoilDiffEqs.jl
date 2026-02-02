using SoilDifferentialEquations
using Parameters, YAML

# z: 全部都采用cm
@with_kw struct Config
  file::String = ""
  scale_factor::Float64 = 1.0
  zs_obs_orgin::Vector{Float64} = Float64[]
  zs_obs::Vector{Float64} = Float64[]
  z_bound_top::Float64 = 10.0                # top boundary layer depth [cm]

  soil_type::Int = 7
  same_layer::Bool = false
  method_retention::String = "van_Genuchten"
  method_solve::String = "Bonan"
  dt::Float64 = 3600.0
  zs_center::Vector{Float64} = Float64[]

  enable::Bool = false
  maxn::Int = 100
  objective::String = "NSE"
  of_fun::Function = of_NSE

  plot_file::String = ""
end

function load_config(cfg_file)
  cfg = YAML.load_file(cfg_file)
  data_cfg = get(cfg, "data", Dict())
  model_cfg = get(cfg, "model", Dict())
  opt_cfg = get(cfg, "optimization", Dict())
  output_cfg = get(cfg, "output", Dict())

  zs_obs_orgin = Float64.(get(data_cfg, "zs_obs_orgin", Float64[]))
  zs_obs = Float64.(get(data_cfg, "zs_obs", zs_obs_orgin))
  zs_center = Float64.(get(model_cfg, "zs_center", Float64[]))
  
  Config(
    file=get(data_cfg, "file", ""),
    scale_factor=Float64(get(data_cfg, "scale_factor", 1.0)),
    zs_obs_orgin=zs_obs_orgin,
    zs_obs=zs_obs,
    z_bound_top=Float64(get(model_cfg, "z_bound_top", 10.0)),
    soil_type=Int(get(model_cfg, "soil_type", 7)),
    same_layer=Bool(get(model_cfg, "same_layer", false)),
    method_retention=String(get(model_cfg, "method_retention", "van_Genuchten")),
    method_solve=String(get(model_cfg, "method_solve", "Bonan")),
    dt=Float64(get(model_cfg, "dt", 3600.0)),
    zs_center=zs_center,
    enable=Bool(get(opt_cfg, "enable", false)),
    maxn=Int(get(opt_cfg, "maxn", 100)),
    objective=String(get(opt_cfg, "objective", "NSE")),
    plot_file=String(get(output_cfg, "plot_file", ""))
  )
end
