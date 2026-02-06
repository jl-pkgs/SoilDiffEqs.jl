using Serialization: serialize
using Ipaper: approx

export Soil_predict, Soil_goal, Soil_main
export setup, model_main
export Tsoil_param2theta, Tsoil_UpdateParam!, Tsoil_paramBound


"""
    observe_to_state(state_obs, zs_obs, zs_sim)

将观测数据通过深度插值映射到模拟网格。
使用 `approx` 进行线性插值，将 `zs_obs` 深度的观测值插值到 `zs_sim` 深度。

# Arguments
- `state_obs`: 观测数据向量 (长度 = length(zs_obs))
- `zs_obs`: 观测深度 [cm]
- `zs_sim`: 模拟网格深度 [cm]

# Returns
- `state`: 插值到模拟网格的状态向量 (长度 = length(zs_sim))
"""
function observe_to_state(state_obs::AbstractVector{T}, zs_obs::AbstractVector, zs_sim::AbstractVector) where {T<:Real}
  # 使用 approx 进行深度插值
  state = approx(zs_obs, state_obs, zs_sim)
  return state
end

function set_state!(soil, state, model_type="SM")
  N = soil.N
  if model_type == "SM"
    soil.θ[1:N] .= state
    # 赋予初始值之后，更新ψ/K
    cal_ψ!(soil)
    cal_K!(soil)
  elseif model_type == "Tsoil"
    soil.Tsoil[1:N] .= state
  end
end


function init_SM(config)
  (; dt, soil_type, method_retention, same_layer, soil_params) = config
  (; N, ibeg, z, z₊ₕ, Δz, Δz₊ₕ,) = config.grid

  # 使用自定义参数或标准参数
  if soil_params !== nothing && method_retention == "van_Genuchten"
    @unpack θ_sat, θ_res, Ksat, α, n = soil_params
    par = VanGenuchten{Float64}(; θ_sat, θ_res, Ksat, α, n)
  else
    par = get_soilpar(soil_type; method_retention)
  end

  param = SoilParam(N, par; use_m=false, method_retention, same_layer)
  soil = Soil{Float64}(; dt, N, ibeg, z, z₊ₕ, Δz, Δz₊ₕ,
    method_retention, param)
  return soil
end


function init_Tsoil(config; Tsoil_init=nothing)
  (; dt, soil_type, same_layer) = config
  (; N, ibeg, z, z₊ₕ, Δz, Δz₊ₕ) = config.grid
  inds_obs = config.inds_obs

  # 使用 soil_properties_thermal 物理公式计算 κ 和 cv（与原版一致）
  # 假设 80% 饱和含水量（原版使用的值）
  θ_sat = θ_S[soil_type]
  m_sat = θ_sat * ρ_wat .* Δz  # 每层的饱和含水量 [kg/m²]
  m_ice = zeros(N)
  m_liq = 0.8 .* m_sat  # 80% 饱和（原版使用的值）

  # 如果没有提供初始温度，使用 20°C 作为默认值
  Tsoil_default = isnothing(Tsoil_init) ? fill(20.0, N) : Tsoil_init

  # 使用物理公式计算热力参数
  # κ, cv = soil_properties_thermal(Δz, Tsoil_default, m_liq, m_ice; soil_type)

  κ = fill(1.0, N)  # 热导率 [W/m/K]
  cv = fill(2.0e6, N)  # 热容量 [J/m³/K]

  param = SoilParam{Float64}(; N, κ, cv, same_layer)
  soil = Soil{Float64}(; dt, N, ibeg, z, z₊ₕ, Δz, Δz₊ₕ, param, inds_obs)
  return soil
end


# Model functions
function Soil_predict(config::Config, theta, state, θ_top)
  (; model_type) = config

  if model_type == "SM"
    soil = init_SM(config)
    set_state!(soil, state, model_type)
    SM_UpdateParam!(soil, theta)

    if config.method_solve == "Bonan"
      ysim = solve_SM_Bonan(soil, θ_top)
    elseif config.method_solve == "ODE"
      ysim = solve_SM_ODE(soil, θ_top; solver=Tsit5())
    elseif config.method_solve == "BEPS"
      ysim = soil_moisture_BEPS(soil, θ_top)
    else
      error("Unknown method_solve: $(config.method_solve)")
    end
  elseif model_type == "Tsoil"
    # 使用 state 作为初始温度来计算 κ/cv
    soil = init_Tsoil(config; Tsoil_init=state)
    set_state!(soil, state, model_type)
    Tsoil_UpdateParam!(soil, theta)

    if config.method_solve == "Bonan"
      ysim = solve_Tsoil_Bonan(soil, θ_top)
    elseif config.method_solve == "ODE"
      ysim = solve_Tsoil_ODE(soil, θ_top; solver=Tsit5())
    else
      error("Unknown method_solve: $(config.method_solve)")
    end
  end
  ysim
end


function Soil_goal(config::Config, theta, state, θ_top, yobs)
  of_fun = config.of_fun
  ysim = Soil_predict(config, theta, state, θ_top)
  ncol = size(yobs, 2)
  loss = 0.0
  for i in 1:ncol
    obs = view(yobs, :, i)
    sim = view(ysim, :, i)
    loss += -of_fun(obs, sim)
  end
  loss / ncol
end


function setup(config::Config, data_obs::AbstractMatrix{T}) where {T<:Real}
  (; N, ibeg, itop) = config.grid
  (; zs_obs, zs_center) = config

  # 使用 approx 进行深度插值：将观测深度的数据插值到模拟网格深度
  state = observe_to_state(data_obs[1, :], zs_obs, zs_center)

  θ_top = data_obs[:, itop]         # 边界层（ibeg-1层）的数据作为地表输入
  yobs = data_obs[:, (itop+1):end]  # 从模拟起始层（ibeg层）开始的观测数据 

  if config.model_type == "SM"
    soil = init_SM(config)
  elseif config.model_type == "Tsoil"
    # 使用插值后的初始温度来计算 κ/cv（与原版一致）
    soil = init_Tsoil(config; Tsoil_init=state)
  else
    error("Unknown model_type: $(config.model_type)")
  end
  set_state!(soil, state, config.model_type)

  soil, state, θ_top, yobs
end


function Soil_main(config::Config, data_obs::AbstractMatrix{T}, SITE::AbstractString, dates::AbstractVector;
  maxn::Int=1000,
  method_retention=nothing,
  outdir=nothing, plot_fun=nothing, plot_initial=false, plot_kw...) where {T<:Real}

  ## parameter
  outdir = guess_outdir(config, outdir)

  !isempty(outdir) && mkpath(outdir)
  !isnothing(method_retention) && (config.model_type == "SM" && (config.method_retention = method_retention))
  method_retention = config.model_type == "SM" ? config.method_retention : ""

  ## state and forcing 
  soil, state, θ_top, yobs = setup(config, data_obs)

  ## optimization
  if config.model_type == "SM"
    lower, upper = SM_paramBound(soil)
    theta0 = SM_param2theta(soil)
  elseif config.model_type == "Tsoil"
    lower, upper = Tsoil_paramBound(soil)
    theta0 = Tsoil_param2theta(soil)
  else
    error("Unknown model_type: $(config.model_type)")
  end

  theta_opt, feval, _ = sceua(theta -> Soil_goal(config, theta, state, θ_top, yobs),
    theta0, lower, upper; maxn)

  serialize("$outdir/theta_$(SITE)_$(method_retention)", theta_opt)

  ## visualization
  (; plot_file) = config
  if !isempty(plot_file) && !isnothing(plot_fun)
    depths = round.(Int, -soil.z[soil.inds_obs] .* 100)

    fout = "$outdir/$(method_retention)_$(plot_file)"
    ysim = Soil_predict(config, theta_opt, state, θ_top)
    plot_fun(; ysim, yobs, dates, depths, fout, plot_kw...)

    if plot_initial
      base, ext = splitext(plot_file)
      initial_plot_file = "$(base)_initial$(ext)"

      fout = "$outdir/$(method_retention)_$(initial_plot_file)"
      ysim0 = Soil_predict(config, theta0, state, θ_top)
      plot_fun(; ysim=ysim0, yobs, dates, depths, fout, plot_kw...)
    end
  end

  # update soil params
  if config.model_type == "SM"
    SM_UpdateParam!(soil, theta_opt)
  elseif config.model_type == "Tsoil"
    Tsoil_UpdateParam!(soil, theta_opt)
  end

  best_cost = Soil_goal(config, theta_opt, state, θ_top, yobs)
  return soil, theta_opt, best_cost
end


## Tsoil 参数处理函数
function Tsoil_param2theta(soil::Soil{T}) where {T<:Real}
  N = soil.N
  (; same_layer) = soil.param
  if same_layer
    return [soil.param.κ[1], soil.param.cv[1]]
  else
    return [soil.param.κ; soil.param.cv]
  end
end

function Tsoil_UpdateParam!(soil::Soil{T}, theta::AbstractVector{T}) where {T<:Real}
  N = soil.N
  (; same_layer) = soil.param
  if same_layer
    soil.param.κ .= theta[1]
    soil.param.cv .= theta[2]
  else
    soil.param.κ .= theta[1:N]
    soil.param.cv .= theta[N+1:2N]
  end
  return nothing
end

function Tsoil_paramBound(soil::Soil{T}) where {T<:Real}
  N = soil.N
  (; same_layer) = soil.param
  # κ (热导率): 0.1 - 10 W/m/K
  # cv (热容量): 1e6 - 5e6 J/m³/K
  LOWER = [0.1, 1.0e6]
  UPPER = [10.0, 5.0e6]

  if same_layer
    return LOWER, UPPER
  else
    return repeat(LOWER; inner=N), repeat(UPPER; inner=N)
  end
end
