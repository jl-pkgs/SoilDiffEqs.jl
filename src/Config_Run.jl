using Serialization: serialize
export Soil_predict, Soil_goal, Soil_main
export setup


function observe_to_state(state_obs, N::Int, ibeg::Int, itop::Int)
  N_obs = length(state_obs)
  N_input = min(N_obs - itop + 1, N) # itop:N_obs
  N_need = N - ibeg + 2 # 从ibeg-1到N

  state = zeros(eltype(state_obs), N)
  state[1:ibeg-2] .= state_obs[1] # 上边界层全部赋值为第一个观测值

  if N_input == N_need
    state[ibeg-1:N] = state_obs[itop:end]
  else
    # N_input = x - (ibeg - 1) + 1 ==> x = N_input + ibeg - 2
    state[ibeg-1:N_input+ibeg-2] .= state_obs[itop:end] # len = N_input
    state[N_input+ibeg-1:N] .= state_obs[end]
  end
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
    # soil = init_Tsoil(config)
    # set_state!(soil, state, model_type)
    # Tsoil_UpdateParam!(soil, theta)
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
  state = observe_to_state(data_obs[1, :], N, ibeg, itop)
  θ_top = data_obs[:, itop]         # 边界层（ibeg-1层）的数据作为地表输入
  yobs = data_obs[:, (itop+1):end]  # 从模拟起始层（ibeg层）开始的观测数据 

  soil = init_SM(config)
  set_state!(soil, state, config.model_type)

  soil, state, θ_top, yobs
end


function Soil_main(config::Config, data_obs::AbstractMatrix{T}, site_name::AbstractString, dates::AbstractVector;
  maxn::Int=1000,
  method_retention=nothing,
  outdir=nothing, log_file=nothing,
  plot_fun=nothing, plot_initial=false, plot_kw...) where {T<:Real}

  ## parameter
  outdir = guess_outdir(config, outdir)
  io = open_log(config, log_file)

  !isempty(outdir) && mkpath(outdir)
  !isnothing(method_retention) && (config.method_retention = method_retention)
  method_retention = config.method_retention
  (; objective, plot_file) = config

  ## state and forcing 
  soil, state, θ_top, yobs = setup(config, data_obs)

  ## optimization
  lower, upper = SM_paramBound(soil)
  theta0 = SM_param2theta(soil)

  loss0 = Soil_goal(config, theta0, state, θ_top, yobs) |> x -> round(x, digits=4)
  log(io, "[$site_name] $method_retention/$(config.method_solve) $objective maxn=$maxn init=$loss0")

  t0 = time()
  theta_opt, feval, _ = sceua(theta -> Soil_goal(config, theta, state, θ_top, yobs),
    theta0, lower, upper; maxn)
  log(io, "done $(round(time()-t0, digits=1))s feval=$feval")

  serialize("$outdir/theta_$(site_name)_$(method_retention)", theta_opt)
  SM_UpdateParam!(soil, theta_opt)

  ## visualization
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

  best_cost = Soil_goal(config, theta_opt, state, θ_top, yobs)
  log(io, "best $objective=$(round(-best_cost, digits=6))")
  !isnothing(io) && close(io)

  return soil, theta_opt, best_cost
end
