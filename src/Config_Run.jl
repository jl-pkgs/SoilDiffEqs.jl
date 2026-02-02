function InitSoil(config::Config, data_obs::AbstractMatrix{T}) where {T<:Real}
  (; zs_obs, z_bound_top,
    zs_center, dt,
    soil_type, method_retention, same_layer, soil_params) = config

  z₊ₕ = center_to_face(zs_center)
  Δz = face_to_thickness(z₊ₕ) ./ 100.0 # [cm] -> [m]
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  N = length(Δz)
  zs_sim = z[1:N] .* 100 # [m] -> [cm]

  ibeg = findfirst(==(z_bound_top), abs.(zs_sim)) + 1 # index of soil model start (simulation layer)
  itop = findfirst(==(z_bound_top), abs.(zs_obs))     # index of obs top boundary layer  

  θ_surf = data_obs[:, itop]        # 边界层（ibeg-1层）的数据作为地表输入
  yobs = data_obs[:, (itop+1):end]  # 从模拟起始层（ibeg层）开始的观测数据 

  θ0 = data_obs[1, :]
  θ = fill(0.2, N)
  n = min(length(θ0), N)
  θ[1:n] .= θ0[1:n]
  θ[n+1:N] .= θ0[end]

  # 使用自定义参数或标准参数
  if soil_params !== nothing && method_retention == "van_Genuchten"
    par = VanGenuchten{Float64}(
      θ_sat=soil_params["theta_sat"],
      θ_res=soil_params["theta_res"],
      Ksat=soil_params["ksat"],
      α=soil_params["alpha"],
      n=soil_params["n"],
      m=1.0 - 1.0 / soil_params["n"]
    )
  else
    par = get_soilpar(soil_type; method_retention)
  end

  param = SoilParam(N, par; use_m=false, method_retention, same_layer)
  soil = Soil{Float64}(; N, ibeg, dt, z, z₊ₕ, Δz, Δz₊ₕ, θ, method_retention, param)
  cal_ψ!(soil)
  cal_K!(soil)
  soil, θ_surf, yobs
end

# Model functions
function SM_simulate(config::Config, data_obs::AbstractMatrix{T}, theta) where {T<:Real}
  soil, θ_surf, yobs = InitSoil(config, data_obs)
  SM_UpdateParam!(soil, theta)

  if config.method_solve == "Bonan"
    ysim = solve_SM_Bonan(soil, θ_surf)
  elseif config.method_solve == "ODE"
    ysim = solve_SM_ODE(soil, θ_surf; solver=Tsit5())
  elseif config.method_solve == "BEPS"
    ysim = soil_moisture_BEPS(soil, θ_surf)
  else
    error("Unknown method_solve: $(config.method_solve)")
  end
  ysim, yobs
end

function SM_goal(config::Config, data_obs::AbstractMatrix{T}, theta) where {T<:Real}
  of_fun = config.of_fun
  ysim, yobs = SM_simulate(config, data_obs, theta)
  ncol = size(yobs, 2)
  loss = 0.0
  for i in 1:ncol
    obs = view(yobs, :, i)
    sim = view(ysim, :, i)
    loss += -of_fun(obs, sim)
  end
  loss / ncol
end



function log(io, msg)
  println(msg)
  if !isnothing(io)
    println(io, msg)
    flush(io)
  end
end

function open_log(config::Config, log_file=nothing)
  (; fileConfig) = config
  isnothing(log_file) && (log_file = replace(fileConfig, r"\.yaml$" => ".log"))
  io = open(log_file, "a")
  return io
end

function SM_main(config::Config, data_obs::AbstractMatrix{T}, site_name::AbstractString, dates::AbstractVector;
  method_retention=nothing, outdir="", log_file=nothing,
  plot_fun=nothing, plot_initial=false, plot_kw...) where {T<:Real}

  !isempty(outdir) && mkpath(outdir)
  !isnothing(method_retention) && (config.method_retention = method_retention)

  method_retention = config.method_retention
  (; maxn, objective, plot_file) = config

  soil, θ_surf, yobs = InitSoil(config, data_obs)
  lower, upper = SM_paramBound(soil)
  theta0 = SM_param2theta(soil)

  loss0 = round(SM_goal(config, data_obs, theta0), digits=4)

  io = open_log(config, log_file)
  log(io, "[$site_name] $method_retention/$(config.method_solve) $objective maxn=$maxn init=$loss0")

  t0 = time()
  theta_opt, feval, _ = sceua(theta -> SM_goal(config, data_obs, theta), theta0, lower, upper; maxn)
  log(io, "done $(round(time()-t0, digits=1))s feval=$feval")

  Serialization.serialize(joinpath(outdir, "theta_$(site_name)_$(method_retention)"), theta_opt)
  SM_UpdateParam!(soil, theta_opt)

  if !isempty(plot_file) && !isnothing(plot_fun)
    depths = round.(Int, -soil.z[soil.inds_obs] .* 100)

    fout = joinpath(outdir, "$(method_retention)_$(plot_file)")
    ysim, yobs = SM_simulate(config, data_obs, SM_param2theta(soil))
    plot_fun(; ysim, yobs, dates, depths, fout, plot_kw...)

    if plot_initial
      base, ext = splitext(plot_file)
      initial_plot_file = "$(base)_initial$(ext)"

      fout = joinpath(outdir, "$(method_retention)_$(initial_plot_file)")
      ysim0, _ = SM_simulate(config, data_obs, theta0)
      plot_fun(; ysim=ysim0, yobs, dates, depths, fout, plot_kw...)
    end
  end

  best_cost = SM_goal(config, data_obs, SM_param2theta(soil))

  log(io, "best $objective=$(round(-best_cost, digits=6))")
  !isnothing(io) && close(io)
  return soil, theta_opt, best_cost
end


export InitSoil, SM_simulate, SM_goal, SM_main
