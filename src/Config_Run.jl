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
      θ_sat = soil_params["theta_sat"],
      θ_res = soil_params["theta_res"],
      Ksat  = soil_params["ksat"],
      α     = soil_params["alpha"],
      n     = soil_params["n"],
      m     = 1.0 - 1.0/soil_params["n"]
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


export InitSoil, SM_simulate, SM_goal
