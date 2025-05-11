using SoilDifferentialEquations, Test

using Plots
using OrdinaryDiffEqTsit5
using Plots
gr(framestyle=:box, legend=:topright)

includet("main_pkgs.jl")

begin
  Δ = 0.02
  dt = 3600
  N = floor(Int, 2 / Δ)
  dz = fill(Δ, N) # 2m
  # θ = fill(0.3, N)
  θ = LinRange(0.36, 0.2, N) |> collect

  ## 求解的均是1个时刻
  @time soil_ode = solve_ode(dt)
  # soil_zeng = solve_zeng(dt)
  error_SM(soil_ode) |> display
  ∂K₊ₕ∂θ = _cal_K(soil_ode)

  soil_zeng = init_soil(; zwt=-2.5, dt);
  (; ∂K₊ₕ∂θ, ∂ψ∂θ, ∂qᵢ∂θᵢ, ∂qᵢ∂θᵢ₊₁) = soil_moisture_Zeng2009(soil_zeng; ∂K₊ₕ∂θ)
  error_SM(soil_zeng) |> display
  "ok"
  # θ_zeng = soil_zeng.θ
  # θ_ode = soil_ode.θ
  # plot(θ_ode, soil_ode.z[1:N], label="θ_ode", legend=:bottomright)
  # plot!(θ_zeng, soil_zeng.z[1:N], label="θ_zeng")
end



function _cal_K(soil)
  θ_prev = soil.θ_prev
  θ_next = soil.θ
  (; N) = soil
  cal_K_CLM5!(soil, soil.θ_prev)
  K₊ₕ_prev = deepcopy(soil.K₊ₕ)

  cal_K_CLM5!(soil, soil.θ)
  K₊ₕ_next = deepcopy(soil.K₊ₕ)

  ∂K₊ₕ∂θ = zeros(N)
  for i in 1:N
    dθ = θ_prev[i] - θ_next[i]
    dK = K₊ₕ_prev[i] - K₊ₕ_next[i]
    ∂K₊ₕ∂θ[i] = dK / dθ
  end
  ∂K₊ₕ∂θ
end



begin
  soil = init_soil(; zwt=-2.5, dt)
  (; ∂K₊ₕ∂θ, ∂ψ∂θ, ∂qᵢ∂θᵢ, ∂qᵢ∂θᵢ₊₁) = soil_moisture_Zeng2009(soil)

  
  plot(obs_∂K₊ₕ∂θ, soil.z[1:N], label="obs_∂K₊ₕ∂θ", legend=:bottomright)
  plot!(∂K₊ₕ∂θ, soil.z[1:N], label="sim_∂K₊ₕ∂θ")
end
