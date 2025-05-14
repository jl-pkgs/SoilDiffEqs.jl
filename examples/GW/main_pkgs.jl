function plot_θ(soil)
  N = soil.N
  z = soil.z[1:N]
  zwt = soil.zwt

  plot(title="zwt = $zwt m", xlabel="θ (m³ m⁻³)", ylabel="z [m]", legend=:topright)
  plot!(θ, z, label="θ_init")
  plot!(soil.θ, z, label="θ_next")
end

function solve_ode(dt)
  soil = init_soil(; dt, zwt=-2.5)
  N = soil.N
  tspan = (0.0, dt)  # Time span for the simulation
  u0 = soil.θ[1:N]
  prob = ODEProblem(RichardsEquation_Zeng2009, u0, tspan, soil)
  sol = solve(prob, Tsit5(), saveat=200)
  @show soil.timestep

  soil.Q = cal_Q_Zeng2009!(soil, soil.θ)
  u = sol.u[end]
  soil.θ_prev[1:N] = soil.θ
  soil.θ = u
  soil
end

function solve_zeng(dt)
  soil = init_soil(; zwt=-2.5, dt)
  soil_moisture_Zeng2009(soil)
  soil
end


# @testset "soil_moisture_zeng2009" 
function init_soil(; zwt=-0.5, dt=3600)
  wa = 4000.0 # [mm]
  par = get_soilpar(1; method_retention="Campbell")
  param = SoilParam(N, par)

  soil = Soil(dz; θ=deepcopy(θ), zwt, wa, dt, 
    param,
    method_retention="Campbell")
  return soil
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
