using SoilDifferentialEquations, Test
using Plots
using OrdinaryDiffEqTsit5
using Plots
gr(framestyle=:box, legend=:topright)

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
  Δ = 0.05
  N = floor(Int, 2 / Δ)
  dz = fill(Δ, N) # 2m
  θ = fill(0.3, N)
  soil = Soil(dz; θ=deepcopy(θ), zwt, wa, dt)
  return soil
end


begin
  dt = 360
  soil_ode = solve_ode(dt)
  soil_zeng = solve_zeng(dt)
  N = soil_ode.N

  error_SM(soil_ode)  |> display
  error_SM(soil_zeng) |> display

  θ_zeng = soil_zeng.θ
  θ_ode = soil_ode.θ

  plot(θ_ode, soil_ode.z[1:N], label="θ_ode")
  plot!(θ_zeng, soil_zeng.z[1:N], label="θ_zeng")
end
