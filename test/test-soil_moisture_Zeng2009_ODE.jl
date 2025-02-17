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

# @testset "soil_moisture_zeng2009" 
function init_soil(; zwt=-0.5, dt=3600)
  wa = 4000.0 # [mm]
  N = 100
  dz = fill(0.02, N) # 2m
  θ = fill(0.3, N)
  soil = Soil(dz; θ=deepcopy(θ), zwt, wa, dt)
  return soil
end

function solve_ode(dt)
  soil = init_soil(; dt, zwt=-2.5)
  tspan = (0.0, dt)  # Time span for the simulation
  u0 = soil.θ[1:N]
  prob = ODEProblem(RichardsEquation_Zeng2009, u0, tspan, soil)
  sol = solve(prob, Tsit5(), saveat=200)
  @show soil.timestep

  u = sol.u[end]
  soil.θ_prev[1:N] = soil.θ
  soil.θ = u
  soil.Q = cal_Q_Zeng2009!(soil, θ)
  soil
end

function solve_zeng(dt)
  soil = init_soil(; zwt=-2.5, dt)
  soil_moisture_Zeng2009(soil)
  soil
end

dt = 360
soil_ode = solve_ode(dt)
soil_zeng = solve_zeng(dt)

error_SM(soil_ode)  # 
error_SM(soil_zeng) # zeng的误差更小

begin
  dt = 60
  θ_zeng2009 = soil.θ 
  plot(θ_ode, soil.z[1:N], label="θ_ode")
  plot!(θ_zeng2009, soil.z[1:N], label="θ_zeng2009")
end
