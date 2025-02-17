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
function init_soil(; zwt=-0.5)
  wa = 4000.0 # [mm]
  dt = 3600 # [s]
  N = 100
  dz = fill(0.02, N) # 2m
  θ = fill(0.3, N)
  soil = Soil(dz; θ=deepcopy(θ), zwt, wa, dt)
  return soil
end


begin
  soil = init_soil(; zwt=-2.5)
  θ_zeng2009 = soil_moisture_Zeng2009(soil)
  # @test maximum(soil.θ) ≈ 0.30996934526428166
  # p1 = plot_θ(soil)
end


begin
  p = init_soil(; zwt=-2.5)
  tspan = (0.0, 3600)  # Time span for the simulation
  u0 = p.θ[1:N]

  prob = ODEProblem(RichardsEquation_Zeng2009, u0, tspan, p)
  sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=200)
  @show p.timestep
  θ_ode = sol.u[end]
end

begin
  (; N) = soil
  plot(θ_ode, soil.z[1:N], label="θ_ode")
  plot!(θ_zeng2009, soil.z[1:N], label="θ_zeng2009")
end
