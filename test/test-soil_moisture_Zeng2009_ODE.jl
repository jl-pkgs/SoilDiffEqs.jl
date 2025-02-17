using SoilDifferentialEquations, Test
using Plots
using OrdinaryDiffEqTsit5

# @testset "soil_moisture_zeng2009" 
function init_soil(; zwt=-0.5)
  wa = 4000.0 # [mm]
  dt = 120 # [s]

  N = 100
  dz = fill(0.02, N) # 2m
  θ = fill(0.3, N)
  soil = Soil(dz; θ=deepcopy(θ), zwt, wa, dt)
  return soil
end

soil = init_soil()
soil_moisture_Zeng2009(soil)


begin
  soil = init_soil(; zwt=-2.5)
  soil_moisture_Zeng2009(soil)
  # @test maximum(soil.θ) ≈ 0.30996934526428166
  p1 = plot_θ(soil)
end



# function solve_ode()
begin
  p = init_soil()
  tspan = (0.0, 3600)  # Time span for the simulation
  u0 = p.θ[1:N]

  prob = ODEProblem(RichardsEquation_Zeng2009, u0, tspan, p)
  sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=200)
  @show p.timestep
  sol.u[end]
end


begin
  (; N) = soil
  plot(soil.θ, soil.z[1:N])
end


begin
  soil_moisture_Zeng2009(soil)
  # @test maximum(soil.θ) ≈ 0.30996934526428166
  # plot(soil.θ, soil.z[1:end])
  # soil = Soil(dz; θ, zwt=-5.0, wa)
  # @time soil_moisture_Zeng2009(soil)
  # @test maximum(soil.θ) ≈ 0.3043210817899786
end
