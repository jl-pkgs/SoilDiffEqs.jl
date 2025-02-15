using SoilDifferentialEquations, Test
# using Plots

@testset "soil_moisture_zeng2009" begin
  wa = 4000.0 # [mm]
  zwt = -0.5
  Δt = 60 # [s]
  # recharge = 1 / 3600 # [mm s-1], namely [1 mm h-1]

  N = 100
  dz = fill(0.02, N) # 2m
  θ = fill(0.3, N)
  # z₊ₕ = soil.z₊ₕ
  soil = Soil(dz; θ=deepcopy(θ), zwt=-0.5, wa)
  soil_moisture_zeng2009(soil)
  @test maximum(soil.θ) ≈ 0.30996934526428166
  # plot(soil.θ, soil.z[1:end])

  soil = Soil(dz; θ, zwt=-5.0, wa)
  @time soil_moisture_zeng2009(soil)
  @test maximum(soil.θ) ≈ 0.3043210817899786
end
