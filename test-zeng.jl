# using Pkg
# Pkg.activate(".")
using SoilDifferentialEquations, Test
using Plots

begin
  wa = 4000.0 # [mm]
  zwt = -0.5
  Δt = 60 # [s]
  # recharge = 1 / 3600 # [mm s-1], namely [1 mm h-1]

  N = 100
  dz = fill(0.02, N) # 2m
  θ = fill(0.3, N)
  soil = Soil(dz; θ, zwt, wa)
  # z₊ₕ = soil.z₊ₕ
  soil_moisture_zeng2009(soil)
  plot(soil.θ, soil.z[1:end])
end


begin
  wa = 4000.0 # [mm]
  zwt = -5.0
  Δt = 60 # [s]
  # recharge = 1 / 3600 # [mm s-1], namely [1 mm h-1]

  N = 100
  dz = fill(0.02, N) # 2m
  θ = fill(0.3, N)
  soil = Soil(dz; θ, zwt, wa)
  # z₊ₕ = soil.z₊ₕ
  @run soil_moisture_zeng2009(soil)
  plot(soil.θ, soil.z[1:end])
end

