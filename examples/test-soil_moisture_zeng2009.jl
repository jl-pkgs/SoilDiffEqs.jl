using Plots
gr(framestyle=:box, legend=:topright)
function plot_θ(soil)
  zwt = soil.zwt
  plot(title="zwt = $zwt m", xlabel="θ (m³ m⁻³)", ylabel="z [m]", legend=:topright)
  plot!(θ, soil.z[1:end], label="θ_init")
  plot!(soil.θ, soil.z[1:end], label="θ_next")
end


begin
  wa = 4000.0 # [mm]
  zwt = -0.5
  Δt = 60 # [s]
  # recharge = 1 / 3600 # [mm s-1], namely [1 mm h-1]

  N = 100
  dz = fill(0.02, N) # 2m
  θ = fill(0.3, N)
  # z₊ₕ = soil.z₊ₕ
  soil = Soil(dz; θ=deepcopy(θ), zwt=-0.5, wa)
  soil_moisture_Zeng2009(soil)
  @test maximum(soil.θ) ≈ 0.30996934526428166
  # p1 = plot_θ(soil)

  soil = Soil(dz; θ=deepcopy(θ), zwt=-2.5, wa)
  @time soil_moisture_Zeng2009(soil)
  θ_next3 = soil.θ
  # plot_θ(soil)
  soil = Soil(dz; θ=deepcopy(θ), zwt=-3.0, wa)
  @time soil_moisture_Zeng2009(soil)
  θ_next5 = soil.θ

  plot(xlabel="θ (m³ m⁻³)", ylabel="Depth [m]", legend=:topright)
  plot!(θ, soil.z[1:end], label="θ_init")
  plot!(θ_next3, soil.z[1:end], label="θ_next (zwt = -2.5 m)")
  plot!(θ_next5, soil.z[1:end], label="θ_next (zwt = -3.0 m)")
end
