using SoilDifferentialEquations, Test


function _init_soil(dz; zwt=-0.5, wa=4000.0)
  N = length(dz) # 2.1m
  θ = fill(0.3, N)
  Soil(dz; θ, zwt, wa)
end

begin
  wa = 4000.0 # [mm]
  zwt = -0.5
  Δt = 60 # [s]
  # recharge = 1 # [mm h-1]
  # dz = fill(0.02, 100)
  dz = [fill(0.1, 5); fill(0.2, 3); fill(0.5, 2)]
  N = length(dz)
  soil = _init_soil(dz)
  z₊ₕ = soil.z₊ₕ
end


@testset "GW_Update_ZWT!" begin
  z₊ₕ = soil.z₊ₕ
  zwt = soil.zwt

  ## recharge
  θ = fill(0.3, N)
  r = GW_Update_ZWT!(soil, θ; zwt=-2.0, wa=4000, ∑=105)
  @test r.zwt == -1.05
  @test θ[8] == 0.32499999999999996

  θ = fill(0.3, N)
  r = GW_Update_ZWT!(soil, θ; zwt=-2.5, wa=4000, ∑=105)
  @test r.zwt ≈ -1.25
  @test θ[9] == 0.36999999999999994

  ## discharge
  θ = fill(0.3, N)
  r = GW_Update_ZWT!(soil, θ; zwt=-2.5, wa=4000, ∑=-105)
  @test r.zwt ≈ -3.025

  θ = fill(0.3, N)
  r = GW_Update_ZWT!(soil, θ; zwt=-2.0, wa=4000, ∑=-105)
  @test r.zwt ≈ -2.525
  @test θ[10] == 0.2600000000000001
end

begin
  zwt = -2.5
  specific_yield!(soil, zwt)
  (; Sy_d, Sy_r, Sy_e) = soil
  d = (; Sy_d, Sy_r, Sy_e)
end
