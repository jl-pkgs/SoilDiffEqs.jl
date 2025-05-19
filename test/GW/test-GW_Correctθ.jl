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
  recharge = 1 # [mm h-1]
  # dz = fill(0.02, 100)
end

@testset "cal_θEψE!" begin
  dz = fill(0.02, 100)
  soil = _init_soil(dz)

  soil.zwt = -0.5
  ψE = cal_θEψE!(soil)
  @test ψE[end] == 0.0
  @test -50 <= minimum(ψE) <= -49

  soil.zwt = -5.0
  ψE = cal_θEψE!(soil)
  @test minimum(ψE) >= -500
  @test maximum(ψE) <= -130
  soil.zwt = -0.5
  # soil.zwt = -2.5
  # ψE = cal_θEψE!(soil)
  # gr(framestyle=:box)
  # plot(
  #   plot(soil.θE[1:N], soil.z[1:N], ylabel="depth", xlabel="θE (%)", label=""),
  #   plot(soil.ψE[1:N], soil.z[1:N], ylabel="depth", xlabel="ψE (cm)", label="")
  # )
end

@testset "find_jwt" begin
  dz = fill(0.02, 100)
  N = length(dz)
  soil = _init_soil(dz)
  z₊ₕ = soil.z₊ₕ

  @test find_jwt(z₊ₕ, 0.0; N) == 1
  @test find_jwt(z₊ₕ, -0.01; N) == 1
  @test find_jwt(z₊ₕ, -0.02; N) == 1 # 调整了边界计数
  @test find_jwt(z₊ₕ, -0.03; N) == 2
  @test find_jwt(z₊ₕ, -100.0; N) == N + 1
end

begin
  dz = [fill(0.1, 5); fill(0.2, 3); fill(0.5, 2)] # 2.1m
  N = length(dz)
  soil = _init_soil(dz)
  z₊ₕ = soil.z₊ₕ

  zwt = -1.5
  jwt = find_jwt(z₊ₕ, zwt; N)
  GW = fill(0.0, N)
  GW[jwt] = 0.5
  GW[jwt+1:N] .= 1.0

  Δt = 1 # [h]
  drainage = 0  # [cm h-1]
  wa = 4000 # [mm]

  θ = fill(0.3, N)
  r = GW_Correctθ!(soil, θ; zwt, exceed2surf=true) # 
  @test r.∑ ≈ -50.0

  θ = fill(0.42, N) # 42mm, 2.1m * 0.02
  r = GW_Correctθ!(soil, θ; zwt, exceed2surf=false) 
  @test r.∑ ≈ 42.0

  θ = fill(0.42, N) # 42mm, 2.1m * 0.02
  r = GW_Correctθ!(soil, θ; zwt, exceed2surf=true) # 
  @test r.uex ≈ 42.0
end
