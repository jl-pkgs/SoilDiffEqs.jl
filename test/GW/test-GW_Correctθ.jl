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
  @test find_jwt(z₊ₕ, -0.02; N) == 2
  @test find_jwt(z₊ₕ, -0.03; N) == 2
  @test find_jwt(z₊ₕ, -100.0; N) == N+1
end


@testset "GW_Correctθ!" begin
  # 亏损
  dz = fill(0.02, 100)
  soil = _init_soil(dz)
  N = length(dz)
  θ = fill(0.3, N)
  
  Δt = 1/60 # [h]
  drainage = 60  # [cm h-1]
  wa = 4000 # [mm]

  @test GW_Correctθ!(soil, θ, -0.5, wa, Δt, drainage) == (wa=wa, uex=0.0, drainage=60)

  θ[2:3] .= -0.1
  @test GW_Correctθ!(soil, θ, -0.5, wa, Δt, drainage) == (wa=wa, uex=0.0, drainage=36)

  θ[2:3] .= -10.0 # [m3 m-3]
  @test GW_Correctθ!(soil, θ, -0.5, wa, Δt, drainage) == (wa=3610.0, uex=0.0, drainage=0.0)
  
  # 超饱和
  θ = fill(0.3, N)
  θ[2:3] .= 0.6
  @test GW_Correctθ!(soil, θ, -0.5, wa, Δt, drainage) == (wa=wa, uex=5.999999999999998, drainage=60)

  θ = fill(0.3, N)
  θ[100] = 12.0
  r = GW_Correctθ!(soil, θ, -0.5, wa, Δt, drainage) 
  @test all(θ .== 0.4)
  @test r == (wa=wa, uex=33.999999999999986, drainage=60)
end

## 限制条件
# 1. zwt以下，地下水饱和
