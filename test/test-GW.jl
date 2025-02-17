using SoilDifferentialEquations, Test

begin
  wa = 4000.0 # [mm]
  zwt = -0.5
  Δt = 60 # [s]
  recharge = 1 / 3600 # [mm s-1], namely [1 mm h-1]

  N = 100
  dz = fill(0.02, N) # 2m
  θ = fill(0.3, N)
  soil = Soil(dz; θ, zwt, wa)
  z₊ₕ = soil.z₊ₕ
end

@testset "cal_θEψE!" begin
  soil.zwt = -0.5
  ψE = cal_θEψE!(soil)
  @test ψE[end] == 0.0
  @test -50 <= minimum(ψE) <= -49

  soil.zwt = -5.0
  ψE = cal_θEψE!(soil)
  @test minimum(ψE) >= -500
  @test maximum(ψE) <= -130
  soil.zwt = -0.5
end

# soil.zwt = -0.5
# ψE = cal_θEψE!(soil)
# gr(framestyle=:box)
# plot(
#   plot(soil.θE[1:end-1], soil.z[1:end], ylabel="depth", xlabel="θE (%)"),
#   plot(soil.ψE[1:end-1], soil.z[1:end], ylabel="depth", xlabel="ψE (cm)")
# )

@testset "find_jwt" begin
  @test find_jwt(z₊ₕ, 0.0) == 0
  @test find_jwt(z₊ₕ, -0.01) == 0
  @test find_jwt(z₊ₕ, -0.02) == 0
  @test find_jwt(z₊ₕ, -0.03) == 1
  @test find_jwt(z₊ₕ, -100.0; N) == N
end




@testset "GW_UpdateRecharge!" begin
  @test GW_UpdateRecharge!(soil, wa, -0.5, Δt, recharge) == # 1 mm h-1
        (zwt=-0.49916666666666665, wa=4000.016666666667, uex=0.0)

  @test GW_UpdateRecharge!(soil, wa, -0.5, Δt, 1000 / 3600) == # 1000 mm h-1
        (zwt=0.0, wa=4016.6666666666665, uex=6.666666666666663)

  @test GW_UpdateRecharge!(soil, wa, -0.5, Δt, -recharge) == # -1 mm h-1
        (zwt=-0.5008333333333334, wa=3999.983333333333, uex=0.0)

  @test GW_UpdateRecharge!(soil, wa, -0.5, Δt, -10000 / 3600) == # -10,000 mm h-1
        (zwt=-8.833333333333291, wa=3833.3333333333335, uex=0.0)
end


@testset "GW_UpdateDrainage!" begin
  θ = fill(0.3, N)
  jwt = find_jwt(z₊ₕ, -0.5)
  jwt2 = find_jwt(z₊ₕ, -1.0)
  r = GW_UpdateDrainage!(soil, θ, -0.5, 4000.0, Δt, 600 / 3600)

  @test r == (zwt=-1.0, wa=3990.0)
  @test all(θ[jwt+1:jwt2+1] .< 0.3)

  @test GW_UpdateDrainage!(soil, θ, -2.5, 4000.0, Δt, 1 / 3600) ==
        (zwt=-2.5008333333333335, wa=3999.983333333333)
end


@testset "GW_Correctθ!" begin
  # 亏损
  θ = fill(0.3, N)
  @test GW_Correctθ!(soil, θ, -0.5, 4000.0, 60, 600 / 3600) == (wa=4000.0, uex=0.0, drainage=0.16666666666666666)

  θ[2:3] .= -0.1
  @test GW_Correctθ!(soil, θ, -0.5, 4000.0, 60, 600 / 3600) == (wa=4000.0, uex=0.0, drainage=0.09999999999999999)

  θ[2:3] .= -10.0 # [m3 m-3]
  @test GW_Correctθ!(soil, θ, -0.5, 4000.0, 60, 600 / 3600) == (wa=3610.0, uex=0.0, drainage=0.0)
  
  # 超饱和
  θ = fill(0.3, N)
  θ[2:3] .= 0.6
  @test GW_Correctθ!(soil, θ, -0.5, 4000.0, 60, 600 / 3600) == (wa=4000.0, uex=5.999999999999998, drainage=0.16666666666666666)

  θ = fill(0.3, N)
  θ[100] = 12.0
  r = GW_Correctθ!(soil, θ, -0.5, 4000.0, 60, 600 / 3600) 
  @test all(θ .== 0.4)
  @test r == (wa=4000.0, uex=33.999999999999986, drainage=0.16666666666666666)
end
