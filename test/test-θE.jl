using SoilDifferentialEquations, Test
# using Plots
# gr(framestyle=:box)

# sand: 40%, clay: 40%, silt: 20%
function cal_θE(z₋ₕ, z₊ₕ, zwt, param; fun=cal_θE_campbell)
  N = length(z₋ₕ)
  θE = zeros(N)

  for i = 1:N
    z1 = z₊ₕ[i]
    z0 = z₋ₕ[i]
    θE[i] = fun(z1, z0, zwt, param...)
  end
  θE
end

## param from zeng2009
# param_campbell = (;
#   ψ_sat=-22.7, # cm
#   θ_sat=0.44,
#   B=9.3
# )

@testset "cal_θE" begin
  Δz = fill(10., 50) # 10cm x 50 layers = 3m
  soil = soil_depth_init(Δz) # 向下为负
  (; z₋ₕ, z₊ₕ) = soil

  ## param from Bonan 2019, Table 8.3, sandy clay
  param_campbell = (;
    ψ_sat=-15.3, # cm
    θ_sat=0.426,
    B=10.4
  )
  param_van1980 = (;
    ψ_sat=-15.3, # cm
    # θ_sat=0.38,
    θ_sat=0.426,
    # θ_res=0.325,
    θ_res=0.15,
    α=0.027,
    n=1.23,
  )
  zwt = -3 * 1e2 # cm
  θE_campbell = cal_θE(z₋ₕ, z₊ₕ, zwt, param_campbell; fun=cal_θE_campbell)
  θE_van1980 = cal_θE(z₋ₕ, z₊ₕ, zwt, param_van1980; fun=cal_θE_van1980)

  @test GOF(θE_van1980, θE_campbell).R2 >= 0.995
  # plot(ylabel="Depth (cm)", xlabel="θE")
  # plot!(θE_campbell, soil.z, label="Campbell θE", linewidth=1)
  # plot!(θE_van1980, soil.z, label="Van1980 θE", linewidth=1)
  # hline!([zwt], label="zwt", linewidth=0.4)
end
