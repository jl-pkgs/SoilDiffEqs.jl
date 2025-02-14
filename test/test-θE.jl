using SoilDifferentialEquations, Test
# using Plots
# gr(framestyle=:box)

function main_θE(z₋ₕ, z₊ₕ, zwt, ψ_sat, par)
  N = length(z₋ₕ)
  θE = zeros(N)
  for i = 1:N
    z1 = z₊ₕ[i]
    z0 = z₋ₕ[i]
    θE[i] = cal_θE(z1, z0, zwt, ψ_sat, par)
  end
  θE
end

## param from zeng2009, sand: 40%, clay: 40%, silt: 20%
# param_campbell = (;
#   ψ_sat=-22.7, # cm
#   θ_sat=0.44,
#   B=9.3
# )

# @testset "cal_θE" 
begin
  Δz = fill(10., 50) # 10cm x 50 layers = 3m
  soil = soil_depth_init(Δz) # 向下为负
  (; z₋ₕ, z₊ₕ) = soil

  ## param from Bonan 2019, Table 8.3, sandy clay
  param_campbell = ParamCampbell(;
    ψ_sat=-15.3, # cm
    θ_sat=0.426,
    b=10.4
  )
  ψ_sat = -15.3 # cm
  param_van1980 = ParamVanGenuchten(;
    # θ_sat=0.38,
    θ_sat=0.426,
    # θ_res=0.325,
    θ_res=0.15,
    α=0.027,
    n=1.23,
  )
  zwt = -3 * 1e2 # cm
  θE_campbell = main_θE(z₋ₕ, z₊ₕ, zwt, ψ_sat, param_campbell)
  θE_van1980 = main_θE(z₋ₕ, z₊ₕ, zwt, ψ_sat, param_van1980)

  @test GOF(θE_van1980, θE_campbell).R2 >= 0.995
  # plot(ylabel="Depth (cm)", xlabel="θE")
  # plot!(θE_campbell, soil.z[1:end], label="Campbell θE", linewidth=1)
  # plot!(θE_van1980, soil.z[1:end], label="Van1980 θE", linewidth=1)
  # hline!([zwt], label="zwt", linewidth=0.4)
end
