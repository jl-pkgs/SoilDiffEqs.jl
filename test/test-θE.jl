using SoilDifferentialEquations, Test

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
  N = 50
  Δz = fill(10., N) # 10cm x 50 layers = 3m
  soil = soil_depth_init(Δz) # 向下为负
  (; z₋ₕ, z₊ₕ) = soil

  ## param from Bonan 2019, Table 8.3, sandy clay
  param_campbell = Campbell(;
    ψ_sat=-15.3, # cm
    θ_sat=0.426,
    b=10.4
  )
  ψ_sat = -15.3 # cm
  param_van1980 = VanGenuchten(;
    # θ_sat=0.38,
    # θ_res=0.325,
    θ_sat=0.426,
    θ_res=0.15,
    α=0.027,
    n=1.23,
  )
  zwt = -2.0 * 1e2 # cm
  θE_campbell = main_θE(z₋ₕ, z₊ₕ, zwt, ψ_sat, param_campbell)
  θE_van1980 = main_θE(z₋ₕ, z₊ₕ, zwt, ψ_sat, param_van1980)

  @test GOF(θE_van1980, θE_campbell).R2 >= 0.995
end


# begin
#   zwt = -4.0 * 1e2 # cm
#   θE_campbell_4 = main_θE(z₋ₕ, z₊ₕ, zwt, ψ_sat, param_campbell)
# end


# # 这里最后一层是地下水
# begin
#   using Plots
#   gr(framestyle=:box)

#   plot(ylabel="Depth (cm)", xlabel="θE")
#   inds = 1:N
#   plot!(θE_campbell[inds], soil.z[inds], label="θE (zwt = -2m)", linewidth=1)
#   plot!(θE_campbell_4[inds], soil.z[inds], label="θE (zwt = -4m)", linewidth=1)

#   # plot!(θE_van1980[inds], soil.z[inds], label="Van1980 θE", linewidth=1)
#   hline!([-200.], label="", linewidth=0.4, color=:blue, linestyle=:dash)
#   hline!([-400.], label="", linewidth=0.4, color=:red, linestyle=:dash)
#   savefig("Figure_θE.pdf")
# end
