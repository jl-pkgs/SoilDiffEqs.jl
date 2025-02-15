using SoilDifferentialEquations, Test
@time using OrdinaryDiffEqTsit5
# using RTableTools, Dates

@testset "soil_depth" begin
  z = -[1.25, 5, 10, 20, 50, 100.0]
  dz = cal_Δz(z)
  z2, z₋ₕ, z₊ₕ, dz₊ₕ = soil_depth_init(dz)
  @test z2[1:6] == z
end

include("test-soil_moisture_Zeng2009.jl")
include("test-θE.jl")
include("test-GW.jl")

include("test-soil_temperature.jl")
include("test-soil_temperature_F0.jl")

include("test-soil_moisture.jl")
include("test-soil_moisture_Q0.jl")

include("test-solve_Tsoil.jl")
include("test-solve_SM.jl")

include("SM_uscrn/case_solve_BEPS.jl")
include("SM_uscrn/case_solve_Bonan.jl")
