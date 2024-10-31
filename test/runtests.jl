using SoilDifferentialEquations, OrdinaryDiffEq, Test
# using RTableTools, Dates

@testset "soil_depth" begin
  z = -[1.25, 5, 10, 20, 50, 100.0]
  dz = cal_Δz(z)
  z2, z₊ₕ, dz₊ₕ = soil_depth_init(dz)
  @test z == z2
end

include("test-soil_temperature.jl")
include("test-soil_temperature_F0.jl")


include("test-soil_moisture.jl")
include("test-soil_moisture_Q0.jl")

include("test-solve_SM.jl")
include("test-solve_Tsoil.jl")
