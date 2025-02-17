export find_jwt


function find_jwt(z₊ₕ::AbstractVector, zwt::Real; N::Int=length(z₊ₕ))
  jwt = N
  zwt_abs = abs(zwt)
  for j in 1:N
    if zwt_abs <= abs(z₊ₕ[j])
      jwt = j - 1
      break
    end
  end
  return jwt
end

# include("Soil_depth.jl")
include("soil_ParamTable.jl")
include("Retention_Campbell.jl")
include("Retention_van_Genuchten.jl")
include("Retention.jl")

include("Equilibrium.jl")
include("EquationRichards.jl")
include("Equation_Zeng2009.jl")
include("soil_moisture_Bonan.jl")
include("soil_moisture_Bonan_Q0.jl")
include("soil_moisture_BEPS.jl")
include("soil_moisture_Zeng2009.jl")
include("Solve_SM.jl")
