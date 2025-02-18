export find_jwt
export error_SM

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



function error_SM(soil::Soil{T}) where {T<:Real}
  (; N, Q0) = soil
  θ_prev = soil.θ_prev
  θ_next = soil.θ
  dt = soil.dt / 3600 # [s] -> [h]
  Δz = soil.Δz_cm
  
  Q_prev = cal_Q_Zeng2009!(soil, θ_prev)
  Q_next = cal_Q_Zeng2009!(soil, θ_next)

  Q = zeros(N)
  θ = zeros(N)

  # 中间插的思路（最精确的解法），计算误差
  dθ = 0.0
  for i in 1:N
    Q[i] = 0.5 * (Q_prev[i] + Q_next[i])
    θ[i] = 0.5 * (θ_prev[i] + θ_next[i])
    dθ += (θ_next[i] - θ_prev[i]) * Δz[i] # in cm
  end
  QN = Q[N]
  
  obs = (QN - Q0) * dt # cm
  sim = dθ

  bias = sim - obs
  perc = bias / obs * 100 |> _round

  info = (; obs, sim, bias, perc, dθ, QN, Q0)
  info
end

_round(x::Real; digits=3) = round(x; digits)
