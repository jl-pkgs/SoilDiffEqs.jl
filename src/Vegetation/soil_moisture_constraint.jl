# Function to compute soil water stress factor
function soil_water_factor(soil::Soil, args...; fun=cal_β_BEPS, kw...)
  N = soil.N
  for i in 1:N
    soil.βw[i] = fun(soil.ψ[i] * 0.01, args...; kw...)  # 土壤水限制, [cm] -> [m]
    # soil.βt[i] = cal_β_BEPS(soil.Tsoil[i])           # 土壤温度限制
    soil.w[i] = soil.f_root[i] * soil.βw[i]
  end
  ∑ = sum(soil.w) # 每层的土壤水分限制因子

  β = 0.0
  for i in 1:N
    # isnan(soil.dt[i]) && println(soil.dt[1])
    soil.w[i] = soil.w[i] / ∑
    β += soil.βw[i] * soil.w[i]
  end
  ∑ < 1e-6 && (β = 0.1)
  soil.β = max(0.1, β)
  return β
end

# if landcover == 6 || landcover == 9 # EBF or DBF
#   p.ψ_min = 10.0 # ψ_min
#   p.alpha = 1.5
# else
#   p.ψ_min = 33.0 # ψ_min
#   p.alpha = 0.4
# end
"""
- `ψ`: 越干越负, 饱和为0, in [m] or [kPa]
"""
function cal_βw_BEPS(ψ::T; ψ_min=33.0, α=0.4)::T where {T<:Real}
  abs(ψ) <= ψ_min && return T(1)
  1 / (1 + (abs(ψ) / ψ_min - 1)^α) ## He 2017, JGR-B, Eq.4
end

# ρ g h = Pa
# - `33kPa`  : field capacity, 33 m
# - `1500kPa`: wilting point, 1500m
function cal_βt_BEPS(Tsoil::T; t1=-0.02, t2=2.0)::T where {T<:Real}
  Tsoil <= 0 && return T(0)
  1 - exp(t1 * Tsoil^t2) ## He 2017, JGR-B, Eq.4
end


export cal_βt_BEPS, cal_βw_BEPS
export soil_water_factor
