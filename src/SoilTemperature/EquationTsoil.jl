export soil_HeatFlux!

"""
计算每层的土壤的热通量 [W m-2]

- `ibeg`: 第一层土壤温度观测具有较大的误差，因此不使用第一层土壤温度
"""
function soil_HeatFlux!(soil::Soil{FT}, T::AbstractVector{FT};
  F0::FT=NaN, Tsurf::FT=NaN, method="Tsurf", ibeg::Int=1) where {FT<:Real}
  (; N, F, z, z₊ₕ) = soil
  (; κ) = soil.param
  
  if method == "Tsurf"
    if ibeg > 1
      d1 = z[ibeg-1] - z[ibeg]
      d2 = z[ibeg] - z[ibeg+1]
      _κ₊ₕ = mean_arithmetic(κ[ibeg-1], κ[ibeg], d1, d2)
      _dz = z[ibeg-1] - z[ibeg]
    else
      _κ₊ₕ = κ[1]
      _dz = 0 - z[1]
    end
    F0 = -_κ₊ₕ * (Tsurf - T[ibeg]) / _dz
  end

  @inbounds for i in ibeg:N-1
    d1 = z[i] - z₊ₕ[i]
    d2 = z₊ₕ[i] - z[i+1]
    κ₊ₕ = mean_arithmetic(κ[i], κ[i+1], d1, d2)  # Eq. 5.16, 
    Δz₊ₕ = z[i] - z[i+1]
    F[i] = -κ₊ₕ * (T[i] - T[i+1]) / Δz₊ₕ
  end
  F[N] = 0
  return F0
end


"""
# Arguments
- `method`: 
  + `Tsurf` : Tsurf boundary condition, 第一类边界条件
  + `F0`  : F0 boundary condition, 第二类边界条件
"""
function TsoilEquation(dT, T, soil::Soil, t; method="Tsurf", ibeg::Int=1)
  soil.timestep += 1
  # TODO: 根据t，更新Tsurf
  (; N, Δz, F, Tsurf) = soil
  (; cv) = soil.param
  F0 = soil_HeatFlux!(soil, T; Tsurf, F0=soil.F0, method, ibeg)

  dT[ibeg] = -(F0 - F[ibeg]) / (Δz[ibeg] * cv[ibeg])
  @inbounds for i in ibeg+1:N
    dT[i] = -(F[i-1] - F[i]) / (Δz[i] * cv[i])
  end
end

# ibeg:N
function TsoilEquation_partial(dT, T, p::Soil, t; method="Tsurf", ibeg::Int=1)
  (; N) = p
  p.du[ibeg:N] .= dT
  p.u[ibeg:N] .= T
  TsoilEquation(p.du, p.u, p, t; method, ibeg)
  dT .= p.du[ibeg:N]
  return nothing
end
