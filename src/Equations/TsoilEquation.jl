"""
计算每层的土壤的热通量 [W m-2]

- `ibeg`: 第一层土壤温度观测具有较大的误差，因此不使用第一层土壤温度
"""
function soil_HeatFlux!(F::V, T::V, κ::V, z::V, z₊ₕ::V;
  F0::FT=NaN, TS0::FT=NaN, method="TS0", ibeg::Int=1) where {FT<:Real,V<:AbstractVector{FT}}

  n = length(T)
  if method == "TS0"
    if ibeg == 1
      F0 = -κ[1] * (TS0 - T[1]) / (0 - z[1])
    elseif ibeg > 1
      F0 = -κ[ibeg] * (TS0 - T[ibeg]) / (z[ibeg-1] - z[ibeg])
    end
  end

  @inbounds for i in ibeg:n-1
    d1 = z[i] - z₊ₕ[i]
    d2 = z₊ₕ[i] - z[i+1]
    κ₊ₕ = κ[i] * κ[i+1] * (d1 + d2) / (κ[i] * d2 + κ[i+1] * d1) # Eq. 5.16, 
    # κ₊ₕ = (κ[i] + κ[i+1]) / 2
    Δz₊ₕ = z[i] - z[i+1]
    F[i] = -κ₊ₕ * (T[i] - T[i+1]) / Δz₊ₕ
  end
  F[n] = 0
  return F0
end


"""
# Arguments
- `method`: 
  + `TS0` : TS0 boundary condition, 第一类边界条件
  + `F0`  : F0 boundary condition, 第二类边界条件
"""
function TsoilEquation(dT, T, p::Soil, t; method="TS0", ibeg::Int=1)
  p.timestep += 1
  # TODO: 根据t，更新TS0
  (; n, Δz, z, z₊ₕ, F, κ, cv, TS0) = p
  F0 = soil_HeatFlux!(F, T, κ, z, z₊ₕ; TS0, F0=p.F0, method, ibeg)

  dT[ibeg] = -(F0 - F[ibeg]) / (Δz[ibeg] * cv[ibeg])
  @inbounds for i in ibeg+1:n
    dT[i] = -(F[i-1] - F[i]) / (Δz[i] * cv[i])
  end
end

function TsoilEquation_partial(dT, T, p::Soil, t; method="TS0", ibeg::Int=1)
  p._dTsoil[ibeg:end] .= dT
  p._Tsoil[ibeg:end] .= T
  TsoilEquation(p._dTsoil, p._Tsoil, p, t; method, ibeg)
  dT .= p._dTsoil[ibeg:end]
  return nothing
end
