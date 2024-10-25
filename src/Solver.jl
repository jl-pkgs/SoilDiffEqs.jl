function solve_Tsoil_Bonan(soil::Soil{FT}, TS0::AbstractVector{FT}; ibeg=1) where {FT<:Real}
  ntime = length(TS0)
  (; inds_obs, n) = soil
  R = zeros(ntime, soil.n - ibeg + 1)
  R[1, :] .= soil.Tsoil[ibeg:end]
  
  for i = 2:ntime
    soil_temperature!(soil, TS0[i]; ibeg)
    R[i, :] .= soil.Tsoil[ibeg:end]
  end

  _inds = inds_obs .- ibeg .+ 1
  R[:, _inds]
end
