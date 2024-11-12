"""
边界条件隐含在`soil`中：ψ0, Q0
"""
function soil_moisture_BEPS(soil::Soil{FT}, θ0::AbstractVector{FT}, θ_surf::AbstractVector{FT}; method="ψ0") where {FT}
  (; ibeg, N) = soil
  soil.θ[ibeg:end] .= θ0
  
  ntime = length(θ_surf)
  nlayer = N - ibeg + 1
  res = zeros(ntime, nlayer)

  @inbounds for i = 1:ntime
    r = soil_moisture_BEPS(soil, θ_surf[i]; method)
    res[i, :] .= r[ibeg:end]
  end
  res
end


function soil_moisture_BEPS(soil::Soil{FT}, θ_surf::FT; method="ψ0") where {FT}
  (; ψ0, Q0, ibeg, ψ) = soil
  kstep = soil.dt # 3600s
  ∑t = 0.0
  θ = soil.θ
  θ[ibeg] = θ_surf

  # 这里是一个时间步长
  while ∑t < kstep
    cal_ψ!(soil, θ)
    ψ0 = ψ[ibeg]

    _Q0 = soil_WaterFlux!(soil, θ; ψ0, Q0, method) # update Q
    Qmax = maximum(abs.(soil.Q))
    Qmax = max(Qmax, abs(_Q0))

    Δt = guess_step(Qmax) # this_step
    ∑t += Δt
    ∑t > kstep && (Δt -= (∑t - kstep))

    soil_Updateθ!(soil, Δt) # update θ during Δt
    ## 误差如何计算？
  end
  θ
end

"""
    guess_step(max_Fb)

    如果流速过快，则减小时间步长

# Arguments
- `max_Fb`: [cm s-1]

# Notes
[m s-1] -> 1000*[mm s-1] -> 1000*[kg m-2 s-1]
"""
function guess_step(max_Fb)
  # this constraint is too large
  # 需要调研，暴雨期间土壤的下渗速率，或者说`Ksat`
  max_Fb = max_Fb * 86400 * 10 # convert to [mm/day]
  # dt = 360.0
  if max_Fb > 100
    30
  elseif max_Fb > 10
    120
  else
    360
  end
  # if max_Fb > 86.4 # 86.4 mm/day
  #   dt = 1.0
  # elseif max_Fb > 8.64 # 8.64 mm/day
  #   dt = 30.0 # seconds
  # else
  #   dt = 360.0
  # end
  # dt
end

export soil_moisture_BEPS
