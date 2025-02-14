## 速度太慢，不具有可行性
"""
边界条件隐含在`soil`中：ψ0, Q0
"""
function soil_moisture_BEPS(soil::Soil{FT}, θ0::AbstractVector{FT}, θ_surf::AbstractVector{FT}; method="ψ0") where {FT}
  (; ibeg, N) = soil
  i0 = max(ibeg - 1, 1)
  soil.θ[i0:end] .= θ0
  
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
  (; ψ0, Q, ibeg, ψ) = soil
  kstep = soil.dt # 3600s
  ∑t = 0.0
  θ = soil.θ
  i0 = max(ibeg - 1, 1)
  θ[i0] = θ_surf

  # 这里是一个时间步长
  n = 0
  while ∑t < kstep
    n += 1
    cal_ψ!(soil, θ)
    ψ0 = ψ[i0]

    Q0 = cal_Q!(soil, θ; ψ0, method) # update Q
    
    Qmax = guess_Qmax(Q, Q0; ibeg)
    Δt = guess_step(Qmax) # this_step
    ∑t += Δt
    ∑t > kstep && (Δt -= (∑t - kstep))

    # mod(n, 20) == 0 && println("Qmax = $Qmax, Δt = $Δt")
    # println("Qmax = $Qmax, Δt = $Δt")
    soil_Updateθ!(soil, Δt) # update θ during Δt
  end
  θ
end


function guess_Qmax(Q::AbstractVector{FT}, Q0::FT=0; ibeg::Int=1) where {FT}
  Qmax = 0.0
  @inbounds for i = ibeg:length(Q)
    val = abs(Q[i])
    val > Qmax && (Qmax = val)
  end
  return max(Qmax, abs(Q0))
end

"""
    guess_step(max_Fb)

    如果流速过快，则减小时间步长

# Arguments
- `max_Fb`: [cm h-1]
"""
function guess_step(Qmax::T) where {T}
  Qmax = Qmax * 24 * 10 # [cm h-1] -> [mm d-1]
  # this constraint is too large
  # 需要调研，暴雨期间土壤的下渗速率，或者说`Ksat`
  if Qmax >= 100 # [mm d-1]
    T(30)
  elseif Qmax >= 20
    T(120)
  else
    T(360)
  end
end
# if max_Fb > 86.4 # 86.4 mm/day
#   dt = 1.0
# elseif max_Fb > 8.64 # 8.64 mm/day
#   dt = 30.0 # seconds
# else
#   dt = 360.0
# end

export soil_moisture_BEPS
export guess_Qmax
