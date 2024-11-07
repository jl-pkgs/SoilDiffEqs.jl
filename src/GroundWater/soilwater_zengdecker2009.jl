# Constants
SHR_CONST_TKFRZ = 273.15
SHR_CONST_LATICE = 3.337e5
SHR_CONST_G = 9.80665

# Local variables
subname = "soilwater_zengdecker2009"

z = col.z
zi = col.zi
dz = col.dz

zmm = zeros(1:N+1)
dzmm = zeros(1:N+1)

amx = zeros(1:N+1)
bmx = zeros(1:N+1)
cmx = zeros(1:N+1)
rmx = zeros(1:N+1)

ψ = zeros(1:N)
ψE = zeros(1:N+1)
θE = zeros(1:N+1)

dqidw0 = zeros(1:N+1)
dqidw1 = zeros(1:N+1)
dqodw1 = zeros(1:N+1)
dqodw2 = zeros(1:N+1)
dψdw = zeros(1:N+1)

# Associate variables
origflag = soilhydrology_inst.origflag
qcharge = soilhydrology_inst.qcharge_col
zwt = soilhydrology_inst.zwt_col
fracice = soilhydrology_inst.fracice_col
icefrac = soilhydrology_inst.icefrac_col
ψE_min = soilstate_inst.ψmin_col
θ_sat = soilstate_inst.θ_sat_col
hksat = soilstate_inst.hksat_col
bsw = soilstate_inst.bsw_col
ψ_sat = soilstate_inst.ψ_sat_col
ψ_l = soilstate_inst.ψ_l_col
hk_l = soilstate_inst.hk_l_col
θ_ice = waterstate_inst.h2osoi_ice_col
θ_liq = waterstate_inst.h2osoi_liq_col
θ = waterstate_inst.h2osoi_vol_col
qflx_deficit = waterflux_inst.qflx_deficit_col
qflx_infl = waterflux_inst.qflx_infl_col
qflx_rootsoi_col = waterflux_inst.qflx_rootsoi_col
t_soisno = temperature_inst.t_soisno_col


function find_jwt(zi::AbstractVector, zwt::Real)
  N = length(zi)
  jwt = N
  for j in 1:N
    if zwt <= zi[j]
      jwt = j - 1
      break
    end
  end
  return jwt
end

function soilwater_zengdecker2009(bounds, num_hydrologyc, filter_hydrologyc,
  soilhydrology_inst, soilstate_inst,
  waterflux_inst, waterstate_inst, temperature_inst)

  hk = zeros(1:N)
  dhkdw = zeros(1:N)

  dz = 0.0
  num = 0.0
  qin = zeros(1:N+1)
  qout = zeros(1:N+1)
  se = 0.0
  s1 = 0.0
  s2 = 0.0
  sdamp = 0.0
  ψ1 = 0.0
  dψdw1 = 0.0
  wh = 0.0
  wh_zwt = 0.0
  _K = 0.0
  dwat2 = zeros(1:N+1)
  dψE = 0.0
  zimm = zeros(0:N)
  tempi = 0.0
  temp0 = 0.0
  voleq1 = 0.0
  vwc_zwt = 0.0
  imped = zeros(1:N)
  vol_ice = zeros(1:N)
  vwc_liq = zeros(1:N+1)

  dtime = get_step_size()

  # Convert depths to mm
  for j in 1:N
    zmm[j] = z[j] * 1e3
    dzmm[j] = dz[j] * 1e3
    zimm[j] = zi[j] * 1e3

    vol_ice[j] = min(θ_sat[j], θ_ice[j] / (dz[j] * ρ_ice))
    icefrac[j] = min(1.0, vol_ice[j] / θ_sat[j])
    vwc_liq[j] = max(θ_liq[j], 1.0e-6) / (dz[j] * ρ_h20)
  end

  zimm[0] = 0.0
  zwtmm = zwt * 1e3

  # Compute jwt index
  jwt = find_jwt(zi, zwt)
  vwc_zwt = θ_sat[N]
  if t_soisno[jwt+1] < tfrz
    vwc_zwt = vwc_liq[N]
    for j in N:nlevgrnd
      if zwt <= zi[j]
        ψ1 = hfus * (tfrz - t_soisno[j]) / (grav * t_soisno[j]) * 1000.0
        ψ1 = max(ψ_sat[N], ψ1)
        vwc_zwt = θ_sat[N] * (ψ1 / ψ_sat[N])^(-1.0 / bsw[N])
        vwc_zwt = min(vwc_zwt, 0.5 * (θ_sat[N] + θ[N]))
        break
      end
    end
  end


  # ! 核心参考部分
  # Calculate the equilibrium water content based on the water table depth
  B = bsw[j]
  c = filter_hydrologyc[fc]
  C = ψ_sat[j] + zwtmm
  _Δz = zimm[j] - zimm[j-1]

  if zwtmm <= zimm[j-1]
    θE[j] = θ_sat[j]
  elseif zimm[j-1] < zwtmm < zimm[j]
    # 部分饱和
    tempi = 1.0
    temp0 = ((C - zimm[j-1]) / ψ_sat[j])^(1 - 1 / B)
    d1 = zwtmm - zimm[j-1] # 未饱和
    d2 = zimm[j] - zwtmm   # 饱和

    _θE = -ψ_sat[j] * θ_sat[j] / (1 - 1 / B) / d1 * (tempi - temp0)
    θE[j] = (_θE * d1 + θ_sat[j] * d2) / _Δz
    θE[j] = clamp(θE[j], 0, θ_sat[j])
  else
    # zwt > z₊ₕ[j], 地下水水位在这一层之下，这一层非饱和
    tempi = ((C - zimm[j]) / ψ_sat[j])^(1 - 1 / B)
    temp0 = ((C - zimm[j-1]) / ψ_sat[j])^(1 - 1 / B)

    θE[j] = -ψ_sat[j] * θ_sat[j] / (1 - 1 / B) / (_Δz) * (tempi - temp0) # Zeng 2009, Eq.9
    θE[j] = clamp(θE[j], 0, θ_sat[j])
  end

  ψE[j] = -ψ_sat[j] * (max(θE[j] / θ_sat[j], 0.01))^(-B) # 
  ψE[j] = max(ψE_min, ψE[j])

  # If water table is below soil column calculate ψE for the 11th layer
  j = N
  if jwt == N
    # 积分的过程在：z₊ₕ[N] ~ zwt，因此，需要注意，最后一层的θ_E、ψ_E代表的区间
    tempi = 1.0
    temp0 = ((ψ_sat[j] + zwtmm - zimm[j]) / ψ_sat[j])^(1.0 - 1.0 / bsw[j])
    θE[j+1] = -ψ_sat[j] * θ_sat[j] / (1.0 - 1.0 / bsw[j]) / (zwtmm - zimm[j]) * (tempi - temp0) # 这里没有进行加权，人为
    θE[j+1] = clamp(θE[j+1], 0.0, θ_sat[j])

    ψE[j+1] = -ψ_sat[j] * (max(θE[j+1] / θ_sat[j], 0.01))^(-bsw[j])
    ψE[j+1] = max(ψE_min, ψE[j+1])
  end

  # Hydraulic conductivity and soil matric potential and their derivatives
  sdamp = 0.0
  for j in 1:N
    if origflag == 1
      s1 = 0.5 * (θ[j] + θ[min(N, j + 1)]) / (0.5 * (θ_sat[j] + θ_sat[min(N, j + 1)]))
    else
      s1 = 0.5 * (vwc_liq[j] + vwc_liq[min(N, j + 1)]) / (0.5 * (θ_sat[j] + θ_sat[min(N, j + 1)]))
    end
    s1 = min(1.0, s1)
    s2 = hksat[j] * s1^(2.0 * bsw[j] + 2.0)
    if origflag == 1
      imped[j] = (1.0 - 0.5 * (fracice[j] + fracice[min(N, j + 1)]))
    else
      imped[j] = 10.0^(-e_ice * (0.5 * (icefrac[j] + icefrac[min(N, j + 1)])))
    end
    hk[j] = imped[j] * s1 * s2
    dhkdw[j] = imped[j] * (2.0 * bsw[j] + 3.0) * s2 * (1.0 / (θ_sat[j] + θ_sat[min(N, j + 1)]))

    _θ = origflag == 1 ? θ[j] : vwc_liq[j]
    se = clamp(_θ / θ_sat[j], 0.01, 1.0)
    ψ[j] = -ψ_sat[j] * se^(-bsw[j])
    ψ[j] = max(ψ[j], ψE_min)
    dψdw[j] = -bsw[j] * ψ[j] / (_θ)

    ψ_l[j] = ψ[j]
    hk_l[j] = hk[j]
  end

  # Aquifer (11th) layer
  zmm[N+1] = 0.5 * (1e3 * zwt + zmm[N]) # 动态调整最后一层的深度，高明!
  if jwt < N
    dzmm[N+1] = dzmm[N] # 这里需要核对
  else
    dzmm[N+1] = (1e3 * zwt - zmm[N])
  end

  # Set up r, a, b, and c vectors for tridiagonal solution
  j = 1
  qin[j] = qflx_infl
  dz = (zmm[j+1] - zmm[j])
  dψE = (ψE[j+1] - ψE[j])
  num = (ψ[j+1] - ψ[j]) - dψE
  qout[j] = -hk[j] * num / dz
  dqodw1[j] = -(-hk[j] * dψdw[j] + num * dhkdw[j]) / dz
  dqodw2[j] = -(hk[j] * dψdw[j+1] + num * dhkdw[j]) / dz
  rmx[j] = qin[j] - qout[j] - qflx_rootsoi_col[j]
  amx[j] = 0.0
  bmx[j] = dzmm[j] * (sdamp + 1.0 / dtime) + dqodw1[j]
  cmx[j] = dqodw2[j]

  # Nodes j=2 to j=N-1
  for j in 2:N-1
    dz = (zmm[j] - zmm[j-1])
    dψE = (ψE[j] - ψE[j-1])
    num = (ψ[j] - ψ[j-1]) - dψE
    qin[j] = -hk[j-1] * num / dz
    dqidw0[j] = -(-hk[j-1] * dψdw[j-1] + num * dhkdw[j-1]) / dz
    dqidw1[j] = -(hk[j-1] * dψdw[j] + num * dhkdw[j-1]) / dz

    dz = (zmm[j+1] - zmm[j])
    dψE = (ψE[j+1] - ψE[j])
    num = (ψ[j+1] - ψ[j]) - dψE
    qout[j] = -hk[j] * num / dz
    dqodw1[j] = -(-hk[j] * dψdw[j] + num * dhkdw[j]) / dz
    dqodw2[j] = -(hk[j] * dψdw[j+1] + num * dhkdw[j]) / dz
    rmx[j] = qin[j] - qout[j] - qflx_rootsoi_col[j]
    amx[j] = -dqidw0[j]
    bmx[j] = dzmm[j] / dtime - dqidw1[j] + dqodw1[j]
    cmx[j] = dqodw2[j]
  end

  # Node j=N (bottom)
  j = N
  if j > jwt  # water table is in soil column
    dz = (zmm[j] - zmm[j-1])
    dψE = (ψE[j] - ψE[j-1])
    num = (ψ[j] - ψ[j-1]) - dψE
    qin[j] = -hk[j-1] * num / dz
    dqidw0[j] = -(-hk[j-1] * dψdw[j-1] + num * dhkdw[j-1]) / dz
    dqidw1[j] = -(hk[j-1] * dψdw[j] + num * dhkdw[j-1]) / dz
    qout[j] = 0.0
    dqodw1[j] = 0.0
    rmx[j] = qin[j] - qout[j] - qflx_rootsoi_col[j]
    amx[j] = -dqidw0[j]
    bmx[j] = dzmm[j] / dtime - dqidw1[j] + dqodw1[j]
    cmx[j] = 0.0

    # next set up aquifer layer; hydrologically inactive
    rmx[j+1] = 0.0
    amx[j+1] = 0.0
    bmx[j+1] = dzmm[j+1] / dtime
    cmx[j+1] = 0.0
  else  # water table is below soil column
    # compute aquifer soil moisture as average of layer 10 and saturation
    if origflag == 1
      se = 0.5 * (1.0 + θ[j] / θ_sat[j])
    else
      se = 0.5 * ((vwc_zwt + vwc_liq[j]) / θ_sat[j])
    end
    se = clamp(se, 0.01, 1.0)

    # compute ψ for aquifer layer
    ψ1 = -ψ_sat[j] * se^(-bsw[j])
    ψ1 = max(ψE_min, ψ1)

    # compute dψdw for aquifer layer
    dψdw1 = -bsw[j] * ψ1 / (se * θ_sat[j])

    # first set up bottom layer of soil column
    dz = (zmm[j] - zmm[j-1])
    dψE = (ψE[j] - ψE[j-1])
    num = (ψ[j] - ψ[j-1]) - dψE
    qin[j] = -hk[j-1] * num / dz
    dqidw0[j] = -(-hk[j-1] * dψdw[j-1] + num * dhkdw[j-1]) / dz
    dqidw1[j] = -(hk[j-1] * dψdw[j] + num * dhkdw[j-1]) / dz

    dz = (zmm[j+1] - zmm[j])
    dψE = (ψE[j+1] - ψE[j])
    num = (ψ1 - ψ[j]) - dψE
    qout[j] = -hk[j] * num / dz
    dqodw1[j] = -(-hk[j] * dψdw[j] + num * dhkdw[j]) / dz
    dqodw2[j] = -(hk[j] * dψdw1 + num * dhkdw[j]) / dz

    rmx[j] = qin[j] - qout[j] - qflx_rootsoi_col[j]
    amx[j] = -dqidw0[j]
    bmx[j] = dzmm[j] / dtime - dqidw1[j] + dqodw1[j]
    cmx[j] = dqodw2[j]

    # next set up aquifer layer; dz/num unchanged, qin=qout
    qin[j+1] = qout[j]
    dqidw0[j+1] = -(-hk[j] * dψdw[j] + num * dhkdw[j]) / dz
    dqidw1[j+1] = -(hk[j] * dψdw1 + num * dhkdw[j]) / dz
    qout[j+1] = 0.0  # zero-flow bottom boundary condition
    dqodw1[j+1] = 0.0  # zero-flow bottom boundary condition
    rmx[j+1] = qin[j+1] - qout[j+1]
    amx[j+1] = -dqidw0[j+1]
    bmx[j+1] = dzmm[j+1] / dtime - dqidw1[j+1] + dqodw1[j+1]
    cmx[j+1] = 0.0
  end

  # Solve for dwat
  jtop = 1
  Tridiagonal(bounds, 1, N + 1, jtop, num_hydrologyc, filter_hydrologyc, amx[:], bmx[:], cmx[:], rmx[:], dwat2[:])

  # Renew the mass of liquid water
  # also compute qcharge from dwat in aquifer layer
  # update in drainage for case jwt < N
  for j in 1:N
    θ_liq[j] += dwat2[j] * dzmm[j]
  end

  j = jwt
  # calculate qcharge for case jwt < N
  if jwt < N
    wh_zwt = 0.0  # since wh_zwt = -ψ_sat - ψE_zwt, where ψE_zwt = -ψ_sat

    # Recharge rate qcharge to groundwater (positive to aquifer)
    se = max(θ[jwt+1] / θ_sat[jwt+1], 0.01)
    s1 = min(1.0, se)

    # scs: this is the expression for unsaturated hk
    _K = imped[jwt+1] * hksat[jwt+1] * s1^(2.0 * bsw[jwt+1] + 3.0)

    ψ1 = max(ψE_min, ψ[max(1, jwt)])
    wh = ψ1 - ψE[max(1, jwt)]  # 这里是向地下水的排泄，Zeng2009, Eq.14

    if jwt == 0
      qcharge = -_K * (wh_zwt - wh) / ((zwt + 1.0e-3) * 1000.0)
    else
      qcharge = -_K * (wh_zwt - wh) / ((zwt - z[jwt]) * 1000.0 * 2.0)
    end
    # To limit qcharge (for the first several timesteps)
    qcharge = max(qcharge, -10.0 / dtime, 10.0 / dtime)
  else
    # if water table is below soil column, compute qcharge from dwat2(11)
    qcharge = dwat2[N+1] * dzmm[N+1] / dtime
  end

  # compute the water deficit and reset negative liquid water content
  qflx_deficit = 0.0
  for j in 1:N
    if θ_liq[j] < 0.0
      qflx_deficit -= θ_liq[j]
    end
  end

end
