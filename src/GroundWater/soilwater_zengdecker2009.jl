include("GroundWater.jl")

zimm = zeros(0:N)
zmm = zeros(1:N+1)
dzmm = zeros(1:N+1)

amx = zeros(1:N+1)
bmx = zeros(1:N+1)
cmx = zeros(1:N+1)
rmx = zeros(1:N+1)

ψ = zeros(1:N)
K = zeros(1:N)
ψE = zeros(1:N+1)
θE = zeros(1:N+1)

qin = zeros(1:N+1)
qout = zeros(1:N+1)

dqidw0 = zeros(1:N+1)
dqidw1 = zeros(1:N+1)
dqodw1 = zeros(1:N+1)
dqodw2 = zeros(1:N+1)
dψdw = zeros(1:N+1)
dKdw = zeros(1:N)

# Associate variables
qcharge = soilhydrology_inst.qcharge_col
zwt = soilhydrology_inst.zwt_col
fracice = soilhydrology_inst.fracice_col
icefrac = soilhydrology_inst.icefrac_col
ψmin = soilstate_inst.ψmin_col
θ_sat = soilstate_inst.θ_sat_col
Ksat = soilstate_inst.hksat_col
bsw = soilstate_inst.bsw_col
ψ_sat = soilstate_inst.ψ_sat_col
ψ_l = soilstate_inst.ψ_l_col
K_l = soilstate_inst.hk_l_col
θ_ice = waterstate_inst.h2osoi_ice_col
θ_liq = waterstate_inst.h2osoi_liq_col
θ = waterstate_inst.h2osoi_vol_col
qflx_deficit = waterflux_inst.qflx_deficit_col
qflx_infl = waterflux_inst.qflx_infl_col
qflx_rootsoi_col = waterflux_inst.qflx_rootsoi_col
Tsoil = temperature_inst.t_soisno_col

function soilwater_zengdecker2009()
  sdamp = 0.0
  dψdw1 = 0.0
  dθ = zeros(1:N+1)
  dψE = 0.0
  imped = zeros(1:N)

  vwc_zwt = 0.0
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
  if Tsoil[jwt+1] < tfrz
    vwc_zwt = vwc_liq[N]
    for j in N:nlevgrnd
      if zwt <= zi[j]
        _ψ = hfus * (tfrz - Tsoil[j]) / (grav * Tsoil[j]) * 1000.0
        _ψ = max(ψ_sat[N], _ψ)
        vwc_zwt = θ_sat[N] * (_ψ / ψ_sat[N])^(-1.0 / bsw[N])
        vwc_zwt = min(vwc_zwt, 0.5 * (θ_sat[N] + θ[N]))
        break
      end
    end
  end

  # ! 核心参考部分
  # Calculate the equilibrium water content based on the water table depth
  for j = 1:N
    if zwtmm <= zimm[j-1]
      θE[j] = θ_sat[j]
    elseif zimm[j-1] < zwtmm < zimm[j]
      # 部分饱和
      d1 = zwtmm - zimm[j-1] # 未饱和
      d2 = zimm[j] - zwtmm   # 饱和
      _θE = cal_θE(zwtmm, zimm[j-1], zwtmm, ψ_sat[j], bsw[j]; use_clamp=false)
      θE[j] = (_θE * d1 + θ_sat[j] * d2) / (d1 + d2)
      θE[j] = clamp(θE[j], 0, θ_sat[j])
    else
      # zwt > z₊ₕ[j], 地下水水位在这一层之下，这一层非饱和
      θE[j] = cal_θE(zimm[j], zimm[j-1], zwtmm, ψ_sat[j], bsw[j])
    end
    ψE[j] = cal_ψ(θE[j], θ_sat[j], ψ_sat[j], bsw[j]; ψmin)
  end

  # If water table is below soil column calculate ψE for the 11th layer
  if jwt == N
    # 积分的过程在：z₊ₕ[N] ~ zwt，因此，需要注意，最后一层的θ_E、ψ_E代表的区间
    j = N
    θE[j+1] = cal_θE(zwtmm, zimm[j], ψ_sat[j], zwtmm, bsw[j])
    ψE[j+1] = cal_ψ(θE[j+1], θ_sat[j], ψ_sat[j], bsw[j]; ψmin)
  end

  # Hydraulic conductivity and soil matric potential and their derivatives
  sdamp = 0.0
  for j in 1:N
    if origflag == 1
      se = (θ[j] + θ[min(N, j + 1)]) / ((θ_sat[j] + θ_sat[min(N, j + 1)]))
    else
      se = (vwc_liq[j] + vwc_liq[min(N, j + 1)]) / ((θ_sat[j] + θ_sat[min(N, j + 1)]))
    end
    se = min(1.0, se)

    s2 = Ksat[j] * se^(2.0 * bsw[j] + 2.0)
    if origflag == 1
      imped[j] = (1.0 - 0.5 * (fracice[j] + fracice[min(N, j + 1)]))
    else
      imped[j] = 10.0^(-e_ice * (0.5 * (icefrac[j] + icefrac[min(N, j + 1)])))
    end
    K[j] = imped[j] * se * s2
    dKdw[j] = imped[j] * (2.0 * bsw[j] + 3.0) * s2 * (1.0 / (θ_sat[j] + θ_sat[min(N, j + 1)]))

    _θ = origflag == 1 ? θ[j] : vwc_liq[j]
    ψ[j] = cal_ψ(_θ[j], θ_sat[j], ψ_sat[j], bsw[j]; ψmin)
    dψdw[j] = -bsw[j] * ψ[j] / (_θ)

    ψ_l[j] = ψ[j]
    K_l[j] = K[j]
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
  qout[j] = -K[j] * num / dz
  dqodw1[j] = -(-K[j] * dψdw[j] + num * dKdw[j]) / dz
  dqodw2[j] = -(K[j] * dψdw[j+1] + num * dKdw[j]) / dz
  rmx[j] = qin[j] - qout[j] - qflx_rootsoi_col[j]
  amx[j] = 0.0
  bmx[j] = dzmm[j] * (sdamp + 1.0 / dtime) + dqodw1[j]
  cmx[j] = dqodw2[j]

  # Nodes j=2 to j=N-1
  for j in 2:N-1
    dz = (zmm[j] - zmm[j-1])
    dψE = (ψE[j] - ψE[j-1])
    num = (ψ[j] - ψ[j-1]) - dψE
    qin[j] = -K[j-1] * num / dz
    dqidw0[j] = -(-K[j-1] * dψdw[j-1] + num * dKdw[j-1]) / dz
    dqidw1[j] = -(K[j-1] * dψdw[j] + num * dKdw[j-1]) / dz

    dz = (zmm[j+1] - zmm[j])
    dψE = (ψE[j+1] - ψE[j])
    num = (ψ[j+1] - ψ[j]) - dψE
    qout[j] = -K[j] * num / dz
    dqodw1[j] = -(-K[j] * dψdw[j] + num * dKdw[j]) / dz
    dqodw2[j] = -(K[j] * dψdw[j+1] + num * dKdw[j]) / dz
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
    qin[j] = -K[j-1] * num / dz
    dqidw0[j] = -(-K[j-1] * dψdw[j-1] + num * dKdw[j-1]) / dz
    dqidw1[j] = -(K[j-1] * dψdw[j] + num * dKdw[j-1]) / dz
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

    # compute for aquifer layer
    _ψ = cal_ψ(se, ψ_sat[j], bsw[j]; ψmin)
    dψdw1 = -bsw[j] * _ψ / (se * θ_sat[j])

    # first set up bottom layer of soil column
    dz = (zmm[j] - zmm[j-1])
    dψE = (ψE[j] - ψE[j-1])
    num = (ψ[j] - ψ[j-1]) - dψE
    qin[j] = -K[j-1] * num / dz
    dqidw0[j] = -(-K[j-1] * dψdw[j-1] + num * dKdw[j-1]) / dz
    dqidw1[j] = -(K[j-1] * dψdw[j] + num * dKdw[j-1]) / dz

    dz = (zmm[j+1] - zmm[j])
    dψE = (ψE[j+1] - ψE[j])
    num = (_ψ - ψ[j]) - dψE
    qout[j] = -K[j] * num / dz
    dqodw1[j] = -(-K[j] * dψdw[j] + num * dKdw[j]) / dz
    dqodw2[j] = -(K[j] * dψdw1 + num * dKdw[j]) / dz

    rmx[j] = qin[j] - qout[j] - qflx_rootsoi_col[j]
    amx[j] = -dqidw0[j]
    bmx[j] = dzmm[j] / dtime - dqidw1[j] + dqodw1[j]
    cmx[j] = dqodw2[j]

    # next set up aquifer layer; dz/num unchanged, qin=qout
    qin[j+1] = qout[j]
    dqidw0[j+1] = -(-K[j] * dψdw[j] + num * dKdw[j]) / dz
    dqidw1[j+1] = -(K[j] * dψdw1 + num * dKdw[j]) / dz
    qout[j+1] = 0.0  # zero-flow bottom boundary condition
    dqodw1[j+1] = 0.0  # zero-flow bottom boundary condition
    rmx[j+1] = qin[j+1] - qout[j+1]
    amx[j+1] = -dqidw0[j+1]
    bmx[j+1] = dzmm[j+1] / dtime - dqidw1[j+1] + dqodw1[j+1]
    cmx[j+1] = 0.0
  end

  Tridiagonal(amx, bmx, cmx, rmx, dθ; jtop=1) # Solve for dθ

  # Renew the mass of liquid water also compute qcharge from dθ in aquifer layer
  # update in drainage for case jwt < N
  for j in 1:N
    θ_liq[j] += dθ[j] * dzmm[j]
  end

  # calculate qcharge for case jwt < N
  if jwt < N
    wh_zwt = 0.0  # since wh_zwt = -ψ_sat - ψE_zwt, where ψE_zwt = -ψ_sat

    # Recharge rate qcharge to groundwater (positive to aquifer)
    se = clamp(θ[jwt+1] / θ_sat[jwt+1], 0.01, 1.0)
    # scs: this is the expression for unsaturated K
    _K = imped[jwt+1] * Ksat[jwt+1] * se^(2.0 * bsw[jwt+1] + 3.0)

    _ψ = max(ψmin, ψ[max(1, jwt)])
    wh = _ψ - ψE[max(1, jwt)]  # 这里是向地下水的排泄，Zeng2009, Eq.14

    if jwt == 0
      qcharge = -_K * (wh_zwt - wh) / ((zwt + 1.0e-3) * 1000.0)
    else
      qcharge = -_K * (wh_zwt - wh) / ((zwt - z[jwt]) * 1000.0 * 2.0)
    end
    # To limit qcharge (for the first several timesteps)
    qcharge = max(qcharge, -10.0 / dtime, 10.0 / dtime)
  else
    # if water table is below soil column, compute qcharge from dθ(11)
    qcharge = dθ[N+1] * dzmm[N+1] / dtime
  end

  # compute the water deficit and reset negative liquid water content
  qflx_deficit = 0.0
  for j in 1:N
    if θ_liq[j] < 0.0
      qflx_deficit -= θ_liq[j]
    end
  end

end
