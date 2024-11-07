# Constants
SHR_CONST_TKFRZ = 273.15
SHR_CONST_LATICE = 3.337e5
SHR_CONST_G = 9.80665

# Local variables
subname = "soilwater_zengdecker2009"

function soilwater_zengdecker2009(bounds, num_hydrologyc, filter_hydrologyc,
  num_urbanc, filter_urbanc, soilhydrology_inst, soilstate_inst,
  waterflux_inst, waterstate_inst, temperature_inst,
  canopystate_inst, energyflux_inst, soil_water_retention_curve)

  jtop = zeros(Int, bounds.begc:bounds.endc)
  dtime = 0.0
  hk = zeros(bounds.begc:bounds.endc, 1:N)
  dhkdw = zeros(bounds.begc:bounds.endc, 1:N)
  amx = zeros(bounds.begc:bounds.endc, 1:N+1)
  bmx = zeros(bounds.begc:bounds.endc, 1:N+1)
  cmx = zeros(bounds.begc:bounds.endc, 1:N+1)
  rmx = zeros(bounds.begc:bounds.endc, 1:N+1)
  zmm = zeros(bounds.begc:bounds.endc, 1:N+1)
  dzmm = zeros(bounds.begc:bounds.endc, 1:N+1)
  den = 0.0
  dqidw0 = zeros(bounds.begc:bounds.endc, 1:N+1)
  dqidw1 = zeros(bounds.begc:bounds.endc, 1:N+1)
  dqodw1 = zeros(bounds.begc:bounds.endc, 1:N+1)
  dqodw2 = zeros(bounds.begc:bounds.endc, 1:N+1)
  dsmpdw = zeros(bounds.begc:bounds.endc, 1:N+1)
  num = 0.0
  qin = zeros(bounds.begc:bounds.endc, 1:N+1)
  qout = zeros(bounds.begc:bounds.endc, 1:N+1)
  s_node = 0.0
  s1 = 0.0
  s2 = 0.0
  smp = zeros(bounds.begc:bounds.endc, 1:N)
  sdamp = 0.0
  pi = 0
  temp = zeros(bounds.begc:bounds.endc)
  jwt = zeros(Int, bounds.begc:bounds.endc)
  smp1 = 0.0
  dsmpdw1 = 0.0
  wh = 0.0
  wh_zwt = 0.0
  ka = 0.0
  dwat2 = zeros(bounds.begc:bounds.endc, 1:N+1)
  dzq = 0.0
  zimm = zeros(bounds.begc:bounds.endc, 0:N)
  zq = zeros(bounds.begc:bounds.endc, 1:N+1)
  θ_E = zeros(bounds.begc:bounds.endc, 1:N+1)
  tempi = 0.0
  temp0 = 0.0
  voleq1 = 0.0
  zwtmm = zeros(bounds.begc:bounds.endc)
  imped = zeros(bounds.begc:bounds.endc, 1:N)
  vol_ice = zeros(bounds.begc:bounds.endc, 1:N)
  z_mid = 0.0
  vwc_zwt = zeros(bounds.begc:bounds.endc)
  vwc_liq = zeros(bounds.begc:bounds.endc, 1:N+1)
  smp_grad = zeros(bounds.begc:bounds.endc, 1:N+1)
  dsmpds = 0.0
  dhkds = 0.0
  hktmp = 0.0
  nstep = 0

  # Associate variables
  z = col.z
  zi = col.zi
  dz = col.dz
  origflag = soilhydrology_inst.origflag
  qcharge = soilhydrology_inst.qcharge_col
  zwt = soilhydrology_inst.zwt_col
  fracice = soilhydrology_inst.fracice_col
  icefrac = soilhydrology_inst.icefrac_col
  hkdepth = soilhydrology_inst.hkdepth_col
  smpmin = soilstate_inst.smpmin_col
  θ_sat = soilstate_inst.θ_sat_col
  hksat = soilstate_inst.hksat_col
  bsw = soilstate_inst.bsw_col
  ψ_sat = soilstate_inst.ψ_sat_col
  eff_porosity = soilstate_inst.eff_porosity_col
  smp_l = soilstate_inst.smp_l_col
  hk_l = soilstate_inst.hk_l_col
  h2osoi_ice = waterstate_inst.h2osoi_ice_col
  h2osoi_liq = waterstate_inst.h2osoi_liq_col
  h2osoi_vol = waterstate_inst.h2osoi_vol_col
  qflx_deficit = waterflux_inst.qflx_deficit_col
  qflx_infl = waterflux_inst.qflx_infl_col
  qflx_rootsoi_col = waterflux_inst.qflx_rootsoi_col
  qflx_tran_veg_col = waterflux_inst.qflx_tran_veg_col
  rootr_col = soilstate_inst.rootr_col
  t_soisno = temperature_inst.t_soisno_col

  # Get time step
  nstep = get_nstep()
  dtime = get_step_size()

  # Convert depths to mm
  for j in 1:N
    for fc in 1:num_hydrologyc
      c = filter_hydrologyc[fc]
      zmm[c, j] = z[c, j] * 1e3
      dzmm[c, j] = dz[c, j] * 1e3
      zimm[c, j] = zi[c, j] * 1e3
      vol_ice[c, j] = min(θ_sat[c, j], h2osoi_ice[c, j] / (dz[c, j] * denice))
      icefrac[c, j] = min(1.0, vol_ice[c, j] / θ_sat[c, j])
      vwc_liq[c, j] = max(h2osoi_liq[c, j], 1.0e-6) / (dz[c, j] * denh2o)
    end
  end

  for fc in 1:num_hydrologyc
    c = filter_hydrologyc[fc]
    zimm[c, 0] = 0.0
    zwtmm[c] = zwt[c] * 1e3
  end
  
  # Compute jwt index
  for fc in 1:num_hydrologyc
    c = filter_hydrologyc[fc]
    jwt[c] = N
    for j in 1:N
      if zwt[c] <= zi[c, j]
        jwt[c] = j - 1
        break
      end
    end
    vwc_zwt[c] = θ_sat[c, N]
    if t_soisno[c, jwt[c]+1] < tfrz
      vwc_zwt[c] = vwc_liq[c, N]
      for j in N:nlevgrnd
        if zwt[c] <= zi[c, j]
          smp1 = hfus * (tfrz - t_soisno[c, j]) / (grav * t_soisno[c, j]) * 1000.0
          smp1 = max(ψ_sat[c, N], smp1)
          vwc_zwt[c] = θ_sat[c, N] * (smp1 / ψ_sat[c, N])^(-1.0 / bsw[c, N])
          vwc_zwt[c] = min(vwc_zwt[c], 0.5 * (θ_sat[c, N] + h2osoi_vol[c, N]))
          break
        end
      end
    end
  end

  # Calculate the equilibrium water content based on the water table depth
  for j in 1:N
    for fc in 1:num_hydrologyc
      c = filter_hydrologyc[fc]
      if zwtmm[c] <= zimm[c, j-1]
        θ_E[c, j] = θ_sat[c, j]
      elseif zwtmm[c] < zimm[c, j] && zwtmm[c] > zimm[c, j-1]
        tempi = 1.0
        temp0 = ((ψ_sat[c, j] + zwtmm[c] - zimm[c, j-1]) / ψ_sat[c, j])^(1.0 - 1.0 / bsw[c, j])
        voleq1 = -ψ_sat[c, j] * θ_sat[c, j] / (1.0 - 1.0 / bsw[c, j]) / (zwtmm[c] - zimm[c, j-1]) * (tempi - temp0)
        θ_E[c, j] = (voleq1 * (zwtmm[c] - zimm[c, j-1]) + θ_sat[c, j] * (zimm[c, j] - zwtmm[c])) / (zimm[c, j] - zimm[c, j-1])
        θ_E[c, j] = min(θ_sat[c, j], θ_E[c, j])
        θ_E[c, j] = max(θ_E[c, j], 0.0)
      else
        tempi = ((ψ_sat[c, j] + zwtmm[c] - zimm[c, j]) / ψ_sat[c, j])^(1.0 - 1.0 / bsw[c, j])
        temp0 = ((ψ_sat[c, j] + zwtmm[c] - zimm[c, j-1]) / ψ_sat[c, j])^(1.0 - 1.0 / bsw[c, j])
        θ_E[c, j] = -ψ_sat[c, j] * θ_sat[c, j] / (1.0 - 1.0 / bsw[c, j]) / (zimm[c, j] - zimm[c, j-1]) * (tempi - temp0)
        θ_E[c, j] = max(θ_E[c, j], 0.0)
        θ_E[c, j] = min(θ_sat[c, j], θ_E[c, j])
      end
      zq[c, j] = -ψ_sat[c, j] * (max(θ_E[c, j] / θ_sat[c, j], 0.01))^(-bsw[c, j])
      zq[c, j] = max(smpmin[c], zq[c, j])
    end
  end

  # If water table is below soil column calculate zq for the 11th layer
  j = N
  for fc in 1:num_hydrologyc
    c = filter_hydrologyc[fc]
    if jwt[c] == N
      tempi = 1.0
      temp0 = ((ψ_sat[c, j] + zwtmm[c] - zimm[c, j]) / ψ_sat[c, j])^(1.0 - 1.0 / bsw[c, j])
      θ_E[c, j+1] = -ψ_sat[c, j] * θ_sat[c, j] / (1.0 - 1.0 / bsw[c, j]) / (zwtmm[c] - zimm[c, j]) * (tempi - temp0)
      θ_E[c, j+1] = max(θ_E[c, j+1], 0.0)
      θ_E[c, j+1] = min(θ_sat[c, j], θ_E[c, j+1])
      zq[c, j+1] = -ψ_sat[c, j] * (max(θ_E[c, j+1] / θ_sat[c, j], 0.01))^(-bsw[c, j])
      zq[c, j+1] = max(smpmin[c], zq[c, j+1])
    end
  end

  # Hydraulic conductivity and soil matric potential and their derivatives
  sdamp = 0.0
  for j in 1:N
    for fc in 1:num_hydrologyc
      c = filter_hydrologyc[fc]
      if origflag == 1
        s1 = 0.5 * (h2osoi_vol[c, j] + h2osoi_vol[c, min(N, j + 1)]) / (0.5 * (θ_sat[c, j] + θ_sat[c, min(N, j + 1)]))
      else
        s1 = 0.5 * (vwc_liq[c, j] + vwc_liq[c, min(N, j + 1)]) / (0.5 * (θ_sat[c, j] + θ_sat[c, min(N, j + 1)]))
      end
      s1 = min(1.0, s1)
      s2 = hksat[c, j] * s1^(2.0 * bsw[c, j] + 2.0)
      if origflag == 1
        imped[c, j] = (1.0 - 0.5 * (fracice[c, j] + fracice[c, min(N, j + 1)]))
      else
        imped[c, j] = 10.0^(-e_ice * (0.5 * (icefrac[c, j] + icefrac[c, min(N, j + 1)])))
      end
      hk[c, j] = imped[c, j] * s1 * s2
      dhkdw[c, j] = imped[c, j] * (2.0 * bsw[c, j] + 3.0) * s2 * (1.0 / (θ_sat[c, j] + θ_sat[c, min(N, j + 1)]))
      if origflag == 1
        s_node = max(h2osoi_vol[c, j] / θ_sat[c, j], 0.01)
      else
        s_node = max(vwc_liq[c, j] / θ_sat[c, j], 0.01)
      end
      s_node = min(1.0, s_node)
      smp[c, j] = -ψ_sat[c, j] * s_node^(-bsw[c, j])
      smp[c, j] = max(smpmin[c], smp[c, j])
      if origflag == 1
        dsmpdw[c, j] = -bsw[c, j] * smp[c, j] / (s_node * θ_sat[c, j])
      else
        dsmpdw[c, j] = -bsw[c, j] * smp[c, j] / vwc_liq[c, j]
      end
      smp_l[c, j] = smp[c, j]
      hk_l[c, j] = hk[c, j]
    end
  end

  # Aquifer (11th) layer
  for fc in 1:num_hydrologyc
    c = filter_hydrologyc[fc]
    zmm[c, N+1] = 0.5 * (1e3 * zwt[c] + zmm[c, N])
    if jwt[c] < N
      dzmm[c, N+1] = dzmm[c, N]
    else
      dzmm[c, N+1] = (1e3 * zwt[c] - zmm[c, N])
    end
  end

  # Set up r, a, b, and c vectors for tridiagonal solution
  j = 1
  for fc in 1:num_hydrologyc
    c = filter_hydrologyc[fc]
    qin[c, j] = qflx_infl[c]
    den = (zmm[c, j+1] - zmm[c, j])
    dzq = (zq[c, j+1] - zq[c, j])
    num = (smp[c, j+1] - smp[c, j]) - dzq
    qout[c, j] = -hk[c, j] * num / den
    dqodw1[c, j] = -(-hk[c, j] * dsmpdw[c, j] + num * dhkdw[c, j]) / den
    dqodw2[c, j] = -(hk[c, j] * dsmpdw[c, j+1] + num * dhkdw[c, j]) / den
    rmx[c, j] = qin[c, j] - qout[c, j] - qflx_rootsoi_col[c, j]
    amx[c, j] = 0.0
    bmx[c, j] = dzmm[c, j] * (sdamp + 1.0 / dtime) + dqodw1[c, j]
    cmx[c, j] = dqodw2[c, j]
  end

  # Associate variables
  z = col.z
  zi = col.zi
  dz = col.dz
  origflag = soilhydrology_inst.origflag
  qcharge = soilhydrology_inst.qcharge_col
  zwt = soilhydrology_inst.zwt_col
  fracice = soilhydrology_inst.fracice_col
  icefrac = soilhydrology_inst.icefrac_col
  hkdepth = soilhydrology_inst.hkdepth_col
  smpmin = soilstate_inst.smpmin_col
  θ_sat = soilstate_inst.θ_sat_col
  hksat = soilstate_inst.hksat_col
  bsw = soilstate_inst.bsw_col
  ψ_sat = soilstate_inst.ψ_sat_col
  eff_porosity = soilstate_inst.eff_porosity_col
  smp_l = soilstate_inst.smp_l_col
  hk_l = soilstate_inst.hk_l_col
  h2osoi_ice = waterstate_inst.h2osoi_ice_col
  h2osoi_liq = waterstate_inst.h2osoi_liq_col
  h2osoi_vol = waterstate_inst.h2osoi_vol_col
  qflx_deficit = waterflux_inst.qflx_deficit_col
  qflx_infl = waterflux_inst.qflx_infl_col
  qflx_rootsoi_col = waterflux_inst.qflx_rootsoi_col
  qflx_tran_veg_col = waterflux_inst.qflx_tran_veg_col
  rootr_col = soilstate_inst.rootr_col
  t_soisno = temperature_inst.t_soisno_col

  # Get time step
  nstep = get_nstep()
  dtime = get_step_size()

  # Convert depths to mm
  for j in 1:N
    for fc in 1:num_hydrologyc
      c = filter_hydrologyc[fc]
      zmm[c, j] = z[c, j] * 1e3
      dzmm[c, j] = dz[c, j] * 1e3
      zimm[c, j] = zi[c, j] * 1e3
      vol_ice[c, j] = min(θ_sat[c, j], h2osoi_ice[c, j] / (dz[c, j] * denice))
      icefrac[c, j] = min(1.0, vol_ice[c, j] / θ_sat[c, j])
      vwc_liq[c, j] = max(h2osoi_liq[c, j], 1.0e-6) / (dz[c, j] * denh2o)
    end
  end

  for fc in 1:num_hydrologyc
    c = filter_hydrologyc[fc]
    zimm[c, 0] = 0.0
    zwtmm[c] = zwt[c] * 1e3
  end

  # Compute jwt index
  for fc in 1:num_hydrologyc
    c = filter_hydrologyc[fc]
    jwt[c] = N
    for j in 1:N
      if zwt[c] <= zi[c, j]
        jwt[c] = j - 1
        break
      end
    end
    vwc_zwt[c] = θ_sat[c, N]
    if t_soisno[c, jwt[c]+1] < tfrz
      vwc_zwt[c] = vwc_liq[c, N]
      for j in N:nlevgrnd
        if zwt[c] <= zi[c, j]
          smp1 = hfus * (tfrz - t_soisno[c, j]) / (grav * t_soisno[c, j]) * 1000.0
          smp1 = max(ψ_sat[c, N], smp1)
          vwc_zwt[c] = θ_sat[c, N] * (smp1 / ψ_sat[c, N])^(-1.0 / bsw[c, N])
          vwc_zwt[c] = min(vwc_zwt[c], 0.5 * (θ_sat[c, N] + h2osoi_vol[c, N]))
          break
        end
      end
    end
  end

  # Calculate the equilibrium water content based on the water table depth
  for j in 1:N
    for fc in 1:num_hydrologyc
      c = filter_hydrologyc[fc]
      if zwtmm[c] <= zimm[c, j-1]
        θ_E[c, j] = θ_sat[c, j]
      elseif zwtmm[c] < zimm[c, j] && zwtmm[c] > zimm[c, j-1]
        tempi = 1.0
        temp0 = ((ψ_sat[c, j] + zwtmm[c] - zimm[c, j-1]) / ψ_sat[c, j])^(1.0 - 1.0 / bsw[c, j])
        voleq1 = -ψ_sat[c, j] * θ_sat[c, j] / (1.0 - 1.0 / bsw[c, j]) / (zwtmm[c] - zimm[c, j-1]) * (tempi - temp0)
        θ_E[c, j] = (voleq1 * (zwtmm[c] - zimm[c, j-1]) + θ_sat[c, j] * (zimm[c, j] - zwtmm[c])) / (zimm[c, j] - zimm[c, j-1])
        θ_E[c, j] = min(θ_sat[c, j], θ_E[c, j])
        θ_E[c, j] = max(θ_E[c, j], 0.0)
      else
        tempi = ((ψ_sat[c, j] + zwtmm[c] - zimm[c, j]) / ψ_sat[c, j])^(1.0 - 1.0 / bsw[c, j])
        temp0 = ((ψ_sat[c, j] + zwtmm[c] - zimm[c, j-1]) / ψ_sat[c, j])^(1.0 - 1.0 / bsw[c, j])
        θ_E[c, j] = -ψ_sat[c, j] * θ_sat[c, j] / (1.0 - 1.0 / bsw[c, j]) / (zimm[c, j] - zimm[c, j-1]) * (tempi - temp0)
        θ_E[c, j] = max(θ_E[c, j], 0.0)
        θ_E[c, j] = min(θ_sat[c, j], θ_E[c, j])
      end
      zq[c, j] = -ψ_sat[c, j] * (max(θ_E[c, j] / θ_sat[c, j], 0.01))^(-bsw[c, j])
      zq[c, j] = max(smpmin[c], zq[c, j])
    end
  end

  # If water table is below soil column calculate zq for the 11th layer
  j = N
  for fc in 1:num_hydrologyc
    c = filter_hydrologyc[fc]
    if jwt[c] == N
      tempi = 1.0
      temp0 = ((ψ_sat[c, j] + zwtmm[c] - zimm[c, j]) / ψ_sat[c, j])^(1.0 - 1.0 / bsw[c, j])
      θ_E[c, j+1] = -ψ_sat[c, j] * θ_sat[c, j] / (1.0 - 1.0 / bsw[c, j]) / (zwtmm[c] - zimm[c, j]) * (tempi - temp0)
      θ_E[c, j+1] = max(θ_E[c, j+1], 0.0)
      θ_E[c, j+1] = min(θ_sat[c, j], θ_E[c, j+1])
      zq[c, j+1] = -ψ_sat[c, j] * (max(θ_E[c, j+1] / θ_sat[c, j], 0.01))^(-bsw[c, j])
      zq[c, j+1] = max(smpmin[c], zq[c, j+1])
    end
  end

  # Hydraulic conductivity and soil matric potential and their derivatives
  sdamp = 0.0
  for j in 1:N
    for fc in 1:num_hydrologyc
      c = filter_hydrologyc[fc]
      if origflag == 1
        s1 = 0.5 * (h2osoi_vol[c, j] + h2osoi_vol[c, min(N, j + 1)]) / (0.5 * (θ_sat[c, j] + θ_sat[c, min(N, j + 1)]))
      else
        s1 = 0.5 * (vwc_liq[c, j] + vwc_liq[c, min(N, j + 1)]) / (0.5 * (θ_sat[c, j] + θ_sat[c, min(N, j + 1)]))
      end
      s1 = min(1.0, s1)
      s2 = hksat[c, j] * s1^(2.0 * bsw[c, j] + 2.0)
      if origflag == 1
        imped[c, j] = (1.0 - 0.5 * (fracice[c, j] + fracice[c, min(N, j + 1)]))
      else
        imped[c, j] = 10.0^(-e_ice * (0.5 * (icefrac[c, j] + icefrac[c, min(N, j + 1)])))
      end
      hk[c, j] = imped[c, j] * s1 * s2
      dhkdw[c, j] = imped[c, j] * (2.0 * bsw[c, j] + 3.0) * s2 * (1.0 / (θ_sat[c, j] + θ_sat[c, min(N, j + 1)]))
      if origflag == 1
        s_node = max(h2osoi_vol[c, j] / θ_sat[c, j], 0.01)
      else
        s_node = max(vwc_liq[c, j] / θ_sat[c, j], 0.01)
      end
      s_node = min(1.0, s_node)
      smp[c, j] = -ψ_sat[c, j] * s_node^(-bsw[c, j])
      smp[c, j] = max(smpmin[c], smp[c, j])
      if origflag == 1
        dsmpdw[c, j] = -bsw[c, j] * smp[c, j] / (s_node * θ_sat[c, j])
      else
        dsmpdw[c, j] = -bsw[c, j] * smp[c, j] / vwc_liq[c, j]
      end
      smp_l[c, j] = smp[c, j]
      hk_l[c, j] = hk[c, j]
    end
  end

  # Aquifer (11th) layer
  for fc in 1:num_hydrologyc
    c = filter_hydrologyc[fc]
    zmm[c, N+1] = 0.5 * (1e3 * zwt[c] + zmm[c, N])
    if jwt[c] < N
      dzmm[c, N+1] = dzmm[c, N]
    else
      dzmm[c, N+1] = (1e3 * zwt[c] - zmm[c, N])
    end
  end

  # Set up r, a, b, and c vectors for tridiagonal solution
  j = 1
  for fc in 1:num_hydrologyc
    c = filter_hydrologyc[fc]
    qin[c, j] = qflx_infl[c]
    den = (zmm[c, j+1] - zmm[c, j])
    dzq = (zq[c, j+1] - zq[c, j])
    num = (smp[c, j+1] - smp[c, j]) - dzq
    qout[c, j] = -hk[c, j] * num / den
    dqodw1[c, j] = -(-hk[c, j] * dsmpdw[c, j] + num * dhkdw[c, j]) / den
    dqodw2[c, j] = -(hk[c, j] * dsmpdw[c, j+1] + num * dhkdw[c, j]) / den
    rmx[c, j] = qin[c, j] - qout[c, j] - qflx_rootsoi_col[c, j]
    amx[c, j] = 0.0
    bmx[c, j] = dzmm[c, j] * (sdamp + 1.0 / dtime) + dqodw1[c, j]
    cmx[c, j] = dqodw2[c, j]
  end

  # Nodes j=2 to j=N-1
  for j in 2:N-1
    for fc in 1:num_hydrologyc
      c = filter_hydrologyc[fc]
      den = (zmm[c, j] - zmm[c, j-1])
      dzq = (zq[c, j] - zq[c, j-1])
      num = (smp[c, j] - smp[c, j-1]) - dzq
      qin[c, j] = -hk[c, j-1] * num / den
      dqidw0[c, j] = -(-hk[c, j-1] * dsmpdw[c, j-1] + num * dhkdw[c, j-1]) / den
      dqidw1[c, j] = -(hk[c, j-1] * dsmpdw[c, j] + num * dhkdw[c, j-1]) / den
      den = (zmm[c, j+1] - zmm[c, j])
      dzq = (zq[c, j+1] - zq[c, j])
      num = (smp[c, j+1] - smp[c, j]) - dzq
      qout[c, j] = -hk[c, j] * num / den
      dqodw1[c, j] = -(-hk[c, j] * dsmpdw[c, j] + num * dhkdw[c, j]) / den
      dqodw2[c, j] = -(hk[c, j] * dsmpdw[c, j+1] + num * dhkdw[c, j]) / den
      rmx[c, j] = qin[c, j] - qout[c, j] - qflx_rootsoi_col[c, j]
      amx[c, j] = -dqidw0[c, j]
      bmx[c, j] = dzmm[c, j] / dtime - dqidw1[c, j] + dqodw1[c, j]
      cmx[c, j] = dqodw2[c, j]
    end
  end

  # Node j=N (bottom)

  j = N
  for fc in 1:num_hydrologyc
    c = filter_hydrologyc[fc]
    if j > jwt[c]  # water table is in soil column
      den = (zmm[c, j] - zmm[c, j-1])
      dzq = (zq[c, j] - zq[c, j-1])
      num = (smp[c, j] - smp[c, j-1]) - dzq
      qin[c, j] = -hk[c, j-1] * num / den
      dqidw0[c, j] = -(-hk[c, j-1] * dsmpdw[c, j-1] + num * dhkdw[c, j-1]) / den
      dqidw1[c, j] = -(hk[c, j-1] * dsmpdw[c, j] + num * dhkdw[c, j-1]) / den
      qout[c, j] = 0.0
      dqodw1[c, j] = 0.0
      rmx[c, j] = qin[c, j] - qout[c, j] - qflx_rootsoi_col[c, j]
      amx[c, j] = -dqidw0[c, j]
      bmx[c, j] = dzmm[c, j] / dtime - dqidw1[c, j] + dqodw1[c, j]
      cmx[c, j] = 0.0

      # next set up aquifer layer; hydrologically inactive
      rmx[c, j+1] = 0.0
      amx[c, j+1] = 0.0
      bmx[c, j+1] = dzmm[c, j+1] / dtime
      cmx[c, j+1] = 0.0
    else  # water table is below soil column

      # compute aquifer soil moisture as average of layer 10 and saturation
      if origflag == 1
        s_node = max(0.5 * (1.0 + h2osoi_vol[c, j] / θ_sat[c, j]), 0.01)
      else
        s_node = max(0.5 * ((vwc_zwt[c] + vwc_liq[c, j]) / θ_sat[c, j]), 0.01)
      end
      s_node = min(1.0, s_node)

      # compute smp for aquifer layer
      smp1 = -ψ_sat[c, j] * s_node^(-bsw[c, j])
      smp1 = max(smpmin[c], smp1)

      # compute dsmpdw for aquifer layer
      dsmpdw1 = -bsw[c, j] * smp1 / (s_node * θ_sat[c, j])

      # first set up bottom layer of soil column
      den = (zmm[c, j] - zmm[c, j-1])
      dzq = (zq[c, j] - zq[c, j-1])
      num = (smp[c, j] - smp[c, j-1]) - dzq
      qin[c, j] = -hk[c, j-1] * num / den
      dqidw0[c, j] = -(-hk[c, j-1] * dsmpdw[c, j-1] + num * dhkdw[c, j-1]) / den
      dqidw1[c, j] = -(hk[c, j-1] * dsmpdw[c, j] + num * dhkdw[c, j-1]) / den
      den = (zmm[c, j+1] - zmm[c, j])
      dzq = (zq[c, j+1] - zq[c, j])
      num = (smp1 - smp[c, j]) - dzq
      qout[c, j] = -hk[c, j] * num / den
      dqodw1[c, j] = -(-hk[c, j] * dsmpdw[c, j] + num * dhkdw[c, j]) / den
      dqodw2[c, j] = -(hk[c, j] * dsmpdw1 + num * dhkdw[c, j]) / den

      rmx[c, j] = qin[c, j] - qout[c, j] - qflx_rootsoi_col[c, j]
      amx[c, j] = -dqidw0[c, j]
      bmx[c, j] = dzmm[c, j] / dtime - dqidw1[c, j] + dqodw1[c, j]
      cmx[c, j] = dqodw2[c, j]

      # next set up aquifer layer; den/num unchanged, qin=qout
      qin[c, j+1] = qout[c, j]
      dqidw0[c, j+1] = -(-hk[c, j] * dsmpdw[c, j] + num * dhkdw[c, j]) / den
      dqidw1[c, j+1] = -(hk[c, j] * dsmpdw1 + num * dhkdw[c, j]) / den
      qout[c, j+1] = 0.0  # zero-flow bottom boundary condition
      dqodw1[c, j+1] = 0.0  # zero-flow bottom boundary condition
      rmx[c, j+1] = qin[c, j+1] - qout[c, j+1]
      amx[c, j+1] = -dqidw0[c, j+1]
      bmx[c, j+1] = dzmm[c, j+1] / dtime - dqidw1[c, j+1] + dqodw1[c, j+1]
      cmx[c, j+1] = 0.0
    end
  end

  # Solve for dwat

  jtop[bounds.begc:bounds.endc] .= 1
  Tridiagonal(bounds, 1, N + 1, jtop[bounds.begc:bounds.endc], num_hydrologyc, filter_hydrologyc, amx[bounds.begc:bounds.endc, :], bmx[bounds.begc:bounds.endc, :], cmx[bounds.begc:bounds.endc, :], rmx[bounds.begc:bounds.endc, :], dwat2[bounds.begc:bounds.endc, :])

  # Renew the mass of liquid water
  # also compute qcharge from dwat in aquifer layer
  # update in drainage for case jwt < N

  for fc in 1:num_hydrologyc
    c = filter_hydrologyc[fc]

    for j in 1:N
      h2osoi_liq[c, j] += dwat2[c, j] * dzmm[c, j]
    end

    # calculate qcharge for case jwt < N
    if jwt[c] < N
      wh_zwt = 0.0  # since wh_zwt = -ψ_sat - zq_zwt, where zq_zwt = -ψ_sat

      # Recharge rate qcharge to groundwater (positive to aquifer)
      s_node = max(h2osoi_vol[c, jwt[c]+1] / θ_sat[c, jwt[c]+1], 0.01)
      s1 = min(1.0, s_node)

      # scs: this is the expression for unsaturated hk
      ka = imped[c, jwt[c]+1] * hksat[c, jwt[c]+1] * s1^(2.0 * bsw[c, jwt[c]+1] + 3.0)

      smp1 = max(smpmin[c], smp[c, max(1, jwt[c])])
      wh = smp1 - zq[c, max(1, jwt[c])]

      if jwt[c] == 0
        qcharge[c] = -ka * (wh_zwt - wh) / ((zwt[c] + 1.0e-3) * 1000.0)
      else
        qcharge[c] = -ka * (wh_zwt - wh) / ((zwt[c] - z[c, jwt[c]]) * 1000.0 * 2.0)
      end

      # To limit qcharge (for the first several timesteps)
      qcharge[c] = max(-10.0 / dtime, qcharge[c])
      qcharge[c] = min(10.0 / dtime, qcharge[c])
    else
      # if water table is below soil column, compute qcharge from dwat2(11)
      qcharge[c] = dwat2[c, N+1] * dzmm[c, N+1] / dtime
    end
  end

  # compute the water deficit and reset negative liquid water content
  for fc in 1:num_hydrologyc
    c = filter_hydrologyc[fc]
    qflx_deficit[c] = 0.0
    for j in 1:N
      if h2osoi_liq[c, j] < 0.0
        qflx_deficit[c] -= h2osoi_liq[c, j]
      end
    end
  end
end
