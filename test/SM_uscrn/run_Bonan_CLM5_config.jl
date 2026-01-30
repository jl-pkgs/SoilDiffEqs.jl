using YAML, SoilDifferentialEquations, Test
import RTableTools: fread
import ModelParams: GOF, sceua
using OrdinaryDiffEqTsit5
import SoilDifferentialEquations: linear_interp, Update_SoilParam_Param!

"""
Config-driven test: Soil moisture simulation (Bonan) on CLM5 layers.
Usage:
  julia --project test/SM_uscrn/run_Bonan_CLM5_config.jl
  julia --project test/SM_uscrn/run_Bonan_CLM5_config.jl test/SM_uscrn/config_Bonan_CLM5.yaml
"""

# 1) Load config
cfg_file::String = isempty(ARGS) ? joinpath(@__DIR__, "config_Bonan_CLM5.yaml") : ARGS[1]
cfg::Dict{String, Any} = YAML.load_file(cfg_file)
println("Config: $cfg_file")

data_cfg = get(cfg, "data", Dict{String, Any}())
model_cfg = get(cfg, "model", Dict{String, Any}())
opt_cfg = get(cfg, "optimization", Dict{String, Any}())

# YAML may parse empty sections as `nothing`
(data_cfg === nothing) && (data_cfg = Dict{String, Any}())
(model_cfg === nothing) && (model_cfg = Dict{String, Any}())
(opt_cfg === nothing) && (opt_cfg = Dict{String, Any}())

# 2) Load observation data
file::String = data_cfg["file"]
if !isfile(file)
  file = joinpath(dirname(cfg_file), file)
end
obs = fread(file)

nsteps::Int = get(data_cfg, "time_steps", 0)
if nsteps > 0
  obs = obs[1:nsteps, :]
end

# columns: time, P_CALC, SM_5, SM_10, SM_20, SM_50, SM_100
obs_start_col::Int = get(data_cfg, "obs_start_col", 3)
θ_obs = obs[:, obs_start_col:end] |> Matrix
z_obs = Float64.(get(data_cfg, "depths_cm", [5, 10, 20, 50, 100])) ./ 100.0

# 3) CLM5 soil layers (center depths, meters)
layers = get_clm5_layers()
_z = layers.z
z₊ₕ = layers.z₊ₕ
Δz = layers.Δz
Δz₊ₕ = layers.Δz₊ₕ
N = length(Δz)

# CLM5 layer center depths (positive)
z_clm = abs.(_z[1:length(DZ_CLM5)])

# 4) Boundary and initial condition
surf_idx::Int = get(model_cfg, "surface_obs_index", 1)
θ_surf = θ_obs[:, surf_idx]
# initial profile from observed depths -> CLM5 layers
θ_clm = interp_obs_to_layer(θ_obs, z_obs, z_clm)
θ0 = θ_clm[1, :]

method_solve::String = get(model_cfg, "method_solve", "Bonan")
same_layer::Bool = get(model_cfg, "same_layer", true)
method_retention_cfg = get(model_cfg, "method_retention", "Campbell")

# Calibration uses observed depths (exclude surface by default)
ibeg = Int(get(model_cfg, "ibeg", 2))  # 2 -> exclude 5cm (surface)
yobs = θ_obs[:, ibeg:end]
set_option!(; yobs, θ_surf, ibeg=1, same_layer, method_solve)

function init_soil_clm5(; θ0, dt=3600.0, soil_type=7, method_retention="Campbell", same_layer=true)
  θ = fill(0.2, N)
  θ .= θ0

  par = get_soilpar(soil_type; method_retention)
  param = SoilParam(N, par; method_retention, same_layer)
  Soil{Float64}(; N, ibeg=1, dt, z=_z, z₊ₕ, Δz, Δz₊ₕ, θ, method_retention, param)
end

# 6) Run simulation
dt::Float64 = Float64(get(model_cfg, "dt", 3600.0))
soil_type::Int = get(model_cfg, "soil_type", 7)
method_default = method_retention_cfg isa AbstractVector ? method_retention_cfg[1] : method_retention_cfg
scheme::String = get(model_cfg, "layer_scheme", "full")
two_step::Bool = get(model_cfg, "two_step", false)

# soil for baseline simulation
soil0 = init_soil_clm5(; θ0, dt, soil_type, method_retention=method_default, same_layer)

# --- helpers for parameter schemes ---
# scheme = "full" | "two_layer" | "trend" | "exp_k" | "same_layer"
function expand_params(theta; scheme::String, mret::String, N::Int)
  if scheme == "two_layer"
    # split at ~30 cm: define top layers as z_clm <= 0.3 m
    ntop = count(z -> z <= 0.30, z_clm)
    nbot = N - ntop

    if mret == "Campbell"
      # [θ_sat_t, Ksat_t, ψ_sat_t, b_t, θ_sat_b, Ksat_b, ψ_sat_b, b_b]
      θt, Kt, ψt, bt, θb, Kb, ψb, bb = theta
      return [fill(θt, ntop); fill(θb, nbot)],
             [fill(Kt, ntop); fill(Kb, nbot)],
             [fill(ψt, ntop); fill(ψb, nbot)],
             [fill(bt, ntop); fill(bb, nbot)],
             nothing, nothing, nothing
    else
      # van_Genuchten: [θs_t, θr_t, K_t, α_t, n_t, θs_b, θr_b, K_b, α_b, n_b]
      θs_t, θr_t, Kt, αt, nt, θs_b, θr_b, Kb, αb, nb = theta
      return [fill(θs_t, ntop); fill(θs_b, nbot)],
             [fill(θr_t, ntop); fill(θr_b, nbot)],
             [fill(Kt, ntop); fill(Kb, nbot)],
             [fill(αt, ntop); fill(αb, nbot)],
             [fill(nt, ntop); fill(nb, nbot)],
             nothing, nothing
    end
  elseif scheme == "trend"
    # linear trend with depth: param = a + b*z (z in m, use center depth)
    z = z_clm
    if mret == "Campbell"
      # [θ0, k0, ψ0, b0, θ1, k1, ψ1, b1] where p = p0 + p1*z
      θ0, K0, ψ0, b0, θ1, K1, ψ1, b1 = theta
      return θ0 .+ θ1 .* z, K0 .+ K1 .* z, ψ0 .+ ψ1 .* z, b0 .+ b1 .* z, nothing, nothing, nothing
    else
      # van_Genuchten: [θs0, θr0, K0, α0, n0, θs1, θr1, K1, α1, n1]
      θs0, θr0, K0, α0, n0, θs1, θr1, K1, α1, n1 = theta
      return θs0 .+ θs1 .* z, θr0 .+ θr1 .* z, K0 .+ K1 .* z, α0 .+ α1 .* z, n0 .+ n1 .* z, nothing, nothing
    end
  else
    return nothing, nothing, nothing, nothing, nothing, nothing, nothing
  end
end

function apply_params!(soil, theta; scheme::String, mret::String)
  if scheme == "full"
    SM_UpdateParam!(soil, theta)
    return
  elseif scheme == "same_layer"
    # apply same-layer parameter vector to all layers
    if mret == "Campbell"
      θsat, Ksat, ψsat, b = theta
      soil.param.θ_sat .= θsat
      soil.param.Ksat .= Ksat
      soil.param.ψ_sat .= ψsat
      soil.param.b .= b
    else
      θsat, θres, Ksat, α, n = theta
      soil.param.θ_sat .= θsat
      soil.param.θ_res .= θres
      soil.param.Ksat .= Ksat
      soil.param.α .= α
      soil.param.n .= n
      soil.param.m .= 1 .- 1.0 ./ n
    end
    Update_SoilParam_Param!(soil.param)
    return
  elseif scheme == "exp_k"
    # exponential decay of Ksat with depth, others constant
    if mret == "Campbell"
      θsat, ψsat, b, K0, zK = theta
      soil.param.θ_sat .= θsat
      soil.param.ψ_sat .= ψsat
      soil.param.b .= b
      soil.param.Ksat .= K0 .* exp.(-z_clm ./ zK)
    else
      θsat, θres, α, n, K0, zK = theta
      soil.param.θ_sat .= θsat
      soil.param.θ_res .= θres
      soil.param.α .= α
      soil.param.n .= n
      soil.param.m .= 1 .- 1.0 ./ n
      soil.param.Ksat .= K0 .* exp.(-z_clm ./ zK)
    end
    Update_SoilParam_Param!(soil.param)
    return
  end

  # bounds for clamping (same-layer bounds)
  lower, upper = base_bounds(mret)

  if mret == "Campbell"
    θsat, Ksat, ψsat, b, _, _, _ = expand_params(theta; scheme, mret, N)
    θsat = clamp.(θsat, lower[1], upper[1])
    Ksat = clamp.(Ksat, lower[2], upper[2])
    ψsat = clamp.(ψsat, lower[3], upper[3])
    b    = clamp.(b,    lower[4], upper[4])

    soil.param.θ_sat .= θsat
    soil.param.Ksat .= Ksat
    soil.param.ψ_sat .= ψsat
    soil.param.b .= b
  else
    θsat, θres, Ksat, α, n, _, _ = expand_params(theta; scheme, mret, N)
    θsat = clamp.(θsat, lower[1], upper[1])
    θres = clamp.(θres, lower[2], upper[2])
    Ksat = clamp.(Ksat, lower[3], upper[3])
    α    = clamp.(α,    lower[4], upper[4])
    n    = clamp.(n,    lower[5], upper[5])

    soil.param.θ_sat .= θsat
    soil.param.θ_res .= θres
    soil.param.Ksat .= Ksat
    soil.param.α .= α
    soil.param.n .= n
    soil.param.m .= 1 .- 1.0 ./ n
  end
  Update_SoilParam_Param!(soil.param)
end

function base_theta0(soil, mret::String)
  if mret == "Campbell"
    return [soil.param.θ_sat[1], soil.param.Ksat[1], soil.param.ψ_sat[1], soil.param.b[1]]
  else
    return [soil.param.θ_sat[1], soil.param.θ_res[1], soil.param.Ksat[1], soil.param.α[1], soil.param.n[1]]
  end
end

function base_bounds(mret::String)
  # use same_layer=true to get base parameter bounds
  soil_tmp = init_soil_clm5(; θ0, dt, soil_type, method_retention=mret, same_layer=true)
  lower, upper = SM_paramBound(soil_tmp)
  return lower, upper
end

function bounds_for_scheme(mret::String; scheme::String)
  lower, upper = base_bounds(mret)
  if scheme == "full"
    # full uses layer-specific bounds from same_layer=false soil
    soil_full = init_soil_clm5(; θ0, dt, soil_type, method_retention=mret, same_layer)
    return SM_paramBound(soil_full)
  elseif scheme == "two_layer"
    return vcat(lower, lower), vcat(upper, upper)
  elseif scheme == "trend"
    # use same bounds for intercepts; slopes constrained to small range
    l = vcat(lower, fill(-0.5, length(lower)))
    u = vcat(upper, fill(0.5, length(upper)))
    return l, u
  elseif scheme == "exp_k"
    # bounds: base params + K0 + zK
    if mret == "Campbell"
      # [θsat, Ksat, ψsat, b] -> use θsat, ψsat, b; K0, zK
      l = [lower[1], lower[3], lower[4], lower[2], 0.05]
      u = [upper[1], upper[3], upper[4], upper[2], 2.00]
      return l, u
    else
      # [θsat, θres, Ksat, α, n] -> use θsat, θres, α, n; K0, zK
      l = [lower[1], lower[2], lower[4], lower[5], lower[3], 0.05]
      u = [upper[1], upper[2], upper[4], upper[5], upper[3], 2.00]
      return l, u
    end
  end
end

function two_layer_to_full(theta; mret::String)
  ntop = count(z -> z <= 0.30, z_clm)
  nbot = N - ntop
  if mret == "Campbell"
    θt, Kt, ψt, bt, θb, Kb, ψb, bb = theta
    θsat = [fill(θt, ntop); fill(θb, nbot)]
    Ksat = [fill(Kt, ntop); fill(Kb, nbot)]
    ψsat = [fill(ψt, ntop); fill(ψb, nbot)]
    b    = [fill(bt, ntop); fill(bb, nbot)]
    return vcat(θsat, Ksat, ψsat, b)
  else
    θs_t, θr_t, Kt, αt, nt, θs_b, θr_b, Kb, αb, nb = theta
    θsat = [fill(θs_t, ntop); fill(θs_b, nbot)]
    θres = [fill(θr_t, ntop); fill(θr_b, nbot)]
    Ksat = [fill(Kt, ntop); fill(Kb, nbot)]
    α    = [fill(αt, ntop); fill(αb, nbot)]
    n    = [fill(nt, ntop); fill(nb, nbot)]
    return vcat(θsat, θres, Ksat, α, n)
  end
end

function model_sim(theta; mret=method_default, scheme=scheme)
  soil = init_soil_clm5(; θ0, dt, soil_type, method_retention=mret, same_layer)
  apply_params!(soil, theta; scheme, mret)
  if method_solve == "Bonan"
    return solve_SM_Bonan(soil, θ_surf)
  else
    return solve_SM_ODE(soil, θ_surf; solver=Tsit5())
  end
end

function goal(theta; mret=method_default, scheme=scheme)
  ysim = model_sim(theta; mret, scheme)

  # interpolate CLM5-layer simulation back to observed depths
  ncol = size(yobs, 2)
  loss = 0.0
  for i in 1:ncol
    obs = @view yobs[:, i]
    # observed depth (m), skip surface by ibeg
    z_tgt = z_obs[ibeg - 1 + i]
    sim = [linear_interp(z_clm, collect(@view(ysim[t, 1:length(z_clm)])), [z_tgt])[1] for t in axes(ysim, 1)]
    loss += -GOF(obs, sim).NSE
  end
  return loss / ncol
end

# optional optimization (SCE-UA)
if get(opt_cfg, "enable", false)
  maxn::Int = get(opt_cfg, "maxn", 10_000)
  method_list = method_retention_cfg isa AbstractVector ? method_retention_cfg : [method_retention_cfg]
  for mret in method_list
    set_option!(; method_retention=mret)
    soil_opt = init_soil_clm5(; θ0, dt, soil_type, method_retention=mret, same_layer)

    # initial guess and bounds by scheme
    if scheme == "full"
      theta0 = SM_param2theta(soil_opt)
      lower, upper = SM_paramBound(soil_opt)
    elseif scheme == "two_layer"
      base = base_theta0(soil_opt, mret)
      theta0 = vcat(base, base)
      lower, upper = bounds_for_scheme(mret; scheme)
    elseif scheme == "trend"
      base = base_theta0(soil_opt, mret)
      theta0 = vcat(base, zeros(length(base)))
      lower, upper = bounds_for_scheme(mret; scheme)
    elseif scheme == "exp_k"
      lower, upper = bounds_for_scheme(mret; scheme)
      if mret == "Campbell"
        # [θsat, ψsat, b, K0, zK]
        θsat, Ksat, ψsat, b = base_theta0(soil_opt, mret)
        theta0 = [θsat, ψsat, b, Ksat, 0.3]
      else
        # [θsat, θres, α, n, K0, zK]
        θsat, θres, Ksat, α, n = base_theta0(soil_opt, mret)
        theta0 = [θsat, θres, α, n, Ksat, 0.3]
      end
    end

    # two-step warm start: optimize same_layer first, then expand
    if two_step
      println("[Step1] Optimizing same_layer (warm start) $mret ...")
      soil_step1 = init_soil_clm5(; θ0, dt, soil_type, method_retention=mret, same_layer=true)
      θ0_step1 = base_theta0(soil_step1, mret)
      lower1, upper1 = base_bounds(mret)
      θopt_step1, feval1, exitflag1 = sceua(θ -> goal(θ; mret, scheme="same_layer"), θ0_step1, lower1, upper1; maxn)
      println("[Step1] done ($mret). feval=$feval1 exit=$exitflag1")

      # expand step1 params as initial guess for step2
      if scheme == "two_layer"
        theta0 = vcat(θopt_step1, θopt_step1)
      elseif scheme == "trend"
        theta0 = vcat(θopt_step1, zeros(length(θopt_step1)))
      elseif scheme == "exp_k"
        if mret == "Campbell"
          θsat, Ksat, ψsat, b = θopt_step1
          theta0 = [θsat, ψsat, b, Ksat, 0.3]
        else
          θsat, θres, Ksat, α, n = θopt_step1
          theta0 = [θsat, θres, α, n, Ksat, 0.3]
        end
      elseif scheme == "full"
        # multi-step: use two_layer as intermediate, then full (see below)
        theta0 = SM_param2theta(soil_opt)
      end
    end

    println("Optimizing (SCE-UA) $mret, maxn=$maxn, scheme=$scheme ...")
    @time theta_opt, feval, exitflag = sceua(theta -> goal(theta; mret, scheme), theta0, lower, upper; maxn)
    println("Optimization done ($mret). feval=$feval exit=$exitflag")

    # multi-step refinement: two_layers -> full (optional)
    if two_step && scheme == "two_layer"
      println("[Step3] Refinement two_layer -> full ($mret) ...")
      soil_full = init_soil_clm5(; θ0, dt, soil_type, method_retention=mret, same_layer=false)
      theta0_full = two_layer_to_full(theta_opt; mret)
      lower_full, upper_full = SM_paramBound(soil_full)
      @time theta_full, feval_full, exitflag_full = sceua(theta -> goal(theta; mret, scheme="full"), theta0_full, lower_full, upper_full; maxn)
      println("[Step3] done ($mret). feval=$feval_full exit=$exitflag_full")
    end
  end
end

ysim = solve_SM_Bonan(soil0, θ_surf)

@testset "CLM5 Bonan quick case (config)" begin
  @test size(ysim, 1) == size(θ_obs, 1)
  @test size(ysim, 2) == N
  @test !any(isnan.(ysim))
end
