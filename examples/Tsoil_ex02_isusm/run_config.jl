using YAML, Random, Printf
using SoilDifferentialEquations, Plots, Test, RTableTools, Dates
using OrdinaryDiffEqTsit5
import ModelParams: sceua, GOF, of_KGE, of_NSE
import NetCDFTools: approx

# --- Helpers ---
function get_clm_layers(n::Int=10; fs=0.025)
    diff([fs * (exp(0.5 * i) - 1) for i in 0:n])
end

function interp_profile(z_grid, profile, z_tgt)
    z_tgt <= z_grid[1] && return profile[1]
    z_tgt >= z_grid[end] && return profile[end]
    for i in 1:length(z_grid)-1
        z1, z2 = z_grid[i], z_grid[i+1]
        z1 <= z_tgt <= z2 && return profile[i] + (z_tgt - z1)/(z2 - z1)*(profile[i+1] - profile[i])
    end
    NaN
end

function solve_full(soil, Tsurf; solver, kw...)
    (; N, ibeg, dt) = soil
    prob = SoilDifferentialEquations._ODEProblem((dT, T, p, t) -> TsoilEquation_partial(dT, T, p, t; ibeg), 
        soil.Tsoil[ibeg:N], (0, dt), soil)
    
    mapreduce(vcat, 1:length(Tsurf)) do i
        if i > 1
            soil.Tsurf = Tsurf[i]
            prob.u0 .= soil.Tsoil[ibeg:N]
            sol = SoilDifferentialEquations._solve(prob, solver; saveat=dt, kw...)
            soil.Tsoil[ibeg:N] .= sol.u[end]
        end
        permutedims(soil.Tsoil[ibeg:N])
    end
end

function model_sim(soil, Tsurf, theta; solver)
    soil.param.κ .= theta[1:length(theta)÷2]
    soil.param.cv .= theta[length(theta)÷2+1:end]
    solve_full(soil, Tsurf; solver)
end

function init_soil_custom(cfg, Tsoil0, Δz, ibeg_grid, inds_obs)
    z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
    m_sat = θ_S[cfg["model"]["soil_type"]] * 1000.0 * Δz # ρ_wat=1000
    κ, cv = soil_properties_thermal(Δz, deepcopy(Tsoil0), 0.8*m_sat, 0*m_sat; soil_type=cfg["model"]["soil_type"])
    Soil{Float64}(; N=length(Δz), ibeg=ibeg_grid, inds_obs=inds_obs, dt=Float64(cfg["simulation"]["dt"]), 
        z, z₊ₕ, Δz, Δz₊ₕ, κ, cv, Tsurf=0.0, Tsoil=deepcopy(Tsoil0))
end

# --- Main ---
function main()
    cfg_file = isempty(ARGS) ? joinpath(@__DIR__, "config.yaml") : ARGS[1]
    cfg = YAML.load_file(cfg_file)
    println("Config: $cfg_file | Data: $(cfg["data"]["file"])")

    # Data
    d = fread(isfile(cfg["data"]["file"]) ? cfg["data"]["file"] : joinpath(dirname(cfg_file), cfg["data"]["file"]))
    t = d[1:cfg["data"]["time_steps"], cfg["data"]["time_col"]]
    yobs_full = Matrix(d[1:cfg["data"]["time_steps"], cfg["data"]["obs_start_col"]:end])

    # Grid
    Δz = haskey(cfg["model"], "layers_thickness") ? Float64.(cfg["model"]["layers_thickness"]) :
         get(cfg["model"], "use_clm_formula", false) ? get_clm_layers(get(cfg["model"], "n_layers", 10)) :
         fill(cfg["model"]["dz"], ceil(Int, 1.35/cfg["model"]["dz"]))
    
    z_nodes = cumsum(Δz) .- 0.5Δz
    z_obs = haskey(cfg["model"], "depths_cm") ? Float64.(cfg["model"]["depths_cm"]) ./ 100 :
            Float64.(cfg["model"]["depths_inch"]) .* 0.0254

    # Map Obs to Grid
    inds_obs = [argmin(abs.(z_nodes .- zo)) for zo in z_obs]
    println("Grid: $(length(Δz)) layers. Obs mapped to layers: $inds_obs")

    # Boundary & Initial
    ibeg_obs = cfg["model"]["Tsurf_layer_index"]
    ibeg_grid = inds_obs[ibeg_obs]
    Tsurf = yobs_full[:, ibeg_obs]
    Tsoil0 = approx(inds_obs, yobs_full[1, :], 1:length(Δz))

    # Optimization Setup
    yobs_tgt = yobs_full[:, ibeg_obs:end]
    obs_flat = yobs_tgt[:, 2:end][:]
    z_tgt = z_obs[ibeg_obs+1:end]
    solver = Tsit5()
    
    soil = init_soil_custom(cfg, Tsoil0, Δz, ibeg_grid, inds_obs[ibeg_obs:end])

    if cfg["optimization"]["enable"]
        println("Optimizing...")
        Random.seed!(1)
        f(θ) = begin
            R = model_sim(init_soil_custom(cfg, Tsoil0, Δz, ibeg_grid, inds_obs[ibeg_obs:end]), Tsurf, θ; solver)
            sim = [interp_profile(soil.z[ibeg_grid:end], R[i,:], z) for i in 1:size(R,1), z in z_tgt]
            -GOF(obs_flat, sim[:]).NSE
        end
        
        n = length(Δz)
        θ_opt, _, nse = sceua(f, [soil.κ; soil.cv], 
            [fill(cfg["optimization"]["bounds"]["kappa"][1], n); fill(cfg["optimization"]["bounds"]["cv"][1], n)],
            [fill(cfg["optimization"]["bounds"]["kappa"][2], n); fill(cfg["optimization"]["bounds"]["cv"][2], n)];
            maxn=cfg["optimization"]["max_iterations"])
        println("NSE: $(-nse)")
        soil.param.κ .= θ_opt[1:n]; soil.param.cv .= θ_opt[n+1:end]
    end

    # Final Run & Plot
    R = solve_full(soil, Tsurf; solver)
    if cfg["output"]["plot"]
        plts = map(1:size(yobs_tgt, 2)) do i
            z_t = z_obs[ibeg_obs + i - 1]
            sim = [interp_profile(soil.z[ibeg_grid:end], R[t,:], z_t) for t in 1:size(R,1)]
            plot(t, [yobs_tgt[:, i] sim], label=["Obs" "Sim"], title=@sprintf("Depth: %.1f cm", z_t*100))
        end
        savefig(plot(plts..., layout=(ceil(Int, length(plts)/3), 3), size=Tuple(cfg["output"]["plot_size"])), 
            cfg["output"]["save_fig"])
        println("Saved: $(cfg["output"]["save_fig"])")
    end
end

main()
