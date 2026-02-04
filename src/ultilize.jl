export interp_data_depths
using Ipaper: approx


function interp_data_depths(A::M, z::V, zout::V) where {
    T<:Real,M<:AbstractMatrix{T},V<:AbstractVector{T}}

    ntime = size(A, 1)
    yout = zeros(ntime, length(zout))
    @inbounds for i in 1:ntime
        yout[i, :] .= approx(z, view(A, i, :), zout)
    end
    yout
end


function init_grid(zs_center::Vector{FT}) where {FT<:Real}
  z₊ₕ = center_to_face(zs_center)
  Δz = face_to_thickness(z₊ₕ) ./ 100.0  # [cm] -> [m]
  z, z₋ₕ, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
  N = length(Δz)
  zs_sim = abs.(z[1:N]) .* 100  # [m] -> [cm]
  return (; z, z₊ₕ, Δz, Δz₊ₕ, N, zs_sim)
end

# 查找层索引: 根据边界深度找到对应的模拟层索引
function find_layer_indices(z_bound_top::T, zs_sim::V, zs_obs=nothing) where {T<:Real,V<:AbstractVector{T}}
  isnothing(zs_obs) && (zs_obs = zs_sim)

  ibeg = findfirst(==(z_bound_top), abs.(zs_sim)) + 1 # 模拟开始层
  itop = findfirst(==(z_bound_top), abs.(zs_obs))     # 观测层中上边界层索引
  return (; ibeg, itop)
end


function guess_outdir(config::Config, outdir=nothing)
  isnothing(outdir) ? joinpath(dirname(config.file), "outdir") : outdir
end

function open_log(config::Config, log_file=nothing)
  isnothing(log_file) && (log_file = replace(config.fileConfig, r"\.yaml$" => ".log"))
  io = open(log_file, "a")
  return io
end

function log(io, msg)
  println(msg)
  if !isnothing(io)
    println(io, msg)
    flush(io)
  end
end
