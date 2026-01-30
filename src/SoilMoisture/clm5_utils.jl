
export DZ_CLM5, get_clm5_layers, interp_obs_to_layer

# CLM5 default layer thickness (m)
# Reference: Zeng & Decker (2009), also used in SoilTemperature/soil_temperature_F0.jl
const DZ_CLM5 = [0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360, 0.5539, 0.9133, 1.5058]

"""
    get_clm5_layers()

Get the standard CLM5 soil layer depths and thicknesses.
Returns a NamedTuple with `z`, `z₋ₕ`, `z₊ₕ`, `Δz`, `Δz₊ₕ`.
"""
function get_clm5_layers()
  Δz = DZ_CLM5
  res = soil_depth_init(Δz)
  merge(res, (; Δz))
end

"""
    interp_obs_to_layer(yobs::AbstractMatrix, z_obs::AbstractVector, z_layer::AbstractVector)

Interpolate observed data `yobs` (time x depth) from observed depths `z_obs` to model layer depths `z_layer`.
Assumes linear interpolation.
"""
function interp_obs_to_layer(yobs::AbstractMatrix, z_obs::AbstractVector, z_layer::AbstractVector)
  # yobs: [time, depth]
  ntime = size(yobs, 1)
  nlayer = length(z_layer)
  y_interp = zeros(eltype(yobs), ntime, nlayer)

  # Pre-sort z_obs just in case
  p = sortperm(z_obs)
  z_obs_sorted = z_obs[p]
  
  # Use absolute depths for interpolation (assume z is negative or positive, but distance matters)
  z_obs_abs = abs.(z_obs_sorted)
  z_layer_abs = abs.(z_layer)

  for t in 1:ntime
    y_t = yobs[t, p]
    y_interp[t, :] = linear_interp(z_obs_abs, y_t, z_layer_abs)
  end
  return y_interp
end

function linear_interp(x, y, xi)
  yi = similar(xi)
  for k in eachindex(xi)
    val = xi[k]
    if val <= x[1]
      yi[k] = y[1]
    elseif val >= x[end]
      yi[k] = y[end]
    else
      i = searchsortedlast(x, val)
      if i == length(x)
        yi[k] = y[end]
      else
        x1, x2 = x[i], x[i+1]
        y1, y2 = y[i], y[i+1]
        yi[k] = y1 + (val - x1) * (y2 - y1) / (x2 - x1)
      end
    end
  end
  return yi
end
