"""
    soil_depth_init(Δz::AbstractVector)
    
Soil depth initialization

```julia
z, z₊ₕ, dz₊ₕ = soil_depth_init(Δz)
```
"""
function soil_depth_init(Δz::AbstractVector)
  # Soil depth (m) at i+1/2 interface between layers i and i+1 (negative distance from surface)
  # z_{i+1/2}
  N = length(Δz)

  z = zeros(N)
  z₊ₕ = zeros(N)
  dz₊ₕ = zeros(N)

  z₊ₕ[1] = -Δz[1]
  for i = 2:N
    z₊ₕ[i] = z₊ₕ[i-1] - Δz[i] # on the edge
  end

  # Soil depth (m) at center of layer i (negative distance from surface)
  z[1] = 0.5 * z₊ₕ[1]
  for i = 2:N
    z[i] = 0.5 * (z₊ₕ[i-1] + z₊ₕ[i]) # on the center
  end

  # Thickness between between z(i) and z(i+1)
  for i = 1:N-1
    dz₊ₕ[i] = z[i] - z[i+1]
  end
  dz₊ₕ[N] = 0.5 * Δz[N]

  (; z, z₊ₕ, dz₊ₕ)
end


function cal_Δz(z)
  N = length(z)
  z₊ₕ = zeros(N)
  Δz = zeros(N)
  Δz[1] = 0 - z[1] * 2
  z₊ₕ[1] = -Δz[1]

  for i in 2:N
    Δz[i] = (z₊ₕ[i-1] - z[i]) * 2
    z₊ₕ[i] = z₊ₕ[i-1] - Δz[i]
  end
  Δz
end


# "face to center"
# function C2F(z::AbstractVector)
#   N = length(z)
#   z₊ₕ = zeros(N)
#   d = zeros(N)

#   z₊ₕ[1] = z[1] * 2
#   d[1] = z[1] * 2

#   @inbounds for i = 2:N
#     d[i] = 2(z[i] - z[i-1]) - d[i-1]
#     z₊ₕ[i] = z₊ₕ[i-1] + d[i]
#   end
#   z₊ₕ
# end

# "center to face"
# function F2C(z₊ₕ::AbstractVector)
#   N = length(z₊ₕ)
#   z = zeros(N)
#   z[1] = z₊ₕ[1] / 2
#   @inbounds for i = 2:N
#     z[i] = 0.5 * (z₊ₕ[i] + z₊ₕ[i-1])
#   end
#   z
# end
# export cal_Δz₊ₕ, C2F, F2C
export cal_Δz

# Δz = [2, 4, 6, 10]
# z, z₊ₕ, Δz₊ₕ = soil_depth_init(Δz)
# cal_Δz(z₊ₕ) == Δz
# cal_Δz₊ₕ(z, z₊ₕ) == Δz₊ₕ
