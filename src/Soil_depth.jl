function cal_Δz(z)
  n = length(z)
  z₊ₕ = zeros(n)
  Δz = zeros(n)
  Δz[1] = 0 - z[1] * 2
  z₊ₕ[1] = -Δz[1]

  for i in 2:n
    Δz[i] = (z₊ₕ[i-1] - z[i]) * 2
    z₊ₕ[i] = z₊ₕ[i-1] - Δz[i]
  end
  Δz
end


# "face to center"
# function C2F(z::AbstractVector)
#   n = length(z)
#   z₊ₕ = zeros(n)
#   d = zeros(n)

#   z₊ₕ[1] = z[1] * 2
#   d[1] = z[1] * 2

#   @inbounds for i = 2:n
#     d[i] = 2(z[i] - z[i-1]) - d[i-1]
#     z₊ₕ[i] = z₊ₕ[i-1] + d[i]
#   end
#   z₊ₕ
# end

# "center to face"
# function F2C(z₊ₕ::AbstractVector)
#   n = length(z₊ₕ)
#   z = zeros(n)
#   z[1] = z₊ₕ[1] / 2
#   @inbounds for i = 2:n
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
