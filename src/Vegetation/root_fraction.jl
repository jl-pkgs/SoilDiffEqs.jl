export root_fraction, root_fraction_sum
export root_dist_par

# 通过WUD，推求β
function root_dist_par(WUD_cm::T=200.0; p=0.95) where {T<:Real}
  # β^(WUD) = 1-p # 0.95
  # WUD ln(β) = ln(1-p)
  exp(log(1 - p) / (WUD_cm))
end

"""
- z: in cm, 向下为正
"""
root_fraction_sum(z_cm::AbstractVector; β=0.94) = 1 .- β .^ z_cm

# function root_fraction(soil::Soil; β=0.94)
#   N = soil.N
#   z₊ₕ = -[0; soil.z₊ₕ[1:N]] .* 100
#   cum = root_fraction(z₊ₕ; β)
#   -diff(cum)
# end
function root_fraction(z₊ₕ_cm::AbstractVector; β=0.94)
  N = length(z₊ₕ_cm)
  f = zeros(N)
  # z₊ₕ = -soil.z₊ₕ[1:N] .* 100 # to cm
  for i in 2:N
    f[i] = β^z₊ₕ_cm[i-1] - β^z₊ₕ_cm[i]
  end
  f[1] = 1 - β^z₊ₕ_cm[1]
  f
end

function root_fraction(soil::Soil; β=0.94)
  N = soil.N
  z₊ₕ = -soil.z₊ₕ[1:N] .* 100 # to cm
  root_fraction(z₊ₕ; β)
end


function Root_Water_Uptake(soil::Soil, T::Float64, E::Float64)
  soil.sink[1] = E # convert [mm/h] to [cm/h]
  for i in 2:p.n_layer
    p.sink[i] = T * fraction[i] # 加权之后的fraction, root, θ
    ## 将ET分配给每一层的土壤
  end
end
