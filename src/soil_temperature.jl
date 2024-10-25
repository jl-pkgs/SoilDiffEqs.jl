"""
    soil_temperature!(soil::Soil, Tsurf_next::Real;
        solution="implicit", method="apparent-heat-capacity")

> 这里采用的是方案1的边界条件, sum(F) = G, G = - k1 * (T1 - T0) / dz1

# Arguments: 
- `κ`          : thermal conductivity, [W/m/K]
- `cv`         : volumetric heat capacity, [J/m3/K]
- `Tsoil`      : Tsoil_cur
- `Tsurf_next` : Tsurf_next_next, T0_{n+1}
- `solution`   : 
  + `implicit`       : 
  + `crank-nicolson` : 

# Examples
```julia
Tsoil_next, G = soil_temperature!(Soil, Real; solution="implicit", method="apparent-heat-capacity")
```
"""
function soil_temperature!(soil::Soil, Tsurf_next::Real;
  solution="implicit", method="apparent-heat-capacity", ibeg::Int=1)
  (; n, dt, Δz, z, z₊ₕ, Δz₊ₕ,
    κ, cv, Tsoil, 
    κ₊ₕ, u, a, b, c, d, f) = soil

  # Thermal conductivity at interface (W/m/K)
  # κ₊ₕ = zeros(1, n - 1)
  @inbounds for i = 1:n-1
    κ₊ₕ[i] = κ[i] * κ[i+1] * (z[i] - z[i+1]) /
             (κ[i] * (z₊ₕ[i] - z[i+1]) + κ[i+1] * (z[i] - z₊ₕ[i])) # Eq. 5.16
  end

  ## Set up tridiagonal matrix
  @inbounds if solution == "implicit"
    for i = ibeg:n
      m = cv[i] * Δz[i] / dt
      if i == ibeg
        a[i] = 0
        c[i] = -κ₊ₕ[i] / Δz₊ₕ[i]
        b[i] = m - c[i] + κ[i] / (0 - z[i])  # κ_(1/2) / dz_(1/2) = κ₁ / (0 - z₁)
        d[i] = m * Tsoil[i] + κ[i] / (0 - z[i]) * Tsurf_next

      elseif i < n
        a[i] = -κ₊ₕ[i-1] / Δz₊ₕ[i-1]
        c[i] = -κ₊ₕ[i] / Δz₊ₕ[i]
        b[i] = m - a[i] - c[i]
        d[i] = m * Tsoil[i]

      elseif i == n
        a[i] = -κ₊ₕ[i-1] / Δz₊ₕ[i-1]
        c[i] = 0
        b[i] = m - a[i]
        d[i] = m * Tsoil[i]
      end
    end

  elseif solution == "crank-nicolson"
    # --- Heat flux at time n (W/m2) of each layer
    # f = zeros(n)
    for i = ibeg:n-1
      f[i] = -κ₊ₕ[i] / Δz₊ₕ[i] * (Tsoil[i] - Tsoil[i+1]) # Eq. 5.15
    end

    for i = ibeg:n
      m = cv[i] * Δz[i] / dt
      if i == ibeg
        a[i] = 0
        c[i] = -0.5 * κ₊ₕ[i] / Δz₊ₕ[i]
        b[i] = m - c[i] + κ[i] / (0 - z[i])
        d[i] = m * Tsoil[i] + κ[i] / (0 - z[i]) * Tsurf_next + 0.5 * f[i]

      elseif i < n
        a[i] = -0.5 * κ₊ₕ[i-1] / Δz₊ₕ[i-1]
        c[i] = -0.5 * κ₊ₕ[i] / Δz₊ₕ[i]
        b[i] = m - a[i] - c[i]
        d[i] = m * Tsoil[i] + 0.5 * (f[i] - f[i-1])

      elseif i == n
        a[i] = -0.5 * κ₊ₕ[i-1] / Δz₊ₕ[i-1]
        c[i] = 0
        b[i] = m - a[i]
        d[i] = m * Tsoil[i] - 0.5 * f[i-1]
      end
    end
  end

  _a = @view a[ibeg:end]
  _b = @view b[ibeg:end]
  _c = @view c[ibeg:end]
  _d = @view d[ibeg:end]
  u[ibeg:end] .= tridiagonal_solver(_a, _b, _c, _d) # the updated Tsoil

  # --- Derive energy flux into soil (W/m2)
  z_above = ibeg == 1 ? 0 : z[ibeg-1]
  G = κ[ibeg] * (Tsurf_next - u[ibeg]) / (z_above - z[ibeg]) # gound heat flux, G

  ## --- Check for energy conservation
  if method == "apparent-heat-capacity"
    LE_f = 0.0
    # elseif method == "excess-heat"
    #   phase_change(physcon, soilvar, dt);
  end

  edif = 0.0 # Sum change in energy (W/m2)
  @inbounds for i = ibeg:n
    edif += cv[i] * Δz[i] * (u[i] - Tsoil[i]) / dt
  end

  # Error check
  err = edif - G - LE_f
  # abs(err) > 1e-03 && error("Soil temp erature energy conservation error")

  Tsoil .= u # update Tsoil
  soil.G = G
  u, G
end
