# Richards Equation Implementation Analysis

File: `src/SoilMoisture/Equation_Richards.jl`

## Main Functions

### 1. **`cal_Q!`** (Lines 3-34) - Calculate Water Flux
Computes the water flux `Q` between soil layers using Darcy's law:
- **Two boundary condition methods**:
  - `ψ0`: First type boundary (prescribed matric potential at surface)
  - `Q0`: Second type boundary (prescribed flux at surface)
- **Flux calculation** (line 29):
  ```julia
  Q[i] = -K₊ₕ[i] * ((ψ[i] - ψ[i+1]) / Δz₊ₕ[i] + 1.0)
  ```
  This is Darcy's law: Q = -K * (dψ/dz + 1), where the "+1" accounts for gravitational potential
- **Bottom boundary** (line 31): Gravity drainage (`Q[N] = -K[N]`)

### 2. **`soil_Updateθ!`** (Lines 38-53) - Update Soil Moisture
Updates soil water content θ using a finite difference mass balance:
- Converts time from seconds to hours (line 41)
- **Mass balance equation** (lines 47-50):
  ```julia
  θ[i] += ((-Q[i-1] + Q[i]) - sink[i]) * dt / Δz[i]
  ```
  Change in θ = (inflow - outflow - sink) * dt / layer_thickness
- Clamps θ to valid physical ranges

### 3. **`clamp_θ!`** (Lines 55-73) - Constrain θ Values
Two versions for different soil models:
- **Van Genuchten**: Bounds θ between `θ_res + δ` and `θ_sat` (lines 55-63)
  - `δ = (θ_sat - θ_res) * 0.01` provides a small buffer above residual water content
- **Campbell**: Bounds θ between 0.01 and `θ_sat` (lines 65-73)
  - Uses a fixed minimum `_θ_res = 0.01` instead of parameter-based residual content

### 4. **`RichardsEquation`** (Lines 83-95) - ODE Right-Hand Side
Main function for ODE solvers (DifferentialEquations.jl):
- Calculates `dθ/dt` for each soil layer
- Divides by 3600 to convert from hourly to per-second rates (line 91)
- Increments timestep counter (line 84)

### 5. **`RichardsEquation_partial`** (Lines 98-105)
Wrapper for solving only active soil layers (`ibeg:N`) instead of full profile.

## Key Physics

- **Conservation of mass**: ∂θ/∂t = -∂Q/∂z - S(z,t)
- **Darcy's law**: Q = -K(θ) * (∂ψ/∂z + 1)
- **Gravity drainage** at bottom boundary
- **Sink term** represents water extraction (e.g., root uptake)

## Units

- θ: [m³ m⁻³] volumetric water content
- Q: [cm h⁻¹] water flux
- z: [cm] depth
- dt: converts between [s] and [h]

## Implementation Details

### Boundary Conditions

**Top boundary (surface)**:
- Method "ψ0": Prescribed matric potential (line 25)
  ```julia
  Q0 = -_K₊ₕ * ((ψ0 - ψ[ibeg]) / _dz + 1)
  ```
- Method "Q0": Prescribed flux (direct assignment)

**Bottom boundary**:
- Free drainage: `Q[N] = -K[N]` (gravity-driven flow only)

### Numerical Scheme

The implementation uses:
- Finite difference method for spatial discretization
- Variable layer thickness (Δz)
- Staggered grid: K₊ₕ calculated at layer interfaces
- Explicit/implicit time integration via DifferentialEquations.jl

### Active Layer Management

- `ibeg`: First active layer index (allows for snow/surface layers)
- `N`: Total number of soil layers
- Only layers `ibeg:N` are actively solved

---

# Soil Data Structure

File: `src/Soil.jl`

## Structure Overview

The `Soil` struct (lines 3-78) is a mutable container holding all state variables, parameters, and temporary variables for soil modeling. It uses the `@with_kw_noshow` macro for keyword argument construction.

## Key Components

### 1. Grid Configuration (Lines 4-16)
```julia
N::Int = 10                    # Number of soil layers
ibeg::Int = 1                  # First active layer (allows surface layers)
inds_obs::Vector{Int}          # Indices of observed layers
z::OffsetVector                # Depth at layer center [m] (negative downward)
Δz::Vector                     # Layer thickness [m]
z₊ₕ, Δz₊ₕ                      # Depths/thickness at layer interfaces
```
- **Staggered grid**: Variables at centers (`z`), interfaces (`z₊ₕ`)
- **Dual units**: Both `[m]` and `[cm]` versions maintained (e.g., `z_cm = z * 100`)

### 2. Soil Moisture Variables (Lines 18-31)
```julia
θ::Vector                      # Volumetric water content [m³/m³]
ψ::Vector                      # Matric potential [cm] (negative when dry)
K::Vector                      # Hydraulic conductivity [cm/h]
K₊ₕ::Vector                    # K at interfaces [cm/h]
Q::Vector                      # Water flux [cm/h]
∂θ∂ψ::Vector                   # Specific moisture capacity dθ/dψ [cm⁻¹]
sink::Vector                   # Water extraction [cm/time]
```
- **Boundary values**: `θ0`, `ψ0`, `Q0` at surface
- **Backup arrays**: `θ_prev`, `ψ_prev` for time stepping
- **Note**: ψ becomes more negative as soil dries

### 3. Equilibrium States (Lines 33-34)
```julia
θE::Vector                     # Equilibrium water content [m³/m³]
ψE::Vector                     # Equilibrium potential [cm]
```
Used for groundwater equilibrium calculations.

### 4. Groundwater Module (Lines 36-45)
```julia
zwt::FT                        # Water table depth [m]
jwt::Int                       # Non-saturated zone index
wa::FT = 5000.0                # Aquifer storage [mm]
uex::FT                        # Surface excess water [mm]
recharge::FT                   # Recharge rate [cm/h]
drainage::FT                   # Drainage rate [cm/h]
Sy_r, Sy_d, Sy_e::Vector      # Specific yield [m³/m³]
```
- **Sy_r**: Specific yield for recharge (water table rising)
- **Sy_d**: Specific yield for discharge (water table falling)
- **Sy_e**: Equilibrium specific yield
- **uex**: Excess water at surface becomes runoff

### 5. Soil Temperature Variables (Lines 47-53)
```julia
Tsoil::Vector                  # Soil temperature [°C]
κ₊ₕ::Vector                    # Thermal conductivity at interfaces [W/m/K]
F::Vector                      # Heat flux [W/m²]
Tsurf::FT                      # Surface temperature [°C]
F0::FT                         # Surface heat flux [W/m²] (negative downward)
G::FT                          # Soil heat flux [W/m²]
```

### 6. Vegetation Module (Lines 55-59)
```julia
f_root::Vector                 # Root distribution fraction
βi::Vector                     # Water stress per layer (dryness)
w::Vector                      # Weighted stress (f_root * βi)
β::FT                          # Total water stress
```
Implements root water uptake with stress functions.

### 7. Soil Parameters (Lines 61-63)
```julia
method_retention::String       # "van_Genuchten" or "Campbell"
param::SoilParam               # Hydraulic + thermal parameters
```
Links to parameter structs (Van Genuchten or Campbell hydraulic models).

### 8. Numerical Solver Variables (Lines 65-77)
```julia
u, du::Vector                  # ODE state vectors
a, b, c, d, e, f::Vector      # Tridiagonal matrix solver temporaries
timestep::Int                  # Iteration counter
```

## Constructor Functions

### `Soil{FT}()` (Lines 80-87)
Creates soil with specified retention curve method:
```julia
soil = Soil{Float64}(; method_retention="van_Genuchten", N=10, ...)
```
- `"van_Genuchten"` → uses `VanGenuchten{FT}`
- `"Campbell"` → uses `Campbell{FT}`

### `Soil(Δz)` (Lines 89-100)
Convenience constructor from layer thicknesses:
```julia
soil = Soil([0.02, 0.04, 0.06, 0.10])  # 4 layers with varying thickness
```
- Automatically calls `soil_depth_init(Δz)` for grid geometry
- Initializes K and ψ via `cal_K!(soil)` and `cal_ψ!(soil)`

## Grid Initialization Functions

### `soil_depth_init(Δz)` (Lines 155-187)
Calculates staggered grid geometry from layer thicknesses.

**Returns**: `(z, z₋ₕ, z₊ₕ, Δz₊ₕ)`

**Grid conventions**:
- `z[i]`: Center depth of layer i [m] (negative, downward)
- `z₊ₕ[i]`: Lower interface depth at i+1/2 [m]
- `z₋ₕ[i]`: Upper interface depth at i-1/2 [m]
- `Δz₊ₕ[i]`: Distance between centers: `z[i] - z[i+1]` [m]

**Example** for `Δz = [0.02, 0.04, 0.06, 0.10]` m:
```julia
z₊ₕ = [-0.02, -0.06, -0.12, -0.22]     # Cumulative depth to interfaces
z[1] = -0.01                            # Center of first layer
z[2] = -0.04                            # Center of second layer
Δz₊ₕ[1] = z[1] - z[2] = 0.03           # Distance between centers 1 and 2
```

### `cal_Δz(z)` (Lines 190-202)
Inverse function: calculates layer thicknesses from center depths.

## Display Method (Lines 110-142)

Pretty-prints soil state showing:
1. Grid configuration (`N`, `ibeg`, `inds_obs`)
2. Vegetation (root distribution `f_root`)
3. Temperature fields (`Tsoil`, `Tsurf`)
4. Moisture fields (`K`, `ψ`, `θ`, `sink`, `Q`, boundary values)
5. Groundwater (`zwt`, `wa`)
6. Soil parameters

## Key Design Features

1. **Staggered Grid**: State variables at layer centers, fluxes at interfaces
2. **Flexible Boundaries**: `ibeg` allows inactive surface layers (e.g., for snow)
3. **Multi-Physics**: Coupled water flow, heat transfer, and vegetation
4. **Groundwater Coupling**: Dynamic water table with aquifer storage
5. **Type Parameterization**: `{FT, P}` for flexible precision and parameter types
6. **Offset Arrays**: `OffsetVector` for convenient 0-based indexing
7. **Dual Units**: Maintains both SI [m] and practical [cm] units

## Relationship to Richards Equation

The `Soil` struct is the data container that `RichardsEquation` operates on:
- `RichardsEquation` updates `dθ` using fields from `Soil`
- `cal_Q!` calculates fluxes stored in `soil.Q`
- `soil_Updateθ!` advances the solution forward in time
- Parameters in `soil.param` define hydraulic properties K(θ) and ψ(θ)
