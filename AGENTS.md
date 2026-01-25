# SoilDiffEqs.jl Agent Guidelines

This document provides instructions for AI agents working on the SoilDiffEqs.jl repository.

## 代码编写原则

- linux极简主义；同时排版符合代码规范，排版美观

## 1. Project Overview
SoilDiffEqs.jl implements solvers for soil water movement (Richards equation) and heat transfer.
- **Core Physics**: Richards equation, Darcy's law, Van Genuchten / Campbell retention curves.
- **Key Dependencies**: `DifferentialEquations.jl`, `OrdinaryDiffEqTsit5`.
- **Language**: Julia (v1.10+).

## 2. Build and Test Commands

### Running Tests
This project uses the standard Julia `Test` suite.

*   **Run All Tests**:
    ```bash
    julia --project -e "using Pkg; Pkg.test()"
    ```

*   **Run a Single Test File**:
    To run a specific test file (e.g., `test/test-soil_moisture.jl`), you must activate the project environment first:
    ```bash
    julia --project -e "using Pkg; Pkg.activate(\".\"); using SoilDifferentialEquations; include(\"test/test-soil_moisture.jl\")"
    ```

### Development Environment
*   Ensure you are using the project environment: `julia --project`
*   No specific linter/formatter config is present, but follow existing indentation (2 spaces).

## 3. Code Style Guidelines

### Naming Conventions
*   **Functions**: `snake_case` (e.g., `cal_Q!`, `soil_Updateθ!`). Use `!` suffix for mutating functions.
*   **Types/Structs**: `CamelCase` (e.g., `Soil`, `SoilParam`, `VanGenuchten`).
*   **Variables**:
    *   Math symbols are strictly encouraged: `θ` (theta), `ψ` (psi), `Δz` (delta z), `K₊ₕ` (K at half level).
    *   Mix `snake_case` with math symbols: `θ_res`, `θ_sat`.
    *   `ibeg` used for the starting active layer index.

### Type System
*   **Parametric Types**: heavily used for performance and precision control.
    ```julia
    function cal_Q!(soil::Soil{T}; ...) where {T<:Real}
    ```
*   **Explicit Typing**: Annotate function arguments and return types where helpful for dispatch or clarity.

### Formatting
*   **Indentation**: 2 spaces.
*   **Performance**:
    *   Use `@inbounds` for hot loops.
    *   Avoid allocations in inner loops (pre-allocate arrays in `Soil` struct).
    *   Use `(; var1, var2) = struct` for unpacking.

### Documentation
*   Use Markdown docstrings above functions.
*   Comments can be in **Chinese** or **English**. Existing code has mixed comments.

## 4. Specific Patterns & Rules

### Copilot / AI Instructions (from `.github/copilot-instructions.md`)
1.  **Response Language**: Use **Chinese** for conversation/explanations.
2.  **Code Output**: Provide **ONLY** the modified parts of the code. Do not dump the entire file unless necessary.
3.  **Typset**:
    *   Use Typst format for formulas in documentation if requested.
    *   Greek letters in Typst do not need a slash (e.g., `theta`, not `\theta`).

### Physics Implementation
*   **Grid**: Staggered grid. State variables (`θ`, `ψ`) at centers `z`, fluxes (`Q`, `K`) at interfaces `z₊ₕ`.
*   **Units**:
    *   Water content `θ`: [m³ m⁻³]
    *   Flux `Q`: [cm h⁻¹]
    *   Potential `ψ`: [cm] (negative)
    *   Time `dt`: converted to hours in calculations (input often in seconds).

### Common Anti-Patterns to Avoid
*   Do not replace Unicode math symbols with ASCII names (e.g., keep `θ`, don't change to `theta` in code).
*   Do not remove `@inbounds` without verification.
*   Do not change the `Soil` struct layout without understanding the memory implications.

## 5. File Structure
*   `src/Soil.jl`: Main data structure.
*   `src/SoilMoisture/Equation_Richards.jl`: Core physics (flux calculation, update logic).
*   `test/`: Test files. `runtests.jl` aggregates them.
