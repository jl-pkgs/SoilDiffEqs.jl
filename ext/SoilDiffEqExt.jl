module SoilDiffEqExt

import OrdinaryDiffEqTsit5: ODEProblem, solve
import SoilDifferentialEquations:_ODEProblem, _solve

_ODEProblem(f, u0, tspan, soil) = ODEProblem(f, u0, tspan, soil)
_solve(arg...; kw...) = solve(arg...; kw...)

end
