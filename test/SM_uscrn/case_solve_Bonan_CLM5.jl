using SoilDifferentialEquations, Test

"""
Quick demo: Soil moisture simulation (Bonan) on CLM5 layers.
Config-driven version (see test/SM_uscrn/config_Bonan_CLM5.yaml).
"""

include(joinpath(@__DIR__, "run_Bonan_CLM5_config.jl"))
