
"""
Placing the checkpointed integraton, loop, experiment functions, etc. here in one place
"""

using Enzyme
Enzyme.Compiler.VERBOSE_ERRORS[] = true
using Checkpointing, HDF5, Serialization
using NetCDF, JLD2, CairoMakie
using Lux, Random

Enzyme.API.looseTypeAnalysis!(true)

using Parameters
using Optim
using LaTeXStrings

include("../ShallowWaters.jl/src/ShallowWaters.jl")
using .ShallowWaters

include("eddy_paper_integration.jl")
include("eddy_paper_optim.jl")
include("eddy_paper_experiment_functions.jl")

include("eddy_paper_run_experiments.jl")