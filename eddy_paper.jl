
"""
Placing the checkpointed integraton, loop, experiment functions, etc. here in one place
"""

include("../ShallowWaters.jl/src/ShallowWaters.jl")
using .ShallowWaters
using Enzyme
using Checkpointing, HDF5, Serialization
using NetCDF, JLD2, CairoMakie
using Lux

Enzyme.API.looseTypeAnalysis!(true)

using Parameters
using Optim
using LaTeXStrings

include("eddy_paper_integration.jl")
include("eddy_paper_optim.jl")