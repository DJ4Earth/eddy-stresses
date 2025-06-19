using Enzyme
using Checkpointing, HDF5, Serialization
using NetCDF, JLD2, CairoMakie
using Lux, Random
using Reactant
using DSP, FFTW, AbstractFFTs
using ChainRules
Enzyme.EnzymeRules.inactive(::typeof(plan_fft), args...; kwargs...) = true
Enzyme.EnzymeRules.inactive(::typeof(plan_rfft), args...; kwargs...) = true
Enzyme.@import_rrule(typeof(*), AbstractFFTs.Plan, AbstractArray)
Enzyme.@import_rrule(typeof(*), AbstractFFTs.ScaledPlan, AbstractArray)

using Parameters
using Optim
using LaTeXStrings

if !Base.isdefined(@__MODULE__, :ShallowWaters)
    include("../ShallowWaters.jl/src/ShallowWaters.jl")
    using .ShallowWaters
end

include("timeavgenergy_10dayintegration_nnrun.jl")
include("kespectrum_10dayintegration_nnrun.jl")
include("kespectrum_weightsandbias_10dayintegration_nnrun.jl")
include("multipletimeseries_kespectrum_10dayintegration_weightandbias.jl")
