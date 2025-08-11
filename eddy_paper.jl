using Enzyme
using Checkpointing, ImageFiltering
using NetCDF, JLD2, CairoMakie
using Lux, Random
using DSP, FFTW, AbstractFFTs
using ChainRules, LinearAlgebra
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

include("save_states.jl")
include("timeavgenergy_10dayintegration_nnrun.jl")
include("kespectrum_weightsandbias_10dayintegration_nnrun.jl")
include("kespectrum_percentdiff_weightsandbias_10dayintegration_nnrun.jl")
include("multipletimeseries_kespectrum_10dayintegration_weightandbias.jl")
include("multipletimeseries_states_10dayintegration_weightandbias.jl")
include("multipletimeseries_states_witheta_10dayintegration_weightsandbias.jl")
include("kespectrum_percentdiff_mycg_weightsandbias_nnrun.jl")
include("multipletimeseries_kespecpd_10dayintegration_weightandbias.jl")
