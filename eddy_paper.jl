using Enzyme
using Checkpointing, ImageFiltering
using NetCDF, JLD2, CairoMakie
using Lux, Random
using DSP, FFTW, AbstractFFTs
using ChainRules, LinearAlgebra
using NLPModels, MadNLP
using KernelDensity
Enzyme.EnzymeRules.inactive(::typeof(plan_fft), args...; kwargs...) = true
Enzyme.EnzymeRules.inactive(::typeof(plan_rfft), args...; kwargs...) = true
Enzyme.@import_rrule(typeof(*), AbstractFFTs.Plan, AbstractArray)
Enzyme.@import_rrule(typeof(*), AbstractFFTs.ScaledPlan, AbstractArray)

using Parameters
using LaTeXStrings

Random.seed!(8)

if !Base.isdefined(@__MODULE__, :ShallowWaters)
    include("../ShallowWaters.jl/src/ShallowWaters.jl")
    using .ShallowWaters
end

include("helper_functions.jl")