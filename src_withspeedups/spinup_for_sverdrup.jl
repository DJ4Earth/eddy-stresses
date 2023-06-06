using Plots, SparseArrays, Parameters, UnPack
using JLD2, LinearAlgebra

include("init_structs.jl")
include("init_params.jl")
include("build_grid.jl")
include("build_discrete_operators.jl")
include("advance.jl")
include("cost_func.jl")
include("compute_time_deriv.jl")
include("temp.jl")

@time states = integrate(2, 128, 128)
# states_nx128_ny128_10yr_060523 = gyre_vector(states.u, states.v, states.eta)
# @save "states_nx128_ny128_10yr_060523.jld2" states_nx128_ny128_10yr_060523