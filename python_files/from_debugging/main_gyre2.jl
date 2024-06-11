# This is my attempt at a Julia equivalent of Milan Kloewer's shallow water model (code found here: https://github.com/milankl/swm)
# We'll solve the following equations on an Arakawa C-grid setup:
#       u_t = f v - g eta_x + A_h (u_xx + u_yy) + F_x 
#       v_t = -f u - g eta_y + A_h (v_xx + v_yy)
#       eta_t = - (H u_x + H v_y)
# using a fully explicit solver. 
# Will add more about how the code is structured at a future point....

using Enzyme, Plots, SparseArrays, Parameters
using InteractiveUtils

import Profile.Allocs: @profile
import PProf

include("init_structs.jl")
include("init_params.jl")
include("build_grid.jl")
include("build_discrete_operators.jl")
include("allocate_rhs.jl")
include("setup.jl")
include("advance_c_grid.jl")
include("compute_time_deriv.jl")
include("checking2.jl")
include("temp.jl")

# Zanna / Bolton setup 
# Lx = 3840e3                    # E-W length of the domain [meters]
# Ly = 3840e3                    # N-S length of the domain [meters]
# nx = 100                        # number of cells in the x-direction
# ny = 100                        # number of cells in the y-direction

# MIT GCM setup 
# Lx = 1200e3 
# Ly = 1200e3
# nx = 62 
# ny = 62 

# grid_params, gyre_params, grad_ops, interp_ops, advec_ops, rhs_terms = setup(Lx, Ly, nx, ny)

function main(Tspinup, Trun)

    Lx = 3840e3                    # E-W length of the domain [meters]
    Ly = 3840e3                    # N-S length of the domain [meters]
    nx = 20                        # number of cells in the x-direction
    ny = 20                        # number of cells in the y-direction

    grid_params, gyre_params, grad_ops, interp_ops, advec_ops, rhs_terms = setup(Lx, Ly, nx, ny)

    # starting from rest ---> all initial conditions are zero 

    # how long to spinup the model for 
    Tspinup_days = Tspinup # [days] 

    # how long to run the model for after spinup
    Trun_days = Trun     # [days] 

    Tspinup, Trun = days_to_seconds(Tspinup_days, Trun_days, gyre_params.dt)

    uout = zeros(grid_params.Nu) 
    vout = zeros(grid_params.Nv) 
    etaout = zeros(grid_params.NT)

    # uout = (collect(LinRange(0, grid_params.Nu - 1, grid_params.Nu)) ./ 1000).^2
    # vout = (collect(LinRange(0, grid_params.Nv - 1, grid_params.Nv)) ./ 1000).^2
    # etaout = (collect(LinRange(0, grid_params.NT - 1, grid_params.NT)) ./ 1000).^2

    u_v_eta = gyre_vector(
        copy(uout), 
        copy(vout), 
        copy(etaout)
    )

    # @show @code_lowered advance_check2(u_v_eta, rhs_terms, gyre_params, grid_params, interp_ops, grad_ops, advec_ops) 
    # @show @code_typed debuginfo=:source advance_check2(u_v_eta, rhs_terms, gyre_params, grid_params, interp_ops, grad_ops, advec_ops) 
    # @code_llvm advance_check2(u_v_eta, rhs_terms, gyre_params, grid_params, interp_ops, grad_ops, advec_ops) 
    # exit()

    # @profile sample_rate=0.1 for t in 1:Tspinup
    #     advance_check2(u_v_eta, rhs_terms, gyre_params, grid_params, interp_ops, grad_ops, advec_ops) 
    # end

    @time for t in 1:Tspinup
        advance_check2(u_v_eta, rhs_terms, gyre_params, grid_params, interp_ops, grad_ops, advec_ops) 
    end

    # PProf.Allocs.pprof(;web=false, out="allocs.pb.gz")

    u_v_eta_mat = vec_to_mat(u_v_eta.u, u_v_eta.v, u_v_eta.eta, grid_params)

    return u_v_eta, u_v_eta_mat

end

# heatmap(u_v_eta_mat.eta)

# u_v_eta_start = deepcopy(u_v_eta)

# @time for t in 1:Tspinup
#     advance(u_v_eta, gyre_params, interp_ops, grad_ops, advec_ops) 
# end

# states, states_mat = main(20, 20);