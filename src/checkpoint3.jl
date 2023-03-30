# ****this script is just for debugging, I don't like editing main code directly until I know what I'm doing works
# This script is almost identical to main_barotropic_gyre, except here I'm attempting to add in
# the Julia package Checkpointing, as I was running into memory issues applying Enzyme

using Plots, SparseArrays, Parameters, UnPack
using JLD2, LinearAlgebra
using Enzyme, Checkpointing, Zygote

include("init_structs.jl")
include("init_params.jl")
include("build_grid.jl")
include("build_discrete_operators.jl")
include("advance.jl")
include("cost_func.jl")
include("compute_time_deriv.jl")


function setup(; days = 5, nx = 10, ny = 10, Lx = 3840e3, Ly = 3840e3)

    grid = build_grid(Lx, Ly, nx, ny)
    params = def_params(grid)

    # building discrete operators
    grad = build_derivs(grid)            # discrete gradient operators
    interp = build_interp(grid, grad)    # discrete interpolation operators (travels between grids)
    advec = build_advec(grid)

    Nu = grid.Nu
    Nv = grid.Nv
    NT = grid.NT
    Nq = grid.Nq 

    T = days_to_seconds(days, params.dt)

    states_rhs = SWM_pde(Nu = Nu, 
        Nv = Nv,
        NT = NT, 
        Nq = Nq, 
        T = T
    )

    data = create_data(days, nx, ny)

    return grid, params, grad, interp, advec, states_rhs, data
end



function chkpt_maybe(
    chkpt_struct::SWM_pde_debug,
    chkpt_scheme::Scheme,
    grid::Grid,
    params::Params,
    interp::Interps,
    grad::Derivatives,
    advec::Advection
)
    @checkpoint_struct chkpt_scheme chkpt_struct for j in 1:chkpt_struct.T

        advance_debug1(chkpt_struct, grid, params, interp, grad, advec)

        costfunction(state, data(t))
        copyto!(chkpt_struct.u, chkpt_struct.u0)
        copyto!(chkpt_struct.v, chkpt_struct.v0)
        copyto!(chkpt_struct.eta, chkpt_struct.eta0)

    end
    return 
end

grid, gyre_params, grad_ops, interp_ops, advec_ops, chkpt_struct = setup()

snaps = 3
verbose = 0
revolve = Revolve{SWM_pde_debug}(chkpt_struct.T, snaps; verbose=verbose)


# for checking that the forward integration still works 

for t = 1:chkpt_struct.T
    advance_debug1(chkpt_struct, grid, gyre_params, interp_ops, grad_ops, advec_ops)
    copyto!(chkpt_struct.u, chkpt_struct.u0)
    copyto!(chkpt_struct.v, chkpt_struct.v0)
    copyto!(chkpt_struct.eta, chkpt_struct.eta0)
end