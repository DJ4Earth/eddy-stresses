# making an effort at the whole experiment here, done with Enzyme + checkpointing 
# once this is running and working will structure a main file with it. 

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

# This function will setup the structures needed to integrate the model. Comes with default values, but
# these can be specified if desired. 
function setup(; days = 100, nx = 10, ny = 10, Lx = 3840e3, Ly = 3840e3)

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

    data_steps = 2:10:T
    scaling = 4

    data, M = create_data(days, nx, ny, data_steps, scaling = scaling)

    states_rhs = SWM_pde(Nu = Nu, 
        Nv = Nv,
        NT = NT, 
        Nq = Nq, 
        T = T,
        scaling = scaling
    )

    return grid, params, grad, interp, advec, states_rhs, data, data_steps, M
end


function chkpt_maybe(
    chkpt_struct::SWM_pde,
    chkpt_scheme::Scheme,
    data::Matrix{Float64},
    data_steps,
    M, 
    grid::Grid,
    params::Params,
    interp::Interps,
    grad::Derivatives,
    advec::Advection
)

    j = 1
    @checkpoint_struct chkpt_scheme chkpt_struct for t = 1:chkpt_struct.T

        advance(chkpt_struct, grid, params, interp, grad, advec)

        if t in data_steps 
            chkpt_struct.J += cost_func(M, data[:, j], 
                chkpt_struct.u0,
                chkpt_struct.v0,
                chkpt_struct.eta0
            )
            j += 1
        end

        copyto!(chkpt_struct.u, chkpt_struct.u0)
        copyto!(chkpt_struct.v, chkpt_struct.v0)
        copyto!(chkpt_struct.eta, chkpt_struct.eta0)

    end

end

grid, gyre_params, grad_ops, interp_ops, advec_ops, chkpt_struct, data, data_steps, M = setup();

snaps = 3
verbose = 0
revolve = Revolve{SWM_pde}(chkpt_struct.T, snaps; verbose=verbose)

dnu = Zygote.gradient(chkpt_maybe, 
    chkpt_struct, 
    revolve, 
    data,
    data_steps,
    M, 
    grid, 
    gyre_params, 
    interp_ops, 
    grad_ops, 
    advec_ops
)

# # for checking that the forward integration still works 

# for t = 1:chkpt_struct.T
#     advance(chkpt_struct, grid, gyre_params, interp_ops, grad_ops, advec_ops)
#     copyto!(chkpt_struct.u, chkpt_struct.u0)
#     copyto!(chkpt_struct.v, chkpt_struct.v0)
#     copyto!(chkpt_struct.eta, chkpt_struct.eta0)
# end