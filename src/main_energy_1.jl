# making an effort at the whole experiment here, done with Enzyme + checkpointing 
# once this is running and working will structure a main file with it. 

# This one will try to use Enzyme + checkpointing for an energy sensitivity
# similar to the Burgers equation

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
function setup(; days = 10, nx = 20, ny = 20, Lx = 3840e3, Ly = 3840e3)

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

    return grid, params, grad, interp, advec, states_rhs
end


function chkpt_maybe(
    chkpt_struct::SWM_pde,
    chkpt_scheme::Scheme,
    grid::Grid,
    params::Params,
    interp::Interps,
    grad::Derivatives,
    advec::Advection
)

    @checkpoint_struct chkpt_scheme chkpt_struct for t in 1:chkpt_struct.T

        advance(chkpt_struct, grid, params, interp, grad, advec)

        copyto!(chkpt_struct.u, chkpt_struct.u0)
        copyto!(chkpt_struct.v, chkpt_struct.v0)
        copyto!(chkpt_struct.eta, chkpt_struct.eta0)

    end

    energy = sum(chkpt_struct.u.^2 .+ chkpt_struct.v.^2) / (grid.nx * grid.ny)

    return energy

end

grid, gyre_params, grad_ops, interp_ops, advec_ops, chkpt_struct = setup()

snaps = 50
verbose = 0
revolve = Revolve{SWM_pde}(chkpt_struct.T, snaps; verbose=verbose)

denergy = Zygote.gradient(chkpt_maybe, 
    chkpt_struct, 
    revolve, 
    grid, 
    gyre_params, 
    interp_ops, 
    grad_ops, 
    advec_ops
)


# for checking that the forward integration still works 

# for t = 1:chkpt_struct.T
#     advance(chkpt_struct, grid, gyre_params, interp_ops, grad_ops, advec_ops)
#     copyto!(chkpt_struct.u, chkpt_struct.u0)
#     copyto!(chkpt_struct.v, chkpt_struct.v0)
#     copyto!(chkpt_struct.eta, chkpt_struct.eta0)
# end