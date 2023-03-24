# This is my attempt at a Julia equivalent of Milan Kloewer's shallow water model (code found here: https://github.com/milankl/swm)
# We'll solve the following equations on an Arakawa C-grid setup:
#       u_t = f v - g eta_x + A_h (u_xx + u_yy) + F_x 
#       v_t = -f u - g eta_y + A_h (v_xx + v_yy)
#       eta_t = - (H u_x + H v_y)
# using a fully explicit solver. 
# Will add more about how the code is structured at a future point.... 

# This script is almost identical to main_barotropic_gyre, except here I'm attempting to add in 
# the Julia package Checkpointing, as I was running into memory issues applying Enzyme

using Plots, SparseArrays, Parameters
using JLD2
using Enzyme, Checkpointing, Zygote

include("init_structs.jl")
include("init_params.jl")
include("build_grid.jl")
include("build_discrete_operators.jl")
include("compute_time_deriv.jl")

function advance_debug(states_rhs::SWM_pde, 
    grid::Grid, 
    params::Params, 
    interp::Interps, 
    grad::Derivatives, 
    advec::Advection
) 

    nx = grid.nx 
    dt = params.dt

    # we now use RK4 as the timestepper, here I'm storing the coefficients needed for this 
    rk_a = [1/6, 1/3, 1/3, 1/6]
    rk_b = [1/2, 1/2, 1.]

    states_rhs.umid .= states_rhs.u
    states_rhs.vmid .= states_rhs.v
    states_rhs.etamid .= states_rhs.eta

    states_rhs.u0 .= states_rhs.u
    states_rhs.v0 .= states_rhs.v
    states_rhs.eta0 .= states_rhs.eta

    states_rhs.u1 .= states_rhs.u
    states_rhs.v1 .= states_rhs.v
    states_rhs.eta1 .= states_rhs.eta

    # for j in 1:4

    #     # comp_u_v_eta_t(nx, states_rhs, params, interp, grad, advec)

    #     if j < 4
    #         states_rhs.u1 .= states_rhs.umid .+ rk_b[j] .* dt .* states_rhs.u_t
    #         states_rhs.v1 .= states_rhs.vmid .+ rk_b[j] .* dt .* states_rhs.v_t
    #         states_rhs.eta1 .= states_rhs.etamid .+ rk_b[j] .* dt .* states_rhs.eta_t
    #     end

    #     states_rhs.u0 .= states_rhs.u0 .+ rk_a[j] .* dt .* states_rhs.u_t
    #     states_rhs.v0 .= states_rhs.v0 .+ rk_a[j] .* dt .* states_rhs.v_t 
    #     states_rhs.eta0 .= states_rhs.eta0 .+ rk_a[j] .* dt .* states_rhs.eta_t 

    # end

    @assert all(x -> x < 7.0, states_rhs.u0)
    @assert all(x -> x < 7.0, states_rhs.v0)
    @assert all(x -> x < 7.0, states_rhs.eta0)

    return nothing 

end 

function setup(; nx = 10, ny = 10, Lx = 3840e3, Ly = 3840e3)

    grid = build_grid(Lx, Ly, nx, ny)
    params = def_params(grid)

    # building discrete operators 
    grad = build_derivs(grid)                # discrete gradient operators 
    interp = build_interp(grid, grad)    # discrete interpolation operators (travels between grids)
    advec = build_advec(grid)
    #rhs = RHS_terms(Nu = grid.Nu, Nv = grid.Nv, NT = grid.NT, Nq = grid.Nq)

    return grid, params, grad, interp, advec

end

function chkpt_maybe(chkpt_struct::SWM_pde, T::Int)
    grid, gyre_params, grad_ops, interp_ops, advec_ops = setup(nx=20, ny=20)
    @checkpoint_struct revolve chkpt_struct for j in 1:T
        advance_debug(chkpt_struct, grid, gyre_params, interp_ops, grad_ops, advec_ops)
        copyto!(chkpt_struct.u, chkpt_struct.u0)
        copyto!(chkpt_struct.v, chkpt_struct.v0)
        copyto!(chkpt_struct.eta, chkpt_struct.eta0)
    end
    return chkpt_struct.u
end

nx = 20
ny = 20 

deriv_struct = SWM_pde(Nu = (nx - 1) * ny, 
    Nv = (ny - 1) * nx, 
    NT = nx * ny, 
    Nq = (nx + 1) * (ny + 1)
)

Trun = 10
snaps = 3
verbose = 0
revolve = Revolve{SWM_pde}(Trun, snaps; verbose=verbose)

g = Zygote.jacobian(chkpt_maybe, deriv_struct, Trun)

# grid, gyre_params, grad_ops, interp_ops, advec_ops = setup(nx=20,ny=20)
# for j in 1:Trun
#     advance2(deriv_struct, grid, gyre_params, interp_ops, grad_ops, advec_ops)
# end

# state_m = vec_to_mat(deriv_struct.u, deriv_struct.v, deriv_struct.eta, grid)
