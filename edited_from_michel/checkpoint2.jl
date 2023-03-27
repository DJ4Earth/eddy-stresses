# This is my attempt at a Julia equivalent of Milan Kloewer's shallow water model (code found here: https://github.com/milankl/swm)
# We'll solve the following equations on an Arakawa C-grid setup:
#       u_t = f v - g eta_x + A_h (u_xx + u_yy) + F_x
#       v_t = -f u - g eta_y + A_h (v_xx + v_yy)
#       eta_t = - (H u_x + H v_y)
# using a fully explicit solver.
# Will add more about how the code is structured at a future point....

# This script is almost identical to main_barotropic_gyre, except here I'm attempting to add in
# the Julia package Checkpointing, as I was running into memory issues applying Enzyme

using Plots, SparseArrays, Parameters, UnPack
using JLD2
using Enzyme, Checkpointing, Zygote

include("init_structs.jl")
include("init_params.jl")
include("build_grid.jl")
include("build_discrete_operators.jl")

@with_kw mutable struct SWM_pde_debug
    # To initialize struct just need to specify the following, the rest
    # of the entries will follow

    # Nu::Int
    # Nv::Int
    # NT::Int
    # Nq::Int

    u::Vector{Float64}
    v::Vector{Float64}
    eta::Vector{Float64}

    # since everything that matters to the derivative needs to live in a single structure,
    # this will also contain all of the placeholders for terms on the RHS of the system

    # Appear in advance
    umid::Vector{Float64}
    vmid::Vector{Float64}
    etamid::Vector{Float64}
    u0::Vector{Float64}
    v0::Vector{Float64}
    eta0::Vector{Float64}
    u1::Vector{Float64}
    v1::Vector{Float64}
    eta1::Vector{Float64}

    # Appear in comp_u_v_eta_t
    # h::Vector{Float64} = zeros(NT)          # height of water columns [meters]

    # h_u::Vector{Float64} = zeros(Nu)        # height of water columns interpolated to u-grid
    # h_v::Vector{Float64} = zeros(Nv)        # height of water columns interpolated to v-grid
    # h_q::Vector{Float64} = zeros(Nq)

    # U::Vector{Float64} = zeros(Nu)
    # V::Vector{Float64} = zeros(Nv)

    # IuT_u1::Vector{Float64} = zeros(NT)
    # IvT_v1::Vector{Float64} = zeros(NT)
    # kinetic::Vector{Float64} = zeros(NT)

    # kinetic_sq::Vector{Float64} = zeros(NT)

    # Gvx_v1::Vector{Float64} = zeros(Nq)
    # Guy_u1::Vector{Float64} = zeros(Nq)
    # q::Vector{Float64} = zeros(Nq)
    # p::Vector{Float64} = zeros(NT)

    # ITu_ksq::Vector{Float64} = zeros(Nu)
    # ITv_ksq::Vector{Float64} = zeros(Nv)
    # bfric_u::Vector{Float64} = zeros(Nu)
    # bfric_v::Vector{Float64} = zeros(Nv)

    # LLu_u1::Vector{Float64} = zeros(Nu)
    # LLv_v1::Vector{Float64} = zeros(Nv)
    # Mu::Vector{Float64} = zeros(Nu)
    # Mv::Vector{Float64} = zeros(Nv)

    # GTx_p::Vector{Float64} = zeros(Nu)
    u_t::Vector{Float64}
    # GTy_p::Vector{Float64} = zeros(Nv)
    v_t::Vector{Float64}
    # Gux_U::Vector{Float64} = zeros(NT)
    # Gvy_V::Vector{Float64} = zeros(NT)
    eta_t::Vector{Float64}
    # adv_u::Vector{Float64} = zeros(Nu)
    # adv_v::Vector{Float64} = zeros(Nv)

    # # Appear in comp_advection
    # AL1q::Vector{Float64} = zeros(NT)
    # AL2q::Vector{Float64} = zeros(NT)
    T::Int

end


function advance_debug1(states_rhs::SWM_pde_debug,
    grid::Grid,
    params::Params,
    interp::Interps,
    grad::Derivatives,
    advec::Advection
)

    @unpack u, v, eta, umid, vmid, etamid, u0, v0, eta0, u1, v1, eta1, u_t, v_t, eta_t = states_rhs

    # nx = grid.nx
    dt = params.dt

    # we now use RK4 as the timestepper, here I'm storing the coefficients needed for this
    rk_a = [1/6, 1/3, 1/3, 1/6]
    # rk_b = [1/2, 1/2, 1.]

    umid .= u
    vmid .= v
    etamid .= eta

    u0 .= u
    v0 .= v
    eta0 .= eta

    u1 .= u
    v1 .= v
    eta1 .= eta

    for j in 1:4

        # comp_u_v_eta_t(nx, states_rhs, params, interp, grad, advec)

        # if j < 4
        #     states_rhs.u1 .= states_rhs.umid .+ rk_b[j] .* dt .* states_rhs.u_t
        #     states_rhs.v1 .= states_rhs.vmid .+ rk_b[j] .* dt .* states_rhs.v_t
        #     states_rhs.eta1 .= states_rhs.etamid .+ rk_b[j] .* dt .* states_rhs.eta_t
        # end

        u0 .= u0 .+ rk_a[j] .* dt .* u_t
        # states_rhs.v0 .= states_rhs.v0 .+ rk_a[j] .* dt .* states_rhs.v_t
        # states_rhs.eta0 .= states_rhs.eta0 .+ rk_a[j] .* dt .* states_rhs.eta_t

    end

    @assert all(x -> x < 7.0, u0)
    @assert all(x -> x < 7.0, v0)
    @assert all(x -> x < 7.0, eta0)

    return nothing


end

function setup(; nx = 10, ny = 10, Lx = 3840e3, Ly = 3840e3)

    grid = build_grid(Lx, Ly, nx, ny)
    params = def_params(grid)

    # building discrete operators
    grad = build_derivs(grid)            # discrete gradient operators
    interp = build_interp(grid, grad)    # discrete interpolation operators (travels between grids)
    advec = build_advec(grid)

    Nu = grid.Nu
    Nv = grid.Nv
    NT = grid.NT

    u = zeros(Nu)
    v = zeros(Nv)
    eta = zeros(NT)
    umid = zeros(Nu)
    vmid = zeros(Nv)
    etamid = zeros(NT)
    u0 = zeros(Nu)
    v0 = zeros(Nv)
    eta0 = zeros(NT)
    u1 = zeros(Nu)
    v1 = zeros(Nv)
    eta1 = zeros(NT)

    u_t = zeros(Nu)
    v_t = zeros(Nv)
    eta_t = zeros(NT)

    states_rhs = SWM_pde_debug(u, v, eta,
        umid, vmid, etamid,
        u0, v0, eta0,
        u1, v1, eta1,
        u_t, v_t, eta_t,
        0
    )

    return grid, params, grad, interp, advec, states_rhs
end

function chkpt_maybe(
    chkpt_struct::SWM_pde_debug,
    chkpt_scheme::Scheme,
    T::Int,
    grid::Grid,
    params::Params,
    interp::Interps,
    grad::Derivatives,
    advec::Advection,
    i::Int
)
    @checkpoint_struct chkpt_scheme chkpt_struct for j in 1:chkpt_struct.T

        advance_debug1(chkpt_struct, grid, params, interp, grad, advec)

        copyto!(chkpt_struct.u, chkpt_struct.u0)
        copyto!(chkpt_struct.v, chkpt_struct.v0)
        copyto!(chkpt_struct.eta, chkpt_struct.eta0)

    end
    return chkpt_struct.u[i]
end


Trun = 10
snaps = 3
verbose = 0
revolve = Revolve{SWM_pde_debug}(Trun, snaps; verbose=verbose)

grid, gyre_params, grad_ops, interp_ops, advec_ops, chkpt_struct = setup()
chkpt_struct.T = Trun

function my_jacobian(chkpt_maybe,
    chkpt_struct,
    chkpt_scheme,
    Trun,
    grid,
    gyre_params,
    interp_ops,
    grad_ops,
    advec_ops,
)

    jacobian = zeros(grid.Nu, grid.Nu)
    for i in 1:grid.Nu
        Checkpointing.reset(chkpt_scheme)
        g = Zygote.gradient(chkpt_maybe,
            chkpt_struct,
            revolve,
            Trun,
            grid,
            gyre_params,
            interp_ops,
            grad_ops,
            advec_ops,
            i
        )
        jacobian[i,:] .= g[1].u
    end
    return jacobian
end

J = my_jacobian(
    chkpt_maybe,
    chkpt_struct,
    revolve,
    Trun,
    grid,
    gyre_params,
    interp_ops,
    grad_ops,
    advec_ops,
)
