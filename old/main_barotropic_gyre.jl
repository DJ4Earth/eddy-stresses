# This is my attempt at a Julia equivalent of Milan Kloewer's shallow water model (code found here: https://github.com/milankl/swm)
# We'll solve the following equations on an Arakawa C-grid setup:
#       u_t = f v - g eta_x + A_h (u_xx + u_yy) + F_x 
#       v_t = -f u - g eta_y + A_h (v_xx + v_yy)
#       eta_t = - (H u_x + H v_y)
# using a fully explicit solver. 
# Will add more about how the code is structured at a future point.... 


using Plots, SparseArrays, Parameters
using JLD2, LinearAlgebra
using Enzyme 

using InteractiveUtils

# import Profile.Allocs: @profile
import Profile: @profile 
import PProf

include("init_structs.jl")
include("init_params.jl")
include("build_grid.jl")
include("build_discrete_operators.jl")
include("compute_time_deriv.jl")
include("advance.jl")


# function ad_calc(grid, rhs, params, interp, grad, advec)

#     ad_rhs = RHS_terms(grid.Nu, grid.Nv, grid.NT, grid.Nq)

#     T = days_to_seconds(2, params.dt)

#     ad_eta = zeros(grid.NT)
#     ad_eta[13] = 1.0;
#     ad_u_v_eta = gyre_vector(zeros(grid.Nu), zeros(grid.Nv), ad_eta)

#     u_v_eta = gyre_vector(zeros(grid.Nu), zeros(grid.Nv), ad_eta)
#     main(T, u_v_eta, grid, rhs, params, interp, grad, advec)

#     # computing and storing all the states 
#     # u, v, eta = integrate(1, 5, 5)

#     for j in T+1:-1:1 

#         u_v_eta = gyre_vector(u[end], v[end], eta[end])

#         autodiff(main, 
#             T,
#             Duplicated(u_v_eta, ad_u_v_eta),
#             DuplicatedNoNeed(grid, ad_grid),
#             DuplicatedNoNeed(rhs, ad_rhs),
#             DuplicatedNoNeed(params, ad_params),
#             DuplicatedNoNeed(interp, ad_interp),
#             DuplicatedNoNeed(grad, ad_grad),
#             DuplicatedNoNeed(advec, ad_advec)
#         )

#     end

#     return u, v, eta, ad_u_v_eta

# end

nx = 5             # grid resolution in x-direction
ny = 5             # grid resolution in y-direction

Lx = 3840e3                     # E-W length of the domain [meters]
Ly = 3840e3                     # N-S length of the domain [meters]

grid_params = build_grid(Lx, Ly, nx, ny)
gyre_params = def_params(grid_params)

# building discrete operators 
grad_ops = build_derivs(grid_params)                # discrete gradient operators 
interp_ops = build_interp(grid_params, grad_ops)    # discrete interpolation operators (travels between grids)
advec_ops = build_advec(grid_params)
rhs_terms = RHS_terms(Nu = grid_params.Nu, Nv = grid_params.Nv, NT = grid_params.NT, Nq = grid_params.Nq)

# u, v, eta, ad = ad_calc(grid_params, rhs_terms, gyre_params, interp_ops, grad_ops, advec_ops)

# u_v_eta = gyre_vector(zeros(grid_params.Nu), zeros(grid_params.Nv), zeros(grid_params.NT))
