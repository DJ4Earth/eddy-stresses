using Plots, SparseArrays, Parameters
using JLD2
using Enzyme 
# using Enzyme_jll#main

include("init_structs.jl")
include("init_params.jl")
include("build_grid.jl")
include("build_discrete_operators.jl")
include("compute_time_deriv.jl")
include("advance_c_grid.jl")

nx = 10            # grid resolution in x-direction
ny = 10             # grid resolution in y-direction

Lx = 3840e3                     # E-W length of the domain [meters]
Ly = 3840e3                     # N-S length of the domain [meters]

grid_params = build_grid(Lx, Ly, nx, ny)
gyre_params = def_params(grid_params)

# building discrete operators 
grad_ops = build_derivs(grid_params)                # discrete gradient operators 
interp_ops = build_interp(grid_params, grad_ops)    # discrete interpolation operators (travels between grids)
advec_ops = build_advec(grid_params)
rhs_terms = RHS_terms(Nu = grid_params.Nu, Nv = grid_params.Nv, NT = grid_params.NT, Nq = grid_params.Nq)

Trun = days_to_seconds(100, gyre_params.dt)
    
uout = zeros(grid_params.Nu) 
vout = zeros(grid_params.Nv) 
etaout = zeros(grid_params.NT)

u = [zeros(grid_params.Nu)]
v = [zeros(grid_params.Nv)]
eta = [zeros(grid_params.NT)]

u_v_eta = gyre_vector(uout, vout, etaout)

@time for t in 1:Trun
    advance(u_v_eta, grid_params, rhs_terms, gyre_params, interp_ops, grad_ops, advec_ops)
    push!(u, copy(u_v_eta.u))
    push!(v, copy(u_v_eta.v))
    push!(eta, copy(u_v_eta.eta))
end

# statem = vec_to_mat(u[end], v[end], eta[end], grid_params)

ad_eta = zeros(grid_params.NT)
ad_eta[13] = 0.0;
ad_u_v_eta = gyre_vector(zeros(grid_params.Nu), zeros(grid_params.Nv), copy(ad_eta))

# runs, can check the derivative this returns to me, I don't think it's correct
u_v_eta_ = gyre_vector(copy(u[end]), copy(v[end]), copy(eta[end]))
# checking that the ad structure starts out as zeros
@show ad_u_v_eta
@show ad_u_v_eta.eta[13]
@show maximum(u_v_eta_.u), maximum(u_v_eta_.v), maximum(u_v_eta_.eta)
autodiff(advance, 
    Duplicated(u_v_eta_, ad_u_v_eta),
    Const(grid_params),
    Const(rhs_terms), 
    Const(gyre_params),
    Const(interp_ops),
    Const(grad_ops),
    Const(advec_ops)
)
@show ad_u_v_eta.eta[13]

# runs but creates divergences for some reason  
# for j = Trun:-1:1 
#u_v_eta_ = gyre_vector(u[j], v[j], eta[j])
# autodiff(advance, 
#     Duplicated(u_v_eta_, ad_u_v_eta),
#     grid_params,
#     rhs_terms, 
#     gyre_params, 
#     interp_ops, 
#     grad_ops,
#     advec_ops, 
# )
# end

# ad_grid = deepcopy(grid_params)
# ad_grad = deepcopy(grad_ops)
# ad_interp = deepcopy(interp_ops)
# ad_advec = deepcopy(advec_ops)
# ad_params = deepcopy(gyre_params)
# ad_rhs = deepcopy(rhs_terms)

# doesn't run 
# autodiff(advance, 
#     Duplicated(u_v_eta_, ad_u_v_eta),
#     Duplicated(grid_params, ad_grid),
#     Duplicated(rhs_terms, ad_rhs),
#     Duplicated(gyre_params, ad_params),
#     Duplicated(interp_ops, ad_interp),
#     Duplicated(grad_ops, ad_grad),
#     Duplicated(advec_ops, ad_advec)
# )

