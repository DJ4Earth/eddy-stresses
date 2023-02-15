# Mimicking the structure of Milan Kloewer's shallow water model (code found here: https://github.com/milankl/swm)
# We'll solve the equations on an Arakawa C-grid

using Enzyme, Plots, SparseArrays

include("init_structs.jl")
include("init_params.jl")
include("build_grid.jl")
include("build_discrete_operators.jl")
include("advance_c_grid.jl")
include("compute_time_deriv.jl")

# setting up the grid parameters first

Lx = 3840000                    # E-W length of the domain [meters]
Ly = 3840000                    # N-S length of the domain [meters]
nx = 50                         # number of cells in the x-direction
ny = 50                         # number of cells in the y-direction

# based on above values this returns more parameters related to the four grids 
grid_params = build_grid(Lx, Ly, nx, ny)

gyre_params = def_params(grid_params)

# building discrete operators 
grad_ops = build_derivs(grid_params)                # discrete gradient operators 
interp_ops = build_interp(grid_params, grad_ops)    # discrete interpolation operators (travels between grids)
advec_ops = build_advec(grid_params)

# starting from rest ---> all initial conditions are zero 

# how long to spinup the model for 
Tspinup_days = 1*365 # [days] 

# how long to run the model for after spinup
Trun_days = 1 * 365     # [days] 

Tspinup, Trun = days_to_seconds(Tspinup_days, Trun_days, gyre_params.dt)

u = zeros(grid_params.Nu)
v = zeros(grid_params.Nv)
eta = zeros(grid_params.NT)

u_v_eta = gyre_vector(u, v, eta)


@time for t = 1:Tspinup
    advance(u_v_eta, gyre_params, interp_ops, grad_ops, advec_ops) 
end


# heatmap(u_v_eta.η)
