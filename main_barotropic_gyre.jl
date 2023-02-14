# Mimicking the structure of Milan Kloewer's shallow water model (code found here: https://github.com/milankl/swm)
# We'll solve the equations on an Arakawa C-grid

using Enzyme, Plots

include("init_structs.jl")
include("build_structs.jl")
include("build_grid.jl")
include("build_discrete_operators.jl")
include("advance_c_grid.jl")
include("compute_time_deriv.jl")

# setting up the grid parameters first

# Nt will be the number of time steps we take to spin-up the model 
Tspinup = 1000   # [seconds]

# T will be the number of steps to take after spinning up the model 
Tsteps = 5000    # [seconds]

Lx = 3840000                    # E-W length of the domain [meters]
Ly = 3840000                    # N-S length of the domain [meters]
nx = 5                         # number of cells in the x-direction
ny = 5                         # number of cells in the y-direction

# based on above values this returns more parameters related to the four grids 
dx, dy, NT, Nu, Nv, Nq, x, y, xu, yu, xv, yv, xq, yq = build_grid(Lx, Ly, nx, ny)

g = 9.81    # gravity [meters^2 / second]
H = 500.0   # depth of the box [meters]

dt = (0.9 * min(dx, dy)) / (sqrt(g * H))   # CFL condition for dt [seconds]

grid_params = Grid(
    Lx, 
    Ly,
    Tsteps,
    nx, 
    ny,
    NT,
    Nu, 
    Nv,
    Nq,
    dx, 
    dy
)

# initializing parameter values 

# first defining the coriolis and wind stress fields 
omega = 2 * pi / (24 * 3600)
R = 6.371e6
f0 = 2 * omega * sin(30 * pi / 180)
beta = (2 * omega / R) * cos(30 * pi / 180)

A_h = (128*540/(min(grid_params.nx, grid_params.ny))) * max(grid_params.dx, grid_params.dy)^2       # viscosity coefficient [meters^2 / second]
rho_c = 1000.0  

drag = 1e-5

# come back to this, want to understand 
Yq_shift = yq .- grid_params.Ly/2
Yq = vec([k for k in Yq_shift, j in 1:grid_params.nx+1]')
f(y) = f0 + beta * y
coriolis = vec(f.(Yq)')

Yu = vec([k for k in yu, j in 1:length(xu)]')
wind_stress = (0.12 .* (cos.(2 * pi .* ((Yu .- Ly/2)./Ly) ) + 2 .* sin.(2 * pi .* ((Yu .- Ly/2)./Ly)  ) ) ) ./rho_c

gyre_params = Parameters(
    Tspinup,
    Tsteps,
    dt,
    g, 
    f0, 
    beta, 
    H, 
    A_h,
    rho_c, 
    eps_ab,
    drag,  
    wind_stress, 
    coriolis
)

# building discrete operators 
GTx, GTy, Gux, Guy, Gvx, Gvy, Gqy, Gqx, Lu, Lv, LT, Lq = build_derivs(grid_params)

grad_ops = Derivatives(
    GTx, 
    GTy, 
    Gux, 
    Guy, 
    Gvx, 
    Gvy, 
    Gqy, 
    Gqx, 
    Lu, 
    Lv, 
    LT, 
    Lq
)

Ivu, Iuv, IqT, IuT, IvT, ITu, ITv, Iqu, Iuq, Iqv, Ivq, ITq = build_interp(grid_params, grad_ops)

interp_ops = Interps(Ivu, 
    Iuv, 
    IqT, 
    IuT, 
    IvT, 
    ITu, 
    ITv, 
    Iqu, 
    Iuq, 
    Iqv, 
    Ivq, 
    ITq
)

AL1, AL2, index_av, index_bv, index_cv, index_dv, ALeur, ALeul, Seul, Seur, Sau, Sbu, Scu, Sdu, ALpvu, ALpvd, Spvu, Spvd, Sav, Sbv, Scv, Sdv = build_advec(grid_params)

advec_ops = Advection(AL1,
    AL2,
    index_av, 
    index_bv, 
    index_cv, 
    index_dv, 
    ALeur, 
    ALeul, 
    Seul, 
    Seur, 
    Sau, 
    Sbu, 
    Scu, 
    Sdu, 
    ALpvu, 
    ALpvd, 
    Spvu, 
    Spvd, 
    Sav, 
    Sbv, 
    Scv, 
    Sdv
)

# starting from rest ---> all initial conditions are zero 

u = zeros(Nu)
v = zeros(Nv)
eta = zeros(NT)

u_v_eta = gyre_vector(u, v, eta)

for t = 1:gyre_params.Tspinup
    advance(u_v_eta, gyre_params, interp_ops, grad_ops, advec_ops) 
end

# heatmap(u_v_eta.η)
