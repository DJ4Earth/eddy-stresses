# now going to try a sample adjoint calculation. just going to look at the sensitivity of the final 
# displacement in center of grid w.r.t. initial displacement. 

using Enzyme, Plots

include("build_structs.jl")
include("advance.jl")

g = 9.81 
f0 = 1e-4
beta = 1e-11
H = 5000.0
A_h = 400
rho_c = 1000.0
tau0 = 0.1

M = 10000

eps_ab = 0.1 

τ(y) = -tau0 + sin( pi * y / (1200e3) )
f(y) = f0 + beta * y 

Nx = 62
Ny = 62 

dx = 1e3 * 20.0
dy = 1e3 * 20.0

dt = 0.5 * dx / (2 * sqrt(g * H))

x = 0:dx:dx*(Nx - 1)
y = 0:dy:dy*(Ny - 1)

X = [j for j in x, k in 1:length(y)]
Y = [k for k in y, j in 1:length(x)]

wind_stress = τ.(Y)
coriolis = f.(X)

gyre_parameters = gyre_Parameters(M,
    g, 
    f0, 
    beta, 
    tau0, 
    H, 
    A_h,
    rho_c, 
    eps_ab, 
    Nx, 
    Ny, 
    dx, 
    dy, 
    dt, 
    τ, 
    f, 
    wind_stress, 
    coriolis
)

# initializing all fields to zero
u_v_eta = gyre(
    zeros(Nx,Ny),
    zeros(Nx,Ny),
    zeros(Nx,Ny),
    zeros(Nx,Ny),
    zeros(Nx,Ny),
    zeros(Nx,Ny),
    zeros(Nx,Ny),
    zeros(Nx,Ny),
    zeros(Nx,Ny)
)

all_states = []
push!(all_states, deepcopy(u_v_eta))

# run the forward problem 
for t = 2:gyre_parameters.M
    advance(u_v_eta, gyre_parameters)
    push!(all_states, deepcopy(u_v_eta))
end

# running the backward problem 

# structure for adjoint variables 
initial_adjoint = deepcopy(u_v_eta)

# setting the initial adjoint variables within structure 
init_adj_eta = 0.0 .* u_v_eta.η
init_adj_eta[31, 31] = 1.0

initial_adjoint.lastu .= 0.
initial_adjoint.lastv .= 0.
initial_adjoint.lastη .= 0.

initial_adjoint.u .= 0.
initial_adjoint.v .= 0. 
initial_adjoint.η = copy(init_adj_eta)

initial_adjoint.lastu_t .= 0.
initial_adjoint.lastv_t .= 0.
initial_adjoint.lastη_t .= 0.

function ad_calc(initial_ad_struct, states, params)

    adjoint_u_v_eta = deepcopy(initial_ad_struct)

    for j = params.M:-1:1 
        current_state = deepcopy(states[j])
        autodiff(advance, Const, Duplicated(current_state, adjoint_u_v_eta), Const(params))
        adjoint_u_v_eta = deepcopy(adjoint_u_v_eta)
    end

    return adjoint_u_v_eta 

end 

ad_u_v_eta = ad_calc(initial_adjoint, all_states, gyre_parameters);
