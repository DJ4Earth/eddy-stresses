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

M = 5

eps_ab = 0.1 

τ(y) = -tau0 + cos( pi * y / (1200e3) )
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

gyre_parameters = gyre_Parameters(g, 
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

lastu = zeros(Nx, Ny)
lastv = zeros(Nx, Ny)
lastη = zeros(Nx, Ny)

u = lastu 
v = lastv 
η = lastη

u_v_eta = gyre(
    lastu,
    lastv,
    lastη,
    u,
    v,
    η
)

all_states = []
push!(all_states, deepcopy(u_v_eta))

# run the forward problem 
for t = 2:M
    advance(u_v_eta, gyre_parameters)
    push!(all_states, deepcopy(u_v_eta))
end

# now moving on to running the backward problem 

# structure for adjoint variables 
initial_adjoint = deepcopy(u_v_eta)

# setting the initial adjoint variables within structure 
init_adj_eta = 0.0 .* u_v_eta.η
init_adj_eta[31, 31] = 1.0

initial_adjoint.lastu .= 0.
initial_adjoint.lastv .= 0.
initial_adjoint.lastη = copy(init_adj_eta)

initial_adjoint.u .= 0.
initial_adjoint.v .= 0. 
initial_adjoint.η = copy(init_adj_eta)

function ad_calc(initial_ad_struct, states, params)

    @show initial_ad_struct.η[31, 31]
    adjoint_u_v_eta = deepcopy(initial_ad_struct)
    @show adjoint_u_v_eta.η[31, 31]

    for j = M:-1:1 
        @show sum(abs.(adjoint_u_v_eta.lastη)), states[j].lastη[31, 31]
        current_state = deepcopy(states[j])
        autodiff(advance, Duplicated(current_state, adjoint_u_v_eta), params)
        @show sum(abs.(adjoint_u_v_eta.lastη)), states[j].lastη[31, 31]
        # adjoint_u_v_eta = deepcopy(adjoint_u_v_eta)
    end

    return adjoint_u_v_eta 

end 

# function ad_calc(ad_struct, state, params)

#         @show ad_struct.η[31, 31]
#         @show sum(abs.(ad_struct.η)), state.η[31, 31]

#         autodiff(advance, Duplicated(state, ad_struct), params)

#         @show sum(abs.(ad_struct.η)), state.η[31, 31]
#         @show ad_struct.η[31,31]

#         return ad_struct
    
# end 
    
ad_u_v_eta = ad_calc(initial_adjoint, all_states, gyre_parameters);
