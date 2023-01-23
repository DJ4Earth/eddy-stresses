using Enzyme, Plots

mutable struct gyre

    lastu::Matrix{Float64}
    lastv::Matrix{Float64}
    lastη::Matrix{Float64} 
    
    u::Matrix{Float64}
    v::Matrix{Float64}
    η::Matrix{Float64}

end

mutable struct gyre_Parameters 
    g::Float64              # gravity
    f0::Float64             # Coriolis parameter
    beta::Float64           # Coriolis parameter
    τ_0::Float64             # wind stress amplitude  
    H::Float64              # ocean depth 
    A_h::Float64            # horizontal Laplacian viscosity 
    ρ_c::Float64            # reference density 
    eps_ab::Float64           # parameter pertaining to AB2
    Nx::Int 
    Ny::Int 
    dx::Float64 
    dy::Float64 
    dt::Float64 
    τ::Function              # wind stress function 
    f::Function              # coriolis beta-plane approx 
    wind_stress::Matrix{Float64}    # wind stress values on the grid 
    coriolis::Matrix{Float64}       # coriolis values on the grid 
end

# function to differentiate 
function advance(u_v_eta, parameters) 

    Nx = parameters.Nx 
    Ny = parameters.Ny
    dx = parameters.dx 
    dy = parameters.dy 

    g = parameters.g
    H = parameters.H
    A_h = parameters.A_h 

    eps_ab = parameters.eps_ab

    coriolis = parameters.coriolis
    wind_stress = parameters.wind_stress

    ρ_c = parameters.ρ_c 

    lastu = u_v_eta.lastu
    lastv = u_v_eta.lastv
    lastη = u_v_eta.lastη

    u = u_v_eta.u
    v = u_v_eta.v
    η = u_v_eta.η

    nextu = zeros(Nx, Ny)
    nextv = zeros(Nx, Ny)
    nextη = zeros(Nx, Ny)

    # lastu_t = u_v_eta.lastu_t
    # lastv_t = u_v_eta.lastv_t
    # lastη_t = u_v_eta.lastη_t

    lastu_t = zeros(Nx, Ny)
    lastv_t = zeros(Nx, Ny)
    lastη_t = zeros(Nx, Ny)

    u_t = zeros(Nx, Ny)
    v_t = zeros(Nx, Ny)
    η_t = zeros(Nx, Ny)

    # compute time derivatives 

    for j = 2:Nx - 1 
        for k = 2:Ny - 1

            lastu_t[j,k] = -g * (lastη[j+1,k] - lastη[j-1,k])/(2 * dx) - 
                        (lastu[j+1,k] - lastu[j-1,k])/(2 * dx) * lastu[j,k] - 
                        lastv[j,k] * (lastu[j,k+1] - lastu[j,k-1])/(2 * dy) + 
                        coriolis[j,k] * lastv[j,k] + 
                        wind_stress[j,k]/(ρ_c * H) +
                        A_h * (lastu[j,k+1] - 2 * lastu[j,k] + lastu[j,k-1])/dy^2 +
                        A_h * (lastu[j+1,k] - 2 * lastu[j,k] + lastu[j-1,k])/dx^2

            lastv_t[j,k] = -g * (lastη[j,k+1] - lastη[j,k-1])/(2*dy) -
                        lastu[j,k] * (lastv[j+1,k] - lastv[j-1,k])/(2*dx) -
                        lastv[j,k] * (lastv[j,k+1] - lastv[j,k-1])/(2*dy) -
                        coriolis[j,k] * lastu[j,k] + 
                        A_h * (lastv[j,k+1] - 2 * lastv[j,k] + lastv[j,k-1])/dy^2 +
                        A_h * (lastv[j+1,k] - 2 * lastv[j,k] + lastv[j-1,k])/dx^2

            lastη_t[j,k] = -H * ( (lastv[j,k+1] - lastv[j,k-1])/(2*dy) + (lastu[j+1,k] - lastu[j-1,k])/(2 * dx) )

            u_t[j, k] = -g * (η[j+1,k] - η[j-1,k])/(2 * dx) - 
                        (u[j+1,k] - u[j-1,k])/(2 * dx) *u[j,k] - 
                        v[j,k] * (u[j,k+1] - u[j,k-1])/(2 * dy) + 
                        coriolis[j,k] * v[j,k] + 
                        wind_stress[j,k]/(ρ_c * H) +
                        A_h * (u[j,k+1] - 2 * u[j,k] + u[j,k-1])/dy^2 +
                        A_h * (u[j+1,k] - 2 * u[j,k] + u[j-1,k])/dx^2
            
            v_t[j,k] = -g * (η[j,k+1] - η[j,k-1])/(2*dy) -
                       u[j,k] * (v[j+1,k] - v[j-1,k])/(2*dx) -
                       v[j,k] * (v[j,k+1] - v[j,k-1])/(2*dy) -
                       coriolis[j,k] * u[j,k] + 
                       A_h * (v[j,k+1] - 2 * v[j,k] + v[j,k-1])/dy^2 +
                       A_h * (v[j+1,k] - 2 * v[j,k] + v[j-1,k])/dx^2

            η_t[j,k] = -H * ( (v[j,k+1] - v[j,k-1])/(2*dy) + (u[j+1,k] - u[j-1,k])/(2*dx) )

        end
    end

    # update fields 

    for j = 2:Nx - 1
        for k = 2:Ny - 1 
            
            nextu[j,k] = u[j,k] + dt * (3/2 + eps_ab) * u_t[j,k] -
                         dt * (1/2 + eps_ab) *lastu_t[j,k]

            nextv[j,k] = v[j,k] + dt * (3/2 + eps_ab) * v_t[j,k] -
                         dt * (1/2 + eps_ab) *lastv_t[j,k]

            nextη[j,k] = η[j,k] + dt * (3/2 + eps_ab) * η_t[j,k] -
                         dt * (1/2 + eps_ab) * lastη_t[j,k]

        end
    end

    # apply vertical boundary conditions
    for j = 1:Ny 

        nextu[1,j] = nextu[2,j]
        nextv[1,j] = 0.0
        nextη[1,j] = 0.0 

        nextu[Nx,j] = nextu[Nx-1,j]
        nextv[Nx,j] = 0.0
        nextη[Nx,j] = 0.0 

    end

    # apply horizontal boundary conditions 
    for k = 1:Nx 

        nextu[k,1] = 0.0
        nextv[k,1] = nextv[k,2] 
        nextη[k,1] = 0.0

        nextu[k,Ny] = 0.0
        nextv[k,Ny] = nextv[k,Ny-1]
        nextη[k,Ny] = 0.0

    end

    # (hopefully) cycle the states 

    # lastu_t = u_t 
    # lastv_t = v_t 
    # lastη_t = η_t

    u_v_eta.lastu = u
    u_v_eta.lastv = v
    u_v_eta.lastη = η

    u_v_eta.u = nextu
    u_v_eta.v = nextv 
    u_v_eta.η = nextη

    return nothing 

end 

# setting all the different parameters 
g = 9.81 
f0 = 1e-4
beta = 1e-11
H = 5000.0
A_h = 400
rho_c = 1000.0
tau0 = 0.1

M = 10

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

# structure of initial states 
u_v_eta = gyre(
    zeros(Nx, Ny),
    zeros(Nx, Ny),
    zeros(Nx, Ny),
    zeros(Nx, Ny),
    zeros(Nx, Ny),
    zeros(Nx, Ny)
)

all_states = []
push!(all_states, u_v_eta)

# run the forward problem and store intermediate states (not the ideal storage method probably, will come back to this)
for t = 1:M
    advance(u_v_eta, gyre_parameters)
    push!(all_states, u_v_eta)
end

# now moving on to running the backward problem 

# structure for adjoint variables 
adjoint_u_v_eta = deepcopy(u_v_eta)

# setting the initial adjoint variables within structure 
init_adj = 0.0 .* u_v_eta.η
init_adj[31, 31] = 1.0

adjoint_u_v_eta.lastu .= 0.
adjoint_u_v_eta.lastv .= 0.
adjoint_u_v_eta.lastη .= 0. 

adjoint_u_v_eta.u .= 0.
adjoint_u_v_eta.v .= 0. 
adjoint_u_v_eta.η = init_adj

autodiff(advance, Const, Duplicated(all_states[M], adjoint_u_v_eta), Const(gyre_parameters))
adjoint_u_v_eta = copy(adjoint_u_v_eta)
