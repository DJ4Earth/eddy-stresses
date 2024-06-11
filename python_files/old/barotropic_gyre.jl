# forward run with MIT GCM values

using Enzyme, Plots

include("build_structs.jl")
include("advance.jl")

# setting up the grid

# finer grid ---> larger value, courser grid ---> smaller value
# a value of one sets the values to the MIT GCM presets
adjust = 1

# since dt is in seconds this needs to be large 
Nt = 1000           

Nx = 62 * adjust
Ny = 62 * adjust 

dx = (1e3 * 20.0) / adjust     # meters
dy = (1e3 * 20.0) / adjust     # meters

g = 9.81
H = 5000.0
dt = 0.5 * dx / (2 * sqrt(g * H))   # seconds

x = 0:dx:dx*(Nx - 1)
y = 0:dy:dy*(Ny - 1)

X = [j for j in x, k in 1:length(y)]
Y = [k for k in y, j in 1:length(x)]

# initializing parameter values 
f0 = 1e-4
beta = 1e-11
A_h = 400
rho_c = 1000.0
tau0 = 0.1

eps_ab = 0.1 

τ(y) = -tau0 * sin( pi * y / (1200e3) )
f(y) = f0 + beta * y 

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

# starting from rest ---> all initial conditions are zero 

lastu = zeros(Nx, Ny)
lastv = zeros(Nx, Ny)
lastη = zeros(Nx, Ny)


lastu_t = zeros(Nx, Ny)
lastv_t = zeros(Nx, Ny)
lastη_t = zeros(Nx, Ny)

u = zeros(Nx, Ny)
v = zeros(Nx, Ny) 
η = zeros(Nx, Ny)

# creating the structure that will hold the states 
u_v_eta = gyre(
    u,
    v,
    η,
    lastu_t,
    lastv_t, 
    lastη_t
)

for t = 1:M
    advance(u_v_eta, gyre_parameters) 
end

heatmap(u_v_eta.η)
