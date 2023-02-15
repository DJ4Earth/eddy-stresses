# this function will define all of the desired parameters and return them in a structure for use in other locations
function def_params(grid)

    nx = grid.nx
    ny = grid.ny
    dx = grid.dx 
    dy = grid.dx
    Lx = grid.Lx
    Ly = grid.Ly

    x = (dx/2):dx:Lx
    y = dy/2:dy:Ly

    xu = x[1:end-1] .+ dx/2 
    yu = copy(y)
    
    xv = copy(x) 
    yv = y[1:end-1] .+ dy/2 
    
    xq = 0:dx:(Lx + dx/2)
    yq = 0:dy:(Ly + dy/2)

    # our beta-plane approximation is centered at 30degrees (is this the correct way to phrase this?)
    omega = 2 * pi / (24 * 3600)
    R = 6.371e6
    f0 = 2 * omega * sin(30 * pi / 180)
    beta = (2 * omega / R) * cos(30 * pi / 180)

    # viscosity coefficient, depends on the grid size 
    A_h = (128*540/(min(nx, ny))) * max(dx, dy)^2       # viscosity coefficient [meters^2 / second]
    rho_c = 1000.0  

    # bottom drag coefficient
    bottom_drag = 1e-5

    # building the vectors that will contain the coriolis force and wind stress (techincally forcings and not parameters I know)
    # come back to this, want to understand 
    Yq_shift = yq .- Ly/2
    Yq = vec([k for k in Yq_shift, j in 1:nx+1]')
    f(y) = f0 + beta * y
    coriolis = vec(f.(Yq)')

    Yu = vec([k for k in yu, j in 1:length(xu)]')
    ws_cos = cos.(2 * pi .* ((Yu .- Ly/2)./Ly) )
    ws_sin = 2 .* sin.(pi .* ((Yu .- Ly/2)./Ly)  ) 
    wind_stress = 0.12 .* (ws_cos + ws_sin) ./rho_c

    g = 9.81    # gravity [meters^2 / second]
    H = 500.0   # depth of the box [meters]

    dt = (0.9 * min(dx, dy)) / (sqrt(g * H))   # CFL condition for dt [seconds]
    
    gyre_params = Parameters(
    dt,
    g, 
    f0, 
    beta, 
    H, 
    A_h,
    rho_c, 
    bottom_drag,  
    wind_stress, 
    coriolis
    )

    return gyre_params 

end