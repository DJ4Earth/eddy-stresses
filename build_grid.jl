# this script will just contain one function, build_grid, which will initialize
# the structue grid. The structure will contain all the needed information about 
# the Arakawa C-grid that we're solving the equations on 

function build_grid(Lx, Ly, Nx, Ny)

    dx = Lx / Nx 
    dy = Ly / Ny 

    NT = Nx * Ny 

    Nu = (Nx - 1) * Ny 
    Nv = (Ny - 1) * Nx 
    Nq = (Nx + 1) * (Ny + 1)

    x = (dx/2):dx:Lx
    y = dy/2:dy:Ly

    xu = x[1:end-1] .+ dx/2 
    yu = copy(y)
    
    xv = copy(x) 
    yv = y[1:end-1] .+ dy/2 
    
    xq = 0:dx:(Lx + dx/2)
    yq = 0:dy:(Ly + dy/2)

    grid_params = Grid(
    Lx, 
    Ly,
    nx, 
    ny,
    NT,
    Nu, 
    Nv,
    Nq,
    dx, 
    dy
    )
    
    return grid_params

end

function days_to_seconds(Tspinup_days, Trun_days, dt)

    Tspinup_seconds = Int(ceil((Tspinup_days * 24 * 3600) / dt))
    Trun_seconds = Int(ceil((Trun_days * 24 * 3600) / dt))

    return Tspinup_seconds, Trun_seconds

end
