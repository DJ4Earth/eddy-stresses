# This script will contain one function which calculates all of the numerical 
# derivatives we need and then advances the model forward a step.  

# for now I'm leaving the timestep as a second order AB method

function advance(u_v_eta, parameters, grid) 

    Nx = grid.Nx 
    Ny = grid.Ny
    dx = grid.dx 
    dy = grid.dy 
    dt = grid.dt

    g = parameters.g
    H = parameters.H
    A_h = parameters.A_h 

    eps_ab = parameters.eps_ab

    coriolis = parameters.coriolis
    wind_stress = parameters.wind_stress

    ρ_c = parameters.ρ_c 

    u = copy(u_v_eta.u)
    v = copy(u_v_eta.v)
    η = copy(u_v_eta.η)

    lastu_t = copy(u_v_eta.lastu_t)
    lastv_t = copy(u_v_eta.lastv_t)
    lastη_t = copy(u_v_eta.lastη_t)

    nextu = zeros(Nx, Ny)
    nextv = zeros(Nx, Ny)
    nextη = zeros(Nx, Ny)

    u_t = zeros(Nx, Ny)
    v_t = zeros(Nx, Ny)
    η_t = zeros(Nx, Ny)

    # compute time derivatives 

    for j = 2:Nx - 1 
        for k = 2:Ny - 1

            # lastu_t[j,k] = -g * (lastη[j+1,k] - lastη[j-1,k])/(2 * dx) - 
            #             (lastu[j+1,k] - lastu[j-1,k])/(2 * dx) * lastu[j,k] - 
            #             lastv[j,k] * (lastu[j,k+1] - lastu[j,k-1])/(2 * dy) + 
            #             coriolis[j,k] * lastv[j,k] + 
            #             wind_stress[j,k]/(ρ_c * H) +
            #             A_h * (lastu[j,k+1] - 2 * lastu[j,k] + lastu[j,k-1])/dy^2 +
            #             A_h * (lastu[j+1,k] - 2 * lastu[j,k] + lastu[j-1,k])/dx^2

            # lastv_t[j,k] = -g * (lastη[j,k+1] - lastη[j,k-1])/(2*dy) -
            #             lastu[j,k] * (lastv[j+1,k] - lastv[j-1,k])/(2*dx) -
            #             lastv[j,k] * (lastv[j,k+1] - lastv[j,k-1])/(2*dy) -
            #             coriolis[j,k] * lastu[j,k] + 
            #             A_h * (lastv[j,k+1] - 2 * lastv[j,k] + lastv[j,k-1])/dy^2 +
            #             A_h * (lastv[j+1,k] - 2 * lastv[j,k] + lastv[j-1,k])/dx^2

            # lastη_t[j,k] = -H * ( (lastv[j,k+1] - lastv[j,k-1])/(2*dy) + (lastu[j+1,k] - lastu[j-1,k])/(2 * dx) )

            u_t[j, k] = -g * (η[j+1,k] - η[j-1,k])/(2 * dx) - 
                        (u[j+1,k] - u[j-1,k])/(2 * dx) * u[j,k] - 
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

            η_t[j,k] = -H * ( (v[j,k+1] - v[j,k-1])/(2*dy) + (u[j+1,k] - u[j-1,k])/(2 * dx) )

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

    # cycle the states 

    # copyto!(u_v_eta.lastu, u)
    # copyto!(u_v_eta.lastv, v)
    # copyto!(u_v_eta.lastη, η)

    copyto!(u_v_eta.u, nextu)
    copyto!(u_v_eta.v, nextv) 
    copyto!(u_v_eta.η, nextη)

    copyto!(u_v_eta.lastu_t, u_t)
    copyto!(u_v_eta.lastv_t, v_t)
    copyto!(u_v_eta.lastη_t, η_t)

    return nothing 

end 