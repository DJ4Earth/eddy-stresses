# random functions written with the sole purpose being to save states/energy/data etc. 


function save_energy_spectrum(u0, v0, eta0, days, nx, ny; Lx = 3840e3, Ly = 3840e3)

    grid = build_grid(Lx, Ly, nx, ny)
    params = def_params(grid)

    # building discrete operators
    grad = build_derivs(grid)            # discrete gradient operators
    interp = build_interp(grid, grad)    # discrete interpolation operators (travels between grids)
    advec = build_advec(grid)

    Nu = grid.Nu
    Nv = grid.Nv
    NT = grid.NT
    Nq = grid.Nq 

    T = days_to_seconds(days, params.dt)

    states_rhs = SWM_pde(Nu = Nu, 
        Nv = Nv,
        NT = NT, 
        Nq = Nq, 
        T = T,
        u = u0,
        v = v0,
        eta = eta0,
        nu = params.nu
    )

    # energy_spectra_u = []
    # energy_spectra_v = []
    energy_spectra = []
    
    for t in 1:T

        advance(states_rhs, grid, params, interp, grad, advec)

        temp1 = fft(reshape(states_rhs.u.^2, grid.nx-1, grid.ny)')
        temp2 = fft(reshape(states_rhs.v.^2, grid.nx, grid.ny-1)')

        # push!(energy_spectra_u, temp1)
        # push!(energy_spectra_v, temp2)
        push!(energy_spectra, temp1[1:end-1, :] + temp2[:, 1:end-1])

        copyto!(states_rhs.u, states_rhs.u0)
        copyto!(states_rhs.v, states_rhs.v0)
        copyto!(states_rhs.eta, states_rhs.eta0)

    end 
        
    return energy_spectra

end

function save_energy_plusend(days, nx, ny; Lx = 3840e3, Ly = 3840e3)                 

    grid = build_grid(Lx, Ly, nx, ny)
    params = def_params(grid)

    # building discrete operators
    grad = build_derivs(grid)            # discrete gradient operators
    interp = build_interp(grid, grad)    # discrete interpolation operators (travels between grids)
    advec = build_advec(grid)

    Nu = grid.Nu
    Nv = grid.Nv
    NT = grid.NT
    Nq = grid.Nq 

    T = days_to_seconds(days, params.dt)

    states_rhs = SWM_pde(Nu = Nu, 
        Nv = Nv,
        NT = NT, 
        Nq = Nq, 
        T = T
    )

    # all_u_states = zeros(Nu, T)
    # all_v_states = zeros(Nv, T)
    # all_eta_states = zeros(NT, T)

    energy = zeros(T)
    
    @time for t in 1:T

        advance(states_rhs, grid, params, interp, grad, advec)

        energy[t] = sum(states_rhs.u.^2 .+ states_rhs.v.^2) / (grid.nx * grid.ny)

        copyto!(states_rhs.u, states_rhs.u0)
        copyto!(states_rhs.v, states_rhs.v0)
        copyto!(states_rhs.eta, states_rhs.eta0)

    end 
        
    return energy, gyre_vector(states_rhs.u, states_rhs.v, states_rhs.eta)

end

function save_energy(days, nx, ny; Lx = 3840e3, Ly = 3840e3)                 

    grid = build_grid(Lx, Ly, nx, ny)
    params = def_params(grid)

    # building discrete operators
    grad = build_derivs(grid)            # discrete gradient operators
    interp = build_interp(grid, grad)    # discrete interpolation operators (travels between grids)
    advec = build_advec(grid)

    Nu = grid.Nu
    Nv = grid.Nv
    NT = grid.NT
    Nq = grid.Nq 

    T = days_to_seconds(days, params.dt)

    states_rhs = SWM_pde(Nu = Nu, 
        Nv = Nv,
        NT = NT, 
        Nq = Nq, 
        T = T
    )

    # all_u_states = zeros(Nu, T)
    # all_v_states = zeros(Nv, T)
    # all_eta_states = zeros(NT, T)

    energy = zeros(T)
    
    for t in 1:T

        advance(states_rhs, grid, params, interp, grad, advec)

        energy[t] = sum(states_rhs.u.^2 .+ states_rhs.v.^2) / (grid.nx * grid.ny)

        copyto!(states_rhs.u, states_rhs.u0)
        copyto!(states_rhs.v, states_rhs.v0)
        copyto!(states_rhs.eta, states_rhs.eta0)

    end 
        
    return energy

end

function save_energy(u0, v0, eta0, days, nx, ny; Lx = 3840e3, Ly = 3840e3)                 

    grid = build_grid(Lx, Ly, nx, ny)
    params = def_params(grid)

    # building discrete operators
    grad = build_derivs(grid)            # discrete gradient operators
    interp = build_interp(grid, grad)    # discrete interpolation operators (travels between grids)
    advec = build_advec(grid)

    Nu = grid.Nu
    Nv = grid.Nv
    NT = grid.NT
    Nq = grid.Nq 

    T = days_to_seconds(days, params.dt)

    states_rhs = SWM_pde(Nu = Nu, 
        Nv = Nv,
        NT = NT, 
        Nq = Nq, 
        T = T,
        u = u0,
        v = v0,
        eta = eta0
    )

    # all_u_states = zeros(Nu, T)
    # all_v_states = zeros(Nv, T)
    # all_eta_states = zeros(NT, T)

    energy = zeros(T)
    
    for t in 1:T

        advance(states_rhs, grid, params, interp, grad, advec)

        energy[t] = sum(states_rhs.u.^2 .+ states_rhs.v.^2) / (grid.nx * grid.ny)

        copyto!(states_rhs.u, states_rhs.u0)
        copyto!(states_rhs.v, states_rhs.v0)
        copyto!(states_rhs.eta, states_rhs.eta0)

    end 
        
    return energy

end


function save_states_and_energy(u0, v0, eta0, days, nx, ny; Lx = 3840e3, Ly = 3840e3)                 

    grid = build_grid(Lx, Ly, nx, ny)
    params = def_params(grid)

    # building discrete operators
    grad = build_derivs(grid)            # discrete gradient operators
    interp = build_interp(grid, grad)    # discrete interpolation operators (travels between grids)
    advec = build_advec(grid)

    Nu = grid.Nu
    Nv = grid.Nv
    NT = grid.NT
    Nq = grid.Nq 

    T = days_to_seconds(days, params.dt)

    states_rhs = SWM_pde(Nu = Nu, 
        Nv = Nv,
        NT = NT, 
        Nq = Nq, 
        T = T,
        u = u0,
        v = v0,
        eta = eta0
    )

    energy = zeros(T)

    every_hour = 10:10:T
    all_u_states = zeros(grid.Nu, length(every_hour))
    all_v_states = zeros(grid.Nv, length(every_hour))

    k = 1    
    for t in 1:T

        advance(states_rhs, grid, params, interp, grad, advec)

        energy[t] = sum(states_rhs.u.^2 .+ states_rhs.v.^2) / (grid.nx * grid.ny)

        if t in every_hour
            all_u_states[:, k] .= states_rhs.u0
            all_v_states[:, k] .= states_rhs.v0
            k += 1
        end

        copyto!(states_rhs.u, states_rhs.u0)
        copyto!(states_rhs.v, states_rhs.v0)
        copyto!(states_rhs.eta, states_rhs.eta0)

    end 
        
    return energy, all_u_states, all_v_states

end

function save_energy_spectra(u0, v0, eta0, days, nx, ny; Lx = 3840e3, Ly = 3840e3)

    grid = build_grid(Lx, Ly, nx, ny)
    params = def_params(grid)

    # building discrete operators
    grad = build_derivs(grid)            # discrete gradient operators
    interp = build_interp(grid, grad)    # discrete interpolation operators (travels between grids)
    advec = build_advec(grid)

    Nu = grid.Nu
    Nv = grid.Nv
    NT = grid.NT
    Nq = grid.Nq 

    T = days_to_seconds(days, params.dt)

    states_rhs = SWM_pde(Nu = Nu, 
        Nv = Nv,
        NT = NT, 
        Nq = Nq, 
        T = T,
        u = u0,
        v = v0,
        eta = eta0
    )

end

function create_averaging_op(nx_lowres, ny_lowres, data_steps; Lx = 3840e3, Ly = 3840e3, scaling = 4)

    nx_highres = nx_lowres * scaling
    ny_highres = ny_lowres * scaling 
    
    grid_lowres = build_grid(Lx, Ly, nx_lowres, ny_lowres)
    
    grid_highres = build_grid(Lx, Ly, nx_highres, ny_highres)
    
    Trun = days_to_seconds(days, params.dt)
    
    Nu = grid_highres.Nu
    Nv = grid_highres.Nv
    NT = grid_highres.NT
    Nq = grid_highres.Nq 

    diag1 = (1 / scaling^2) .* ones(NT)
    M = spdiagm(NT, NT, 0.0 .* diag1)
    for k = 1:scaling
        for j = 1:scaling
            M += spdiagm(NT, NT, j+(k-1)*nx_highres - 1 => diag1[1:end-j-(k-1)*nx_highres + 1])
        end
    end
    M = M[1:scaling:end, :]
    index1 = 1:nx_lowres*ny_lowres*scaling
    for k in 1:ny_lowres
        index1 = filter(x -> x ∉ [j for j in ((k-1)*nx_lowres*scaling+nx_lowres+1):((k-1)*nx_lowres*scaling+nx_lowres*scaling)], index1)
    end
    M = M[index1, :]

    return M 

end

# This function needs to be given 
#           days - how many days to integrate the model for
#           nx_lowres, ny_lowres - grid resolution (number of cells in the x and y directions
#                    respectively) of the courser grid
#           Lx, Ly - size of the domain, have a default value but can set 
#                    manually if needed 
#           data_steps - which timesteps we want to store data at
# Theoretically, we should only ever run this function *once*, from there 
# the data points will be stored as a JLD2 data file 
function create_data(days, nx_lowres, ny_lowres, data_steps, scaling; Lx = 3840e3, Ly = 3840e3)

    nx_highres = nx_lowres * scaling
    ny_highres = ny_lowres * scaling 
    
    grid_lowres = build_grid(Lx, Ly, nx_lowres, ny_lowres)
    
    grid_highres = build_grid(Lx, Ly, nx_highres, ny_highres)
    params = def_params(grid_highres)
    
    # building discrete operators
    grad = build_derivs(grid_highres)            # discrete gradient operators
    interp = build_interp(grid_highres, grad)    # discrete interpolation operators (travels between grids)
    advec = build_advec(grid_highres)
    
    Trun = days_to_seconds(days, params.dt)
    
    Nu = grid_highres.Nu
    Nv = grid_highres.Nv
    NT = grid_highres.NT
    Nq = grid_highres.Nq 
    
    u_v_eta_rhs = SWM_pde(Nu = Nu, 
        Nv = Nv,
        NT = NT, 
        Nq = Nq
    )
    
    # In order to compare high res data to low res velocities I'm going to 
    # (1) interpolate the velocities to the T-grid (cell centers) 
    # (2) average the high res data points down to the low res grid 
    # (3) interpolate the low res results to the cell centers. 
    # Then I'll be comparing apples to apples (hopefully) (will check with Patrick that this is a valid method)
    
    # Building the averaging operator needed for step (2) above
    diag1 = (1 / scaling^2) .* ones(NT)
    M = spdiagm(NT, NT, 0.0 .* diag1)
    for k = 1:scaling
        for j = 1:scaling
            M += spdiagm(NT, NT, j+(k-1)*nx_highres - 1 => diag1[1:end-j-(k-1)*nx_highres + 1])
        end
    end
    M = M[1:scaling:end, :]
    index1 = 1:nx_lowres*ny_lowres*scaling
    for k in 1:ny_lowres
        index1 = filter(x -> x ∉ [j for j in ((k-1)*nx_lowres*scaling+nx_lowres+1):((k-1)*nx_lowres*scaling+nx_lowres*scaling)], index1)
    end
    M = M[index1, :]
    
    # initializing where to store the data 
    data = zeros(3 * grid_lowres.NT, length(data_steps))
    
    # the steps where we want data in the high res model correspond to (roughly) scaling * t for t 
    # in the low res model. for simplicity I'm going to keep the times in the low res where I want to have 
    # data and then just scale them in the for loop to find corresponding high res data points  
    
    # for an initial effort I'm just going to run pretty course resolution models for both the high and low res 
    
    if 1 in scaling .* data_steps 
        data[:, 1] .= [M * (interp.IuT * u_v_eta_rhs.u); 
            M * (interp.IvT * u_v_eta_rhs.v); 
            M * u_v_eta_rhs.eta
        ]
        j = 2
    else
        j = 1
    end
    
    for t in 2:Trun
    
        advance(u_v_eta_rhs, grid_highres, params, interp, grad, advec) 
    
        if t in scaling .* data_steps 
            data[:, j] .= [M * (interp.IuT * u_v_eta_rhs.u); 
            M * (interp.IvT * u_v_eta_rhs.v); 
            M * u_v_eta_rhs.eta
        ]
            j += 1
        end
    
        copyto!(u_v_eta_rhs.u, u_v_eta_rhs.u0)
        copyto!(u_v_eta_rhs.v, u_v_eta_rhs.v0)
        copyto!(u_v_eta_rhs.eta, u_v_eta_rhs.eta0)
    
    end
    
    return data, M
    
end

function plot_energy(which_step, nx, ny)

    u_states = load_object("../data_files/u_afterspinup_nxny128_2years_hourlysaves_062823.jld2")
    v_states = load_object("../data_files/v_afterspinup_nxny128_2years_hourlysaves_062823.jld2")

    u = u_states[:, which_step]
    v = v_states[:, which_step]

    grid = build_grid(3840e3, 3840e3, nx, ny)
    grad = build_derivs(grid)         
    interp = build_interp(grid, grad)

    u_in_center = interp.IuT * u
    v_in_center = interp.IvT * v 
    
    energy = u_in_center.^2 + v_in_center.^2 

    return reshape(energy, nx, ny)';

end