using Plots, SparseArrays, Parameters, UnPack
using JLD2, LinearAlgebra

include("init_structs.jl")
include("init_params.jl")
include("build_grid.jl")
include("build_discrete_operators.jl")
include("advance.jl")
include("compute_time_deriv.jl")

function create_gif(days, nx, ny; Lx = 3840e3, Ly = 3840e3)                 

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
    
    anim = @animate for t in 1:T

        advance(states_rhs, grid, params, interp, grad, advec)
        heatmap(reshape(states_rhs.eta, ny, nx)', xlabel="x", ylabel="y", title="Timestep = $t", clim=(-2, 2))
        copyto!(states_rhs.u, states_rhs.u0)
        copyto!(states_rhs.v, states_rhs.v0)
        copyto!(states_rhs.eta, states_rhs.eta0)
        
    end every 100

    gif(anim, "test.gif", fps = 10)
        
    return nothing

end

function create_gif(u0, v0, eta0, days, nx, ny; Lx = 3840e3, Ly = 3840e3)                 

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
    
    anim = @animate for t in 1:T

        advance(states_rhs, grid, params, interp, grad, advec)
        heatmap(reshape(states_rhs.eta, ny, nx)', xlabel="x", ylabel="y", title="eta", clim=(-3.75, 3.75))        
        copyto!(states_rhs.u, states_rhs.u0)
        copyto!(states_rhs.v, states_rhs.v0)
        copyto!(states_rhs.eta, states_rhs.eta0)
        
    end every 75

    gif(anim, "eta_integration_6month.gif", fps = 15)
        
    return nothing

end

nx = 128
ny = 128
days_to_integrate = 180

@load "./initcond_plus_data/states_nx128_ny128_10year_060523.jld2" states_nx128_ny128_10year_060523
u0 = states_nx128_ny128_10year_060523.u
v0 = states_nx128_ny128_10year_060523.v
eta0 = states_nx128_ny128_10year_060523.eta

create_gif(u0, v0, eta0, days_to_integrate, nx, ny)

