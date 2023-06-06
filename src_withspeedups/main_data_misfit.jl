# making an effort at the whole experiment here, done with Enzyme + checkpointing 
# once this is running and working will structure a main file with it. 

# Writing this here to keep track of the setup: initially, I'm running a low and lower 
# resolution run to debug. The setups will use 

# using Plots, SparseArrays, Parameters, UnPack
# using JLD2, LinearAlgebra
# using Enzyme, Checkpointing, Zygote
# include("init_structs.jl")
# include("init_params.jl")
# include("build_grid.jl")
# include("build_discrete_operators.jl")
# include("advance.jl")
# include("cost_func.jl")
# include("compute_time_deriv.jl")


function setup_data_misfit(
    days, 
    nx, 
    ny, 
    scaling;
    Lx = 3840e3,    
    Ly = 3840e3
)

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

    data_steps = 2:100:T

    data, M = create_data(days, nx, ny, data_steps, scaling)

    states_rhs = SWM_pde(Nu = Nu, 
        Nv = Nv,
        NT = NT, 
        Nq = Nq, 
        T = T
    )

    return grid, params, grad, interp, advec, states_rhs, data, data_steps, M

end


function chkpt_func(
    chkpt_struct::SWM_pde,
    chkpt_scheme::Scheme,
    data::Matrix{Float64},
    data_steps,
    M, 
    grid::Grid,
    params::Params,
    interp::Interps,
    grad::Derivatives,
    advec::Advection
)

    j = 1
    @checkpoint_struct chkpt_scheme chkpt_struct for t in 1:chkpt_struct.T

        advance(chkpt_struct, grid, params, interp, grad, advec)

        if t in data_steps 
            chkpt_struct.J += data_misfit(data[:, j], 
                interp.IuT * chkpt_struct.u0,
                interp.IvT * chkpt_struct.v0,
                chkpt_struct.eta0
            )
            j += 1
        end

        copyto!(chkpt_struct.u, chkpt_struct.u0)
        copyto!(chkpt_struct.v, chkpt_struct.v0)
        copyto!(chkpt_struct.eta, chkpt_struct.eta0)

    end

end

function run_checkpointing_dataex(
    ;days_to_integrate = 30, 
    nx = 50, 
    ny = 50, 
    snaps = 5, 
    scaling = 4
    )

    grid, gyre_params, grad_ops, interp_ops, advec_ops, chkpt_struct, data, data_steps, M = setup_data_misfit(
        days_to_integrate, 
        nx,
        ny, 
        scaling
    )

    snaps = snaps
    verbose = 0
    revolve = Revolve{SWM_pde}(chkpt_struct.T, snaps; verbose=verbose)

    dnu = Zygote.gradient(chkpt_func, 
    chkpt_struct, 
    revolve, 
    data,
    data_steps,
    M, 
    grid, 
    gyre_params, 
    interp_ops, 
    grad_ops, 
    advec_ops
    )

    return dnu

end

# Gradient check with finite differences 

# # for checking that the forward integration still works 

# for t = 1:chkpt_struct.T
#     advance(chkpt_struct, grid, gyre_params, interp_ops, grad_ops, advec_ops)
#     copyto!(chkpt_struct.u, chkpt_struct.u0)
#     copyto!(chkpt_struct.v, chkpt_struct.v0)
#     copyto!(chkpt_struct.eta, chkpt_struct.eta0)
# end