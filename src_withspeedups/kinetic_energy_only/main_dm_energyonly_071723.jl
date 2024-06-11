# making an effort at the whole experiment here, done with Enzyme + checkpointing 
# once this is running and working will structure a main file with it. 

# Writing this here to keep track of the setup: initially, I'm running a low and lower 
# resolution run to debug. The setups will use 

using Plots, SparseArrays, Parameters, UnPack
using JLD2, LinearAlgebra
using Enzyme, Checkpointing, Zygote
using Optim

include("../init_structs.jl")
include("../init_params.jl")
include("../build_grid.jl")
include("../build_discrete_operators.jl")
include("../advance.jl")
# include("cost_func.jl")
include("../compute_time_deriv.jl")


function setup_data_misfit(
    days 
    ;
    nx = 128, 
    ny = 128,
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

    data_steps = [k for k in 1:T] 

    data = load_object("../../data_files/energy_afterspinup_nxny512_2years_everystep_062823.jld2")

    # data, M = create_data(days, nx, ny, data_steps, scaling)

    states_rhs = SWM_pde(Nu = Nu, 
        Nv = Nv,
        NT = NT, 
        Nq = Nq, 
        T = T,
        nu = copy(params.nu)
    )

    return grid, params, grad, interp, advec, states_rhs, data_steps, data 

end


function setup_data_misfit(u0, 
    v0, 
    eta0,   
    days,
    nu 
    ;
    nx = 128, 
    ny = 128,
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

    data_steps = [k for k in 1:T]

    data = load_object("../../data_files/energy_afterspinup_nxny512_2years_everystep_062823.jld2")

    # data, M = create_data(days, nx, ny, data_steps, scaling)

    states_rhs = SWM_pde(Nu = Nu, 
        Nv = Nv,
        NT = NT, 
        Nq = Nq,
        u = u0,
        v = v0,
        eta = eta0, 
        T = T,
        nu = nu
    )

    return grid, params, grad, interp, advec, states_rhs, data_steps, data 

end


function chkpt_func(
    chkpt_struct::SWM_pde,
    chkpt_scheme::Scheme,
    data,
    data_steps,
    grid::Grid,
    params::Params,
    interp::Interps,
    grad::Derivatives,
    advec::Advection
    )

    @checkpoint_struct chkpt_scheme chkpt_struct for chkpt_struct.t in 1:chkpt_struct.T

        advance(chkpt_struct, grid, params, interp, grad, advec)

        if chkpt_struct.t in data_steps 
            chkpt_struct.J += sum(chkpt_struct.u.^2 .+ chkpt_struct.v.^2) / (grid.nx * grid.ny) - data[4 * chkpt_struct.j]
            chkpt_struct.j += 1
        end

        copyto!(chkpt_struct.u, chkpt_struct.u0)
        copyto!(chkpt_struct.v, chkpt_struct.v0)
        copyto!(chkpt_struct.eta, chkpt_struct.eta0)

    end

    return chkpt_struct.J

end

function run_checkpointing_dataex_energyonly(nu 
    ;Ndays = 1
    )

    # non-zero initial condition 
    init_states = load_object("../../data_files/states_nx128_ny128_10year_060523.jld2")
    # grid, gyre_params, grad_ops, interp_ops, advec_ops, chkpt_struct, data_steps, data  = setup_data_misfit(init_states.u,
    # init_states.v,
    # init_states.eta,
    # Ndays, 
    # nu
    # )

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

    data_steps = [k for k in 1:T]

    data = load_object("../../data_files/energy_afterspinup_nxny512_2years_everystep_062823.jld2")

    # data, M = create_data(days, nx, ny, data_steps, scaling)

    states_rhs = SWM_pde(Nu = Nu, 
        Nv = Nv,
        NT = NT, 
        Nq = Nq,
        u = init_states.u,
        v = init_states.v,
        eta = init_states.eta, 
        T = T,
        nu = nu
    )

    # from rest initial condition 
    # grid, gyre_params, grad_ops, interp_ops, advec_ops, chkpt_struct, data_steps, data  = setup_data_misfit(
    # Ndays
    # )

    snaps = Int(floor(sqrt(chkpt_struct.T)))
    revolve = Revolve{SWM_pde}(chkpt_struct.T, snaps; verbose=1, gc=true, write_checkpoints=false)

    dnu = Zygote.gradient(chkpt_func, 
    chkpt_struct, 
    revolve, 
    data,
    data_steps,
    grid, 
    gyre_params, 
    interp_ops, 
    grad_ops, 
    advec_ops
    )

    return dnu

end

function for_optim!(F, G, x)

    dnu, chkpt_struct = run_checkpointing_dataex_energyonly()

    if G !== nothing 

        G .= dnu[1].nu 

    end 

    if F !== nothing 

        print(chkpt_struct.J)

        return chkpt_struct.J

    end 

end 

# init_states = load_object("../../data_files/states_nx128_ny128_10year_060523.jld2")

# grid, gyre_params, grad_ops, interp_ops, advec_ops, chkpt_struct, data_steps, data  = setup_data_misfit(init_states.u,
# init_states.v,
# init_states.eta,
# 30
# )

# A_h = compute_viscosity(grid.nx, grid.ny, grid.dx, grid.dy)
# nu = A_h .* ones(grid.NT)                                            # placing the viscosity coefficient on the tracer grid (cell centers)

# # res = Optim.optimize(Optim.only_fg!(for_optim!), nu, Optim.LBFGS(), Optim.Options(iterations = 30))

# dnu = run_checkpointing_dataex_energyonly(days_to_integrate = 30)

# Gradient check with finite differences 

function gradient_check()


    days_to_integrate = 45
    @time dnu = run_checkpointing_dataex_energyonly(days_to_integrate = days_to_integrate)

    # du = dnu[1].u
    # dv = dnu[1].v
    # deta = dnu[1].eta

    d_param = copy(dnu[1].nu)

    steps = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]

    want_to_converge_to = d_param[8]

    grid_, gyre_params, grad_ops, interp_ops, advec_ops, chkpt_struct_outer, data_steps, data  = setup_data_misfit(
        days_to_integrate
    )

    T = days_to_seconds(days_to_integrate, gyre_params.dt)

    chkpt_struct_new = SWM_pde(Nu = grid_.Nu, 
        Nv = grid_.Nv,
        NT = grid_.NT,
        Nq = grid_.Nq,
        T = T,
        nu = gyre_params.nu
    )

    for t = 1:T

        advance(chkpt_struct_new, grid_, gyre_params, interp_ops, grad_ops, advec_ops)

        if t in data_steps 
            chkpt_struct_new.J += sum(chkpt_struct_new.u.^2 .+ chkpt_struct_new.v.^2) / (grid_.nx * grid_.ny) - data[4 * chkpt_struct_new.j]
            chkpt_struct_new.j += 1
        end

        copyto!(chkpt_struct_new.u, chkpt_struct_new.u0)
        copyto!(chkpt_struct_new.v, chkpt_struct_new.v0)
        copyto!(chkpt_struct_new.eta, chkpt_struct_new.eta0)

    end

    cost_tocheck = copy(chkpt_struct_new.J)

    @show cost_tocheck

    diffs = []
    for s in steps 

        T = days_to_seconds(days_to_integrate, gyre_params.dt)

        chkpt_struct_new2 = SWM_pde(Nu = grid_.Nu, 
            Nv = grid_.Nv,
            NT = grid_.NT,
            Nq = grid_.Nq,
            T = T,
            nu = gyre_params.nu
        )

        chkpt_struct_new2.nu .= copy(gyre_params.nu)

        @show s 
        @show chkpt_struct_new2.nu[8]

        chkpt_struct_new2.nu[8] += s

        @show chkpt_struct_new2.nu[8]

        for t = 1:T
            advance(chkpt_struct_new2, grid_, gyre_params, interp_ops, grad_ops, advec_ops)

            if t in data_steps 
                chkpt_struct_new2.J += sum(chkpt_struct_new2.u.^2 .+ chkpt_struct_new2.v.^2) / (grid_.nx * grid_.ny) - data[4 * chkpt_struct_new2.j]
                chkpt_struct_new2.j += 1
            end

            copyto!(chkpt_struct_new2.u, chkpt_struct_new2.u0)
            copyto!(chkpt_struct_new2.v, chkpt_struct_new2.v0)
            copyto!(chkpt_struct_new2.eta, chkpt_struct_new2.eta0)
        end

        @show chkpt_struct_new2.J

        push!(diffs, (chkpt_struct_new2.J - cost_tocheck) / s)

    end

    return diffs, want_to_converge_to

end

diffs, want_to_converge_to = gradient_check()

# # for checking that the forward integration still works 

# for t = 1:chkpt_struct.T
#     advance(chkpt_struct, grid, gyre_params, interp_ops, grad_ops, advec_ops)
#     copyto!(chkpt_struct.u, chkpt_struct.u0)
#     copyto!(chkpt_struct.v, chkpt_struct.v0)
#     copyto!(chkpt_struct.eta, chkpt_struct.eta0)
# end