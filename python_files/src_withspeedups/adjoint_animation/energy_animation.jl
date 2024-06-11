# making an effort at the whole experiment here, done with Enzyme + checkpointing 
# once this is running and working will structure a main file with it. 

# This one will try to use Enzyme + checkpointing for an energy sensitivity
# similar to the Burgers equation

using Plots, SparseArrays, Parameters, UnPack
using JLD2, LinearAlgebra
using Enzyme, Checkpointing, Zygote 

include("../init_structs.jl")
include("../init_params.jl")
include("../build_grid.jl")
include("../build_discrete_operators.jl")
include("../advance.jl")
include("../cost_func.jl")
include("../compute_time_deriv.jl")
include("../temp.jl")

Enzyme.API.runtimeActivity!(true)

function setup_energy(
    u0::Vector{Float64}, 
    v0::Vector{Float64}, 
    eta0::Vector{Float64},
    days,
    nx, 
    ny; 
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

    states_rhs = SWM_pde(
        Nu = Nu, 
        Nv = Nv,
        NT = NT,
        Nq = Nq,
        T = T, 
        u = u0,
        v = v0, 
        eta = eta0
    )

    return grid, params, grad, interp, advec, states_rhs

end

function save_states(steps, grid, params, grad, interp, advec, states_rhs)

    u_states = []
    v_states = []
    eta_states = []

    for t in 1:states_rhs.T

        advance(states_rhs, grid, params, interp, grad, advec)

        copyto!(states_rhs.u, states_rhs.u0)
        copyto!(states_rhs.v, states_rhs.v0)
        copyto!(states_rhs.eta, states_rhs.eta0)

        if t in steps 
            push!(u_states, states_rhs.u)
            push!(v_states, states_rhs.v)
            push!(eta_states, states_rhs.eta)
        end

    end

    return u_states, v_states, eta_states, grid, params, grad, interp, advec

end

function chkpt_integration(
    chkpt_struct::SWM_pde,
    chkpt_scheme::Scheme,
    grid::Grid,
    params::Params,
    interp::Interps,
    grad::Derivatives,
    advec::Advection
)

    @checkpoint_struct chkpt_scheme chkpt_struct for t in (chkpt_struct.k * 75 + 1):chkpt_struct.T

        advance(chkpt_struct, grid, params, interp, grad, advec)

        copyto!(chkpt_struct.u, chkpt_struct.u0)
        copyto!(chkpt_struct.v, chkpt_struct.v0)
        copyto!(chkpt_struct.eta, chkpt_struct.eta0)

    end

    energy = sum(chkpt_struct.u.^2 .+ chkpt_struct.v.^2) / (grid.nx * grid.ny)

    return energy

end


function run_checkpointing_energyex(u0, v0, eta0, grid, gyre_params, interp_ops, grad_ops, advec_ops, T, snaps)

    @unpack Nu, Nv, NT, Nq = grid 

    chkpt_struct_outer = SWM_pde(
        Nu = Nu, 
        Nv = Nv,
        NT = NT,
        Nq = Nq,
        T = T, 
        u = u0,
        v = v0, 
        eta = eta0
    )

    snaps = snaps
    verbose = 0
    revolve = Revolve{SWM_pde}(chkpt_struct_outer.T, snaps; verbose=verbose)

    denergy = Zygote.gradient(chkpt_integration, 
        chkpt_struct_outer, 
        revolve, 
        grid, 
        gyre_params, 
        interp_ops, 
        grad_ops, 
        advec_ops
    )

    return denergy

end

function compute_adjoints()

    @load "../initcond_plus_data/states_nx128_ny128_10year_060523.jld2" states_nx128_ny128_10year_060523
    u0 = states_nx128_ny128_10year_060523.u 
    v0 = states_nx128_ny128_10year_060523.v
    eta0 = states_nx128_ny128_10year_060523.eta

    days = 1

    grid, params, grad, interp, advec, states_rhs = setup_energy(
        u0,
        v0,
        eta0,
        days,
        128, 
        128
    )

    total_steps = states_rhs.T

    steps_to_save = 1:75:total_steps

    u_states, v_states, eta_states = save_states(
       steps_to_save, 
       grid,
       params, 
       grad,
       interp, 
       advec,
       states_rhs
    )

    eta_anim = Animation()
    u_anim = Animation()
    v_anim = Animation()

    for j in 1:length(steps_to_save)

        denergy = run_checkpointing_energyex(
            u_states[j],
            v_states[j],
            eta_states[j],
            grid,
            params, 
            interp, 
            grad,
            advec, 
            states_rhs.T,
            5
        )

        frame(eta_anim, heatmap(reshape(denergy[1].eta, 128, 128)', title="Energy sensitivity w.r.t eta"))
        frame(u_anim, heatmap(reshape(denergy[1].u, 127, 128)', title="Energy sensitivity w.r.t u"))
        frame(v_anim, heatmap(reshape(denergy[1].v, 128, 127)', title="Energy sensitivity w.r.t u"))

        states_rhs.k += 1

    end 
    
    gif(eta_anim, "deta_integration_days=$days-060723.gif", fps = 15)
    gif(u_anim, "du_integration_days=$days-060723.gif", fps = 15)
    gif(v_anim, "dv_integration_days=$days-060723.gif", fps = 15)

end


# gradient check with the results from checkpointing - passed

# nx = 128
# ny = 128
# days_to_integrate = 90
# snaps = 1
# @time denergy = run_checkpointing_energyex(days_to_integrate, nx, ny, snaps)

# @load "./initcond_plus_data/states_nx128_ny128_10year_060523.jld2" states_nx128_ny128_10year_060523
# u0 = states_nx128_ny128_10year_060523.u
# v0 = states_nx128_ny128_10year_060523.v
# eta0 = states_nx128_ny128_10year_060523.eta

# @time denergy = run_checkpointing_energyex(u0, v0, eta0, days_to_integrate, nx, ny, snaps)

# du = denergy[1].u
# dv = denergy[1].v
# deta = denergy[1].eta

# steps = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]

# use_to_check = du[88]

# T = days_to_seconds(days_to_integrate, gyre_params.dt)

# grid_, gyre_params, grad_ops, interp_ops, advec_ops, chkpt_struct_outer = setup_energy(
#     days_to_integrate, 
#     nx, 
#     ny
# )

# chkpt_struct_new = SWM_pde(Nu = grid_.Nu, 
#     Nv = grid_.Nv,
#     NT = grid_.NT,
#     Nq = grid_.Nq,
#     T = T
# )

# for t = 1:T
#     advance(chkpt_struct_new, grid_, gyre_params, interp_ops, grad_ops, advec_ops)
#     copyto!(chkpt_struct_new.u, chkpt_struct_new.u0)
#     copyto!(chkpt_struct_new.v, chkpt_struct_new.v0)
#     copyto!(chkpt_struct_new.eta, chkpt_struct_new.eta0)
# end
# energy_to_check = energy(grid_, chkpt_struct_new.u0, chkpt_struct_new.v0)

# diffs = []
# for s in steps 

#     T = days_to_seconds(days_to_integrate, gyre_params.dt)

#     chkpt_struct_new = SWM_pde(Nu = grid_.Nu, 
#         Nv = grid_.Nv,
#         NT = grid_.NT,
#         Nq = grid_.Nq,
#         T = T
#     )

#     chkpt_struct_new.u[88] = s

#     for t = 1:T
#         advance(chkpt_struct_new, grid_, gyre_params, interp_ops, grad_ops, advec_ops)
#         copyto!(chkpt_struct_new.u, chkpt_struct_new.u0)
#         copyto!(chkpt_struct_new.v, chkpt_struct_new.v0)
#         copyto!(chkpt_struct_new.eta, chkpt_struct_new.eta0)
#     end

#     new_energy = energy(grid_, chkpt_struct_new.u0, chkpt_struct_new.v0)

#     push!(diffs, (new_energy - energy_to_check) / s)

# end
