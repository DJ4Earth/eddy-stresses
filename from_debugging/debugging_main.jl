using Enzyme, Plots, SparseArrays, Parameters

include("init_structs.jl")
include("init_params.jl")
include("build_grid.jl")
include("build_discrete_operators.jl")
include("allocate_rhs.jl")

function advance_check(u_v_eta, rhs, params, grid, interp, grad, advec) 

    dt = params.dt

    # we now use RK4 as the timestepper, here I'm storing the coefficients needed for this 
    rk_a = [1/6, 1/3, 1/3, 1/6]
    rk_b = [1/2, 1/2, 1]

    u = copy(u_v_eta.u)
    v = copy(u_v_eta.v)
    eta = copy(u_v_eta.eta)

    rhs.u0 .= u
    rhs.v0 .= v
    rhs.eta0 .= eta

    rhs.u1 .= u
    rhs.v1 .= v
    rhs.eta1 .= eta

    for j in 1:4

        comp_u_v_eta_t_check(rhs, params, grid.nx, interp, grad, advec)

        if j < 4
            rhs.u1 .= u + rk_b[j] * dt * rhs.u_t
            rhs.v1 .= v + rk_b[j] * dt * rhs.v_t
            rhs.eta1 .= eta + rk_b[j] * dt * rhs.eta_t
        end

        rhs.u0 .= rhs.u0 + rk_a[j] * dt * rhs.u_t
        rhs.v0 .= rhs.v0 + rk_a[j] * dt * rhs.v_t 
        rhs.eta0 .= rhs.eta0 + rk_a[j] * dt * rhs.eta_t 

    end

    @assert all(x -> x < 10.0, rhs.u0)
    @assert all(x -> x < 10.0, rhs.v0)
    @assert all(x -> x < 10.0, rhs.eta0)

    copyto!(u_v_eta.u, rhs.u0)
    copyto!(u_v_eta.v, rhs.v0)
    copyto!(u_v_eta.eta, rhs.eta0)

    return nothing 

end 

function comp_u_v_eta_t_check(rhs, params, nx, interp, grad, advec) 

    rhs.h .= rhs.eta1 .+ params.H 

    rhs.h_u .= interp.ITu * rhs.h
    rhs.h_v .= interp.ITv * rhs.h
    rhs.h_q .= interp.ITq * rhs.h

    rhs.U .= rhs.u1 .* rhs.h_u 
    rhs.V .= rhs.v1 .* rhs.h_v 

    rhs.IuT_u1 .= interp.IuT * (rhs.u1.^2)
    rhs.IvT_v1 .= interp.IvT * (rhs.v1.^2)
    rhs.kinetic .= rhs.IuT_u1 + rhs.IvT_v1

    # Kloewer defined new terms q and p corresponding to potential vorticity and 
    # Bernoulli potential respectively. To avoid errors in my mimic I'm following 
    # along and doing the same 
    rhs.q .= (params.coriolis + grad.Gvx * rhs.v1 - grad.Guy * rhs.u1) ./ rhs.h_q 
    rhs.p .= 0.5 .* rhs.kinetic .+ params.g .* rhs.h

    # bottom friction
    rhs.kinetic_sq .= (rhs.kinetic).^(1/2)
    rhs.bfric_u .= params.bottom_drag .* ((interp.ITu * rhs.kinetic_sq) .* rhs.u1) ./ rhs.h_u
    rhs.bfric_v .= params.bottom_drag .* ((interp.ITv * rhs.kinetic_sq) .* rhs.v1) ./ rhs.h_v

    # deal with the advection term 
    comp_advection_check(rhs, nx, advec)

    rhs.Mu .= params.A_h .* (grad.LLu * rhs.u1)
    rhs.Mv .= params.A_h .* (grad.LLv * rhs.v1) 

    rhs.u_t .= rhs.adv_u - grad.GTx * rhs.p + params.wind_stress ./ rhs.h_u - rhs.Mu - rhs.bfric_u

    rhs.v_t .= rhs.adv_v - grad.GTy * rhs.p - rhs.Mv - rhs.bfric_v 

    rhs.eta_t .= - (grad.Gux * rhs.U + grad.Gvy * rhs.V)

    return nothing

end 

function comp_advection_check(rhs, nx, advec)

    rhs.AL1q .= advec.AL1 * rhs.q 
    rhs.AL2q .= advec.AL2 * rhs.q 

    AL1q_au = @view rhs.AL1q[1:end-nx]
    AL2q_bu = @view rhs.AL2q[nx+1:end]
    AL2q_cu = @view rhs.AL2q[1:end-nx]
    AL1q_du = @view rhs.AL1q[nx+1:end]

    AL1q_av = @view rhs.AL1q[advec.index_av]
    AL2q_bv = @view rhs.AL2q[advec.index_bv]
    AL2q_cv = @view rhs.AL2q[advec.index_cv]
    AL1q_dv = @view rhs.AL1q[advec.index_dv]

    # rhs.AL1q_av .= rhs.AL1q[advec.index_av]
    # rhs.AL2q_bv .= rhs.AL2q[advec.index_bv]
    # rhs.AL2q_cv .= rhs.AL2q[advec.index_cv]
    # rhs.AL1q_dv .= rhs.AL1q[advec.index_dv]

    rhs.adv_u .= advec.Seur * (advec.ALeur * rhs.q .* rhs.U) + advec.Seul * (advec.ALeul * rhs.q .* rhs.U) +
    advec.Sau * (AL1q_au .* rhs.V) + advec.Sbu * (AL2q_bu .* rhs.V) + 
    advec.Scu * (AL2q_cu .* rhs.V) + advec.Sdu * (AL1q_du .* rhs.V)

    rhs.adv_v .= advec.Spvu * ((advec.ALpvu * rhs.q) .* rhs.V) + advec.Spvd * ((advec.ALpvd * rhs.q) .* rhs.V) -
    advec.Sav * (AL1q_av .* rhs.U) - advec.Sbv * (AL2q_bv .* rhs.U) - 
    advec.Scv * (AL2q_cv .* rhs.U) - advec.Sdv * (AL1q_dv .* rhs.U)

    # rhs.adv_u .= advec.Seur * (advec.ALeur * rhs.q .* rhs.U) + advec.Seul * (advec.ALeul * rhs.q .* rhs.U) +
    # advec.Sau * (rhs.AL1q_au .* rhs.V) + advec.Sbu * (rhs.AL2q_bu .* rhs.V) + 
    # advec.Scu * (rhs.AL2q_cu .* rhs.V) + advec.Sdu * (rhs.AL1q_du .* rhs.V)

    # rhs.adv_v .= advec.Spvu * ((advec.ALpvu * rhs.q) .* rhs.V) + advec.Spvd * ((advec.ALpvd * rhs.q) .* rhs.V) -
    # advec.Sav * (rhs.AL1q_av .* rhs.U) - advec.Sbv * (rhs.AL2q_bv .* rhs.U) - 
    # advec.Scv * (rhs.AL2q_cv .* rhs.U) - advec.Sdv * (rhs.AL1q_dv .* rhs.U)

    return nothing

end

function debug_main()

Lx = 3840e3                    # E-W length of the domain [meters]
Ly = 3840e3                    # N-S length of the domain [meters]
nx = 100                        # number of cells in the x-direction
ny = 100                       # number of cells in the y-direction

# based on above values this returns more parameters related to the four grids 
grid_params = build_grid(Lx, Ly, nx, ny)

gyre_params = def_params(grid_params)

# building discrete operators 
grad_ops = build_derivs(grid_params)                # discrete gradient operators 
interp_ops = build_interp(grid_params, grad_ops)    # discrete interpolation operators (travels between grids)
advec_ops = build_advec(grid_params)
rhs_terms = RHS_terms(Nu = grid_params.Nu, Nv = grid_params.Nv, NT = grid_params.NT, Nq = grid_params.Nq)
# rhs_terms = allocate(grid_params, grad_ops, interp_ops, advec_ops)

# starting from rest ---> all initial conditions are zero 

# how long to spinup the model for 
Tspinup_days = 90 # [days] 

# how long to run the model for after spinup
Trun_days = 1 * 365     # [days] 

Tspinup, Trun = days_to_seconds(Tspinup_days, Trun_days, gyre_params.dt)

uout = zeros(grid_params.Nu) 
vout = zeros(grid_params.Nv) 
etaout = zeros(grid_params.NT) 

u_v_eta = gyre_vector(uout, vout, etaout)

@time for t in 1:Tspinup
    advance_check(u_v_eta, rhs_terms, gyre_params, grid_params, interp_ops, grad_ops, advec_ops) 
end

end 