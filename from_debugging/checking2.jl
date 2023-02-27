function advance_check2(u_v_eta, rhs, params, grid, interp, grad, advec) 

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

        comp_u_v_eta_t_check2(rhs, params, grid, interp, grad, advec)

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

function comp_u_v_eta_t_check2(rhs, params, grid, interp, grad, advec) 

    rhs.h .= rhs.eta1 .+ params.H 

    rhs.h_u .= interp.ITu * rhs.h
    rhs.h_v .= interp.ITv * rhs.h
    rhs.h_q .= interp.ITq * rhs.h

    rhs.U .= rhs.u1 .* rhs.h_u 
    rhs.V .= rhs.v1 .* rhs.h_v 

    rhs.kinetic .= interp.IuT * (rhs.u1.^2) + interp.IvT * (rhs.v1.^2)

    # Kloewer defined new terms q and p corresponding to potential vorticity and 
    # Bernoulli potential respectively. To avoid errors in my mimic I'm following 
    # along and doing the same 
    rhs.q .= (params.coriolis + grad.Gvx * rhs.v1 - grad.Guy * rhs.u1) ./ rhs.h_q 
    rhs.p .= 0.5 .* rhs.kinetic + params.g * rhs.h

    # bottom friction
    rhs.kinetic_sq .= (rhs.kinetic).^(1/2)
    rhs.bfric_u .= params.bottom_drag * ((interp.ITu * rhs.kinetic_sq) .* rhs.u1) ./ rhs.h_u
    rhs.bfric_v .= params.bottom_drag * ((interp.ITv * rhs.kinetic_sq) .* rhs.v1) ./ rhs.h_v

    # deal with the advection term 
    comp_advection_check2(rhs, grid.nx, advec)

    rhs.Mu .= params.A_h .* (grad.LLu * rhs.u1)
    rhs.Mv .= params.A_h .* (grad.LLv * rhs.v1) 

    rhs.u_t .= rhs.adv_u - grad.GTx * rhs.p + params.wind_stress ./ rhs.h_u - rhs.Mu - rhs.bfric_u

    rhs.v_t .= rhs.adv_v - grad.GTy * rhs.p - rhs.Mv - rhs.bfric_v 

    rhs.eta_t .= - (grad.Gux * rhs.U + grad.Gvy * rhs.V)

    return nothing

end 

function comp_advection_check2(rhs, nx, advec)

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

    rhs.adv_u .= advec.Seur * ((advec.ALeur * rhs.q) .* rhs.U) + advec.Seul * ((advec.ALeul * rhs.q) .* rhs.U) + 
    advec.Sau * (AL1q_au .* rhs.V) + advec.Sbu * (AL2q_bu .* rhs.V) + 
    advec.Scu * (AL2q_cu .* rhs.V) + advec.Sdu * (AL1q_du .* rhs.V)

    rhs.adv_v .= advec.Spvu * ((advec.ALpvu * rhs.q) * rhs.V) + advec.Spvd * ((advec.ALpvd * rhs.q) * rhs.V) - 
    advec.Sav * (AL1q_av .* rhs.U) - advec.Sbv * (AL2q_bv .* rhs.U) - 
    advec.Scv * (AL2q_cv .* rhs.U) - advec.Sdv * (AL1q_dv .* rhs.U)

    return nothing

end