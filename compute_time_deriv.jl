function comp_u_v_eta_t(u, v, eta, params, interp_ops, grad_ops, advec_ops) 

    # @assert all(x -> x < 10, u)
    # if any(x -> x > 5, v)
    #     print(u)
    #     print(v)
    #     print(eta)
    # end
    # @assert all(x -> x < 1, v)
    # @assert all(x -> x < 5, eta)

    h = eta .+ params.H 

    @assert all(x -> x < 550, h)

    h_u = interp_ops.ITu * h
    h_v = interp_ops.ITv * h
    h_q = interp_ops.ITq * h

    @assert all(x -> x < 1e4, h_u)
    @assert all(x -> x < 1e4, h_v)
    @assert all(x -> x < 1e4, h_q)

    U = u .* h_u 
    V = v .* h_v 

    @assert all(x -> x < 1e4, U)
    @assert all(x -> x < 1e4, V)

    kinetic = interp_ops.IuT * (u.^2) + interp_ops.IvT * (v.^2)

    @assert all(x -> x < 1e4, u.^2)
    @assert all(x -> x < 1e4, v.^2)
    @assert all(x -> x < 1e5, kinetic)

    # Kloewer defined new terms q and p corresponding to potential vorticity and 
    # Bernoulli potential respectively. To avoid errors in my mimic I'm following 
    # along and doing the same 
    q = (params.coriolis + grad_ops.Gvx * v - grad_ops.Guy * u) ./ h_q 
    p = 0.5 .* kinetic .+ params.g .* h

    @assert all(x -> x < 1e4, q)
    @assert all(x -> x < 1e4, p)

    # bottom friction
    kinetic_sq = (kinetic).^(1/2)
    bfric_u = params.bottom_drag .* ((interp_ops.ITu * kinetic_sq) .* u) ./ h_u
    bfric_v = params.bottom_drag .* ((interp_ops.ITv * kinetic_sq) .* v) ./ h_v

    @assert all(x -> x < 1e4, bfric_u)
    @assert all(x -> x < 1e4, bfric_v)

    # deal with the advection term 
    adv_u, adv_v = comp_advection(q, U, V, advec_ops)

    # @assert all(isfinite, adv_u)
    # @assert all(isfinite, adv_v)

    @assert all(x -> x < 1e4, adv_u)
    @assert all(x -> x < 1e4, adv_v)

    Mu = params.A_h .* (grad_ops.LLu * u)
    Mv = params.A_h .* (grad_ops.LLv * v) 

    # @assert all(isfinite, params.A_h)
    # @assert all(isfinite, grad_ops.LLu)
    # @assert all(isfinite, grad_ops.LLv)
    # @assert all(isfinite, Mu)
    # @assert all(isfinite, Mv)

    @assert all(x -> x < 1e4, Mu)
    @assert all(x -> x < 1e4, Mv)

    rhs_u = adv_u - grad_ops.GTx * p + params.wind_stress ./ h_u - Mu - bfric_u

    # @assert all(isfinite, grad_ops.GTx)
    # @assert all(isfinite, rhs_u)
    @assert all(x -> x < 1e4, rhs_u)

    rhs_v = adv_v - grad_ops.GTy * p - Mv - bfric_v 

    # @assert all(isfinite, grad_ops.GTy)
    # @assert all(isfinite, rhs_v)

    @assert all(x -> x < 1e4, rhs_v)

    rhs_eta = - (grad_ops.Gux * U + grad_ops.Gvy * V)

    # @assert all(isfinite, grad_ops.Gux)
    # @assert all(isfinite, grad_ops.Gvy)
    # @assert all(isfinite, rhs_eta)

    @assert all(x -> x < 1e4, rhs_eta)

    return rhs_u, rhs_v, rhs_eta

end 

function comp_advection(q, U, V, advec_ops)

    # @assert all(isfinite, q)
    # @assert all(isfinite, U)
    # @assert all(isfinite, V)

    @assert all(x -> x < 1e4, q)
    @assert all(x -> x < 1e4, U)
    @assert all(x -> x < 1e4, V)

    AL1q = advec_ops.AL1 * q 
    AL2q = advec_ops.AL2 * q 

    # @assert all(isfinite, AL1q)
    # @assert all(isfinite, AL2q)

    # adv_u = advec_ops.Seur * ((advec_ops.ALeur * q) .* U) + advec_ops.Seul * ((advec_ops.ALeul * q) .* U) +
    # advec_ops.Sau * (AL1q[1:end-nx] .* V) + advec_ops.Sbu * (AL2q[nx+1:end] .* V) + 
    # advec_ops.Scu * (AL2q[1:end-nx] .* V) + advec_ops.Sdu * (AL1q[nx+1:end] .* V)
    
    # adv_v = advec_ops.Spvu * ((advec_ops.ALpvu * q) .* U) + advec_ops.Spvd * ((advec_ops.ALpvd * q) .* V) -
    # advec_ops.Sav * (AL1q[advec_ops.index_av] .* U) - advec_ops.Sbv * (AL2q[advec_ops.index_bv] .* U) - 
    # advec_ops.Scv * (AL2q[advec_ops.index_cv] .* U) - advec_ops.Sdv * (AL1q[advec_ops.index_dv] .* U)

    adv_u = advec_ops.Seur * (advec_ops.ALeur * q .* U) + advec_ops.Seul * (advec_ops.ALeul * q .* U) +
    advec_ops.Sau * (AL1q[1:end-nx] .* V) + advec_ops.Sbu * (AL2q[nx+1:end] .* V) + 
    advec_ops.Scu * (AL2q[1:end-nx] .* V) + advec_ops.Sdu * (AL1q[nx+1:end] .* V)

    # @assert all(isfinite, adv_u)
    @assert all(x -> x < 1e4, adv_u)
    
    adv_v = advec_ops.Spvu * ((advec_ops.ALpvu * q) .* V) + advec_ops.Spvd * ((advec_ops.ALpvd * q) .* V) -
    advec_ops.Sav * (AL1q[advec_ops.index_av] .* U) - advec_ops.Sbv * (AL2q[advec_ops.index_bv] .* U) - 
    advec_ops.Scv * (AL2q[advec_ops.index_cv] .* U) - advec_ops.Sdv * (AL1q[advec_ops.index_dv] .* U)

    # @assert all(isfinite, adv_v)
    @assert all(x -> x < 1e4, adv_v)

    return adv_u, adv_v

end