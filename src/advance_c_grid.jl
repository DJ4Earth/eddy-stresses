# Contains one function advance, which just takes a single step forward in time. 
# This will differ from advance.jl: here I'm going to try and implement a C-grid scheme 

function advance(u_v_eta, grid, rhs, params, interp, grad, advec) 

    nx = grid.nx 
    dt = params.dt

    # we now use RK4 as the timestepper, here I'm storing the coefficients needed for this 
    rk_a = [1/6, 1/3, 1/3, 1/6]
    rk_b = [1/2, 1/2, 1.]

    rhs.umid .= u_v_eta.u
    rhs.vmid .= u_v_eta.v
    rhs.etamid .= u_v_eta.eta

    rhs.u0 .= u_v_eta.u
    rhs.v0 .= u_v_eta.v
    rhs.eta0 .= u_v_eta.eta

    rhs.u1 .= u_v_eta.u
    rhs.v1 .= u_v_eta.v
    rhs.eta1 .= u_v_eta.eta

    for j in 1:4

        comp_u_v_eta_t(nx, rhs, params, interp, grad, advec)

        if j < 4
            rhs.u1 .= rhs.umid .+ rk_b[j] .* dt .* rhs.u_t
            rhs.v1 .= rhs.vmid .+ rk_b[j] .* dt .* rhs.v_t
            rhs.eta1 .= rhs.etamid .+ rk_b[j] .* dt .* rhs.eta_t
        end

        rhs.u0 .= rhs.u0 .+ rk_a[j] .* dt .* rhs.u_t
        rhs.v0 .= rhs.v0 .+ rk_a[j] .* dt .* rhs.v_t 
        rhs.eta0 .= rhs.eta0 .+ rk_a[j] .* dt .* rhs.eta_t 

    end

    @assert all(x -> x < 10.0, rhs.u0)
    @assert all(x -> x < 10.0, rhs.v0)
    @assert all(x -> x < 10.0, rhs.eta0)

    copyto!(u_v_eta.u, rhs.u0)
    copyto!(u_v_eta.v, rhs.v0)
    copyto!(u_v_eta.eta, rhs.eta0)

    return nothing 

end 