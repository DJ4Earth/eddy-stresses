# minimizing the error I'm getting from Enzyme
using Parameters
using Enzyme_jll#main

mutable struct gyre_vector
    u::Vector{Float64}
    v::Vector{Float64}
    eta::Vector{Float64}
end

@with_kw struct RHS_terms_debug
    Nu::Int
    u0::Vector{Float64} = zeros(Nu) 
end

function advance(u_v_eta, rhs)

    rk_b = [1/2, 1/2, 1]
    rhs.u0 .= u_v_eta.u
    rhs.u0 .= rhs.u0 

    copyto!(u_v_eta.u, rhs.u0)

    return nothing 

end 

nx = 5             # grid resolution in x-direction
ny = 5             # grid resolution in y-direction

NT = nx * ny 
Nu = (nx - 1) * ny 
Nv = (ny - 1) * nx 
Nq = (nx + 1) * (ny + 1)

rhs_terms = RHS_terms_debug(Nu = Nu)
ad_rhs_terms = RHS_terms_debug(Nu = Nu)

u_v_eta = gyre_vector(zeros(Nu), zeros(Nv), zeros(NT))
ad_u_v_eta = gyre_vector(zeros(Nu), zeros(Nv), zeros(NT))


# checking that the function does in fact run
advance(u_v_eta, rhs_terms)

autodiff(advance,
    Duplicated(u_v_eta, ad_u_v_eta),
    Duplicated(rhs_terms, ad_rhs_terms)
)