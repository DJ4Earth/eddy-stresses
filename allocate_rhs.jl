# This function will allocate all the space needed for terms that appear when 
# we compute the RHS of the system 

function allocate(grid, grad, interp, advec)

    u0 = zeros(grid.Nu)
    v0 = zeros(grid.Nv)
    eta0 = zeros(grid.NT)
    u1 = zeros(grid.Nu)
    v1 = zeros(grid.Nv)
    eta1 = zeros(grid.NT)

    h = zeros(grid.NT)
    h_u = zeros(size(interp.ITu)[1])
    h_v = zeros(size(interp.ITv)[1])
    h_q = zeros(size(interp.ITq)[1])
    U = zeros(grid.Nu)
    V = zeros(grid.Nv)
    kinetic = zeros(grid.NT)
    kinetic_sq = zeros(grid.NT)
    q = zeros(size(interp.ITq)[1])
    p = zeros(grid.NT)
    bfric_u = zeros(size(interp.ITu)[1])
    bfric_v = zeros(size(interp.ITv)[1])
    Mu = zeros(size(grad.LLu)[1])
    Mv = zeros(size(grad.LLv)[1])
    rhs_u = zeros(grid.Nu)
    rhs_v = zeros(grid.Nv)
    rhs_eta = zeros(grid.NT)

    AL1q = zeros(size(advec.AL1)[1])
    AL2q = zeros(size(advec.AL2)[1])
    AL1q_au = zeros(grid.Nv)
    AL1q_du = zeros(grid.Nv)
    AL1q_av = zeros(grid.Nu)
    AL1q_dv = zeros(grid.Nu)
    AL2q_bu = zeros(grid.Nv)
    AL2q_cu = zeros(grid.Nv)
    AL2q_bv = zeros(grid.Nu)
    AL2q_cv = zeros(grid.Nu)

    rhs = RHS_terms(
    u0,
    v0,
    eta0,
    u1,
    v1,
    eta1,
    h,          
    h_u,        
    h_v,        
    h_q,
    U,
    V,
    kinetic,
    kinetic_sq,
    q,
    p,
    bfric_u,
    bfric_v,
    Mu,
    Mv,
    rhs_u,
    rhs_v,
    rhs_eta,
    AL1q,
    AL2q,
    AL1q_au,
    AL1q_du,
    AL1q_av,
    AL1q_dv,
    AL2q_bu,
    AL2q_cu,
    AL2q_bv,
    AL2q_cv,
    )

    return rhs

end