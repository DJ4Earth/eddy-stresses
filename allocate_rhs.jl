# This function will allocate all the space needed for terms that appear when 
# we compute the RHS of the system 

function allocate(grid, grad, interp, advec)

    # u0 = zeros(grid.Nu)
    # v0 = zeros(grid.Nv)
    # eta0 = zeros(grid.NT)
    # u1 = zeros(grid.Nu)
    # v1 = zeros(grid.Nv)
    # eta1 = zeros(grid.NT)

    # h = zeros(grid.NT)
    # h_u = zeros(size(interp.ITu)[1])
    # h_v = zeros(size(interp.ITv)[1])
    # h_q = zeros(size(interp.ITq)[1])
    # U = zeros(grid.Nu)
    # V = zeros(grid.Nv)
    # kinetic = zeros(grid.NT)
    # kinetic_sq = zeros(grid.NT)
    # q = zeros(size(interp.ITq)[1])
    # p = zeros(grid.NT)
    # bfric_u = zeros(size(interp.ITu)[1])
    # bfric_v = zeros(size(interp.ITv)[1])
    # Mu = zeros(size(grad.LLu)[1])
    # Mv = zeros(size(grad.LLv)[1])
    # u_t = zeros(grid.Nu)
    # v_t = zeros(grid.Nv)
    # eta_t = zeros(grid.NT)
    # adv_u = zeros(grid.Nu)
    # adv_v = zeros(grid.Nv)

    # AL1q = zeros(size(advec.AL1)[1])
    # AL2q = zeros(size(advec.AL2)[1])
    # AL1q_au = zeros(grid.Nv)
    # AL1q_du = zeros(grid.Nv)
    # AL1q_av = zeros(grid.Nu)
    # AL1q_dv = zeros(grid.Nu)
    # AL2q_bu = zeros(grid.Nv)
    # AL2q_cu = zeros(grid.Nv)
    # AL2q_bv = zeros(grid.Nu)
    # AL2q_cv = zeros(grid.Nu)

    # rhs = RHS_terms(
    # u0,
    # v0,
    # eta0,
    # u1,
    # v1,
    # eta1,
    # h,          
    # h_u,        
    # h_v,        
    # h_q,
    # U,
    # V,
    # kinetic,
    # kinetic_sq,
    # q,
    # p,
    # bfric_u,
    # bfric_v,
    # Mu,
    # Mv,
    # u_t,
    # v_t,
    # eta_t,
    # adv_u,
    # adv_v,
    # AL1q,
    # AL2q,
    # AL1q_au,
    # AL1q_du,
    # AL1q_av,
    # AL1q_dv,
    # AL2q_bu,
    # AL2q_cu,
    # AL2q_bv,
    # AL2q_cv
    # )

    rhs = RHS_terms(
    u0 = zeros(grid.Nu),
    v0 = zeros(grid.Nv),
    eta0 = zeros(grid.NT),
    u1 = zeros(grid.Nu),
    v1 = zeros(grid.Nv),
    eta1 = zeros(grid.NT),
    h = zeros(grid.NT),
    h_u = zeros(grid.Nu),
    h_v = zeros(grid.Nv),
    h_q = zeros(grid.Nq),
    U = zeros(grid.Nu),
    V = zeros(grid.Nv),
    IuT_u = zeros(grid.NT),
    IvT_v = zeros(grid.NT),
    kinetic = zeros(grid.NT),
    kinetic_sq = zeros(grid.NT),
    Gvx_v1 = zeros(grid.Nq),
    Guy_u1 = zeros(grid.Nq),
    q = zeros(grid.Nq),
    p = zeros(grid.NT),
    ITu_ksq = zeros(grid.Nu),
    ITv_ksq = zeros(grid.Nv),
    bfric_u = zeros(grid.Nu),
    bfric_v = zeros(grid.Nv),
    LLu_u1 = zeros(grid.Nu),
    LLv_v1 = zeros(grid.Nv),
    Mu = zeros(grid.Nu),
    Mv = zeros(grid.Nv),
    GTx_p = zeros(grid.Nu),
    u_t = zeros(grid.Nu),
    GTy_p = zeros(grid.Nv),
    v_t = zeros(grid.Nv),
    Gux_U = zeros(grid.NT),
    Gvy_V = zeros(grid.NT),
    eta_t = zeros(grid.NT),
    adv_u = zeros(grid.Nu),
    adv_v = zeros(grid.Nv),
    AL1q = zeros(size(advec.AL1)[1]),
    AL2q = zeros(size(advec.AL2)[1])
    # ALeur_q = zeros(grid.Nu),
    # ALeul_q = zeros(grid.Nu)
    )

    return rhs

end