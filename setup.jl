function setup(Lx, Ly, nx, ny)

    # build grid structure 
    grid_params = build_grid(Lx, Ly, nx, ny)

    # initialize parameters 
    gyre_params = def_params(grid_params)

    # building discrete operators 
    grad_ops = build_derivs(grid_params)                # discrete gradient operators 
    interp_ops = build_interp(grid_params, grad_ops)    # discrete interpolation operators (travels between grids)
    advec_ops = build_advec(grid_params)

    # allocating RHS terms 
    rhs_terms = RHS_terms(Nu = grid_params.Nu, Nv = grid_params.Nv, NT = grid_params.NT, Nq = grid_params.Nq)

    return grid_params, gyre_params, grad_ops, interp_ops, advec_ops, rhs_terms 

end