function setup(nx, ny)

    Lx = 3840e3                     # E-W length of the domain [meters]
    Ly = 3840e3                     # N-S length of the domain [meters]

    grid_params = build_grid(Lx, Ly, nx, ny)
    gyre_params = def_params(grid_params)
    
    # building discrete operators 
    grad_ops = build_derivs(grid_params)                # discrete gradient operators 
    interp_ops = build_interp(grid_params, grad_ops)    # discrete interpolation operators (travels between grids)
    advec_ops = build_advec(grid_params)
    rhs_terms = RHS_terms(Nu = grid_params.Nu, Nv = grid_params.Nv, NT = grid_params.NT, Nq = grid_params.Nq)

    return grid_params, gyre_params, grad_ops, interp_ops, advec_ops, rhs_terms

end