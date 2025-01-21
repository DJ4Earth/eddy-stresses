function cost_eval(param_guess)

    energy_high_resolution = load_object("./spinup_files/1024_postspinup_noslip_5years_061824/energy_post_spinup_1024_noslip_5years_061224.jld2")
    grid_scale = 8

    # aiming to have data the final 7 days
    data_steps = 23*225:100:30*225
    data = energy_high_resolution[data_steps[1]*grid_scale:225*grid_scale:data_steps[end]*grid_scale]

    S = ShallowWaters.model_setup(
    output=false,
    L_ratio=1,
    g=9.81,
    H=500,
    wind_forcing_x="double_gyre",
    Lx=3840e3,
    seasonal_wind_x=false,
    topography="flat",
    bc="nonperiodic",
    bottom_drag="quadratic",
    nn_forcing_dissipation=true,
    α=2,
    # νB=1000,
    nx=128,
    Ndays=30,
    data=data,
    data_steps=data_steps,
    param_guess[1] = weights_diagonal,
    param_guess[2] = weights_offdiagonal
    # initial_cond="ncfile",
    # initpath="./run_0001/"
    )

    # snaps = Int(floor(sqrt(S.grid.nt)))
    # revolve = Revolve{ShallowWaters.ModelSetup}(S.grid.nt,
    # snaps;
    # verbose=1,
    # gc=true,
    # write_checkpoints=false
    # )

    J = integration(S)

    return J

end

function gradient_eval(G, param_guess)

    energy_high_resolution = load_object("./spinup_files/1024_postspinup_noslip_5years_061824/energy_post_spinup_1024_noslip_5years_061224.jld2")
    grid_scale = 8

    # aiming to have data the final 7 days
    data_steps = 23*225:100:30*225
    data = energy_high_resolution[data_steps[1]*grid_scale:225*grid_scale:data_steps[end]*grid_scale]

    S = ShallowWaters.model_setup(
    output=false,
    L_ratio=1,
    g=9.81,
    H=500,
    wind_forcing_x="double_gyre",
    Lx=3840e3,
    seasonal_wind_x=false,
    topography="flat",
    bc="nonperiodic",
    bottom_drag="quadratic",
    nn_forcing_dissipation=true,
    α=2,
    # νB=1000,
    nx=128,
    Ndays=30,
    data=data,
    data_steps=data_steps,
    param_guess[1] = weights_diagonal,
    param_guess[2] = weights_offdiagonal
    # initial_cond="ncfile",
    # initpath="./run_0001/"
    )

    dS = Enzyme.Compiler.make_zero(Core.Typeof(S), IdDict(), S)
    # snaps = Int(floor(sqrt(S.grid.nt)))
    # revolve = Revolve{ShallowWaters.ModelSetup}(S.grid.nt,
    # snaps;
    # verbose=1,
    # gc=true,
    # write_checkpoints=false
    # )

    autodiff(Enzyme.ReverseWithPrimal,
    integration,
    Duplicated(S, dS)
    )

    G[1] = dS.parameters.γ₀

    return nothing

end

function FG(F, G, param_guess)

    G === nothing || gradient_eval(G, param_guess)
    F === nothing || return cost_eval(param_guess)

end
