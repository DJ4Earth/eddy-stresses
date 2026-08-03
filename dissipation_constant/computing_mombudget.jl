# offline check to make sure I'm computing everything correctly. The momentum budget terms should
# sum up to the euler tendencies I compute, du and dv
# check passed
function momentum_budget_offline()
    Plr = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        zb_forcing_momentum=false,
        zb_forcing_dissipation=false,
        zb_filtered=true,
        # nn_forcing_momentum=false,
        # nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=2,
        initial_cond="rest",
        # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
    );
    Slr = ShallowWaters.model_setup(Plr);


    Plr2 = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=1e-20,
        H=500,
        cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        ω=1e-22,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        zb_forcing_momentum=false,
        zb_forcing_dissipation=false,
        zb_filtered=true,
        # nn_forcing_momentum=false,
        # nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=2,
        initial_cond="rest",
        # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
    );
    Slr2 = ShallowWaters.model_setup(Plr2);

    Seuler = deepcopy(Slr);
    Sadv = deepcopy(Slr2);

    t = 225 * Slr.grid.dtint
    n = 1096

    u_, v_, eta_ = ShallowWaters.add_halo(uhrcgall[:,:,n], vhrcgall[:,:,n], etahrcgall[:,:,n], Slr)

    Seuler.Prog.u = u_
    Seuler.Prog.v = v_
    Seuler.Prog.η = eta_

    Sadv.Prog.u = u_
    Sadv.Prog.v = v_
    Sadv.Prog.η = eta_

    du, dv, adv_u, adv_v, fu, fv, detagdx, detagdy, bottomdragu, bottomdragv, viscu, viscv, Fx = compute_mom_budget(Seuler, Sadv, t, n);

    detagdx2 = zeros(127,128)
    detagdy2 = zeros(128,127)

    ShallowWaters.∂x!(detagdx2, (Seuler.constants.g .* etahrcgall[:,:,n])./ Seuler.grid.Δ)
    ShallowWaters.∂y!(detagdy2, (Seuler.constants.g .* etahrcgall[:,:,n])./ Seuler.grid.Δ)

    fig = Figure(fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ((2 .* du) ./ 384) .- (adv_u + fv - detagdx + Fx + viscu + bottomdragu),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Euler tendency"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    fig = Figure(size=(1050, 700), fontsize=15);

    Label(
        fig[0, 3],
        "Offline u budget",
        fontsize = 20,
        tellwidth = false
    )

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ((2 .* du) ./ 384),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Euler tendency"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (adv_u + fv - detagdx + Fx + viscu + bottomdragu),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Sum of budget terms"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax2)
    Colorbar(fig[1,4], hm2)

    ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (adv_u),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Advection"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax3)
    Colorbar(fig[1,6], hm3)

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    fv,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Coriolis"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    Colorbar(fig[2,2], hm4)
    hidexdecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -detagdx,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Gravity"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    hidedecorations!(ax5)
    Colorbar(fig[2,4], hm5)

    ax6, hm6 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Fx,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Wind stress"),
    # colorrange=(-1.5e-5, 1.5e-5)
    );
    hideydecorations!(ax6)
    Colorbar(fig[2,6], hm6)

    ax7, hm7 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    viscu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Viscosity"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    Colorbar(fig[3,2], hm7)

    ax8, hm8 = heatmap(fig[3,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    bottomdragu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Bottom drag"),
    colorrange=(-1.5e-8, 1.5e-8)
    );
    hideydecorations!(ax8)
    Colorbar(fig[3,4], hm8)

end

function momentum_budget_multi2()

    Plr = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        zb_forcing_momentum=false,
        zb_forcing_dissipation=false,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=2,
        initial_cond="rest",
        # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
    );
    Slr = ShallowWaters.model_setup(Plr);


    Plr2 = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=1e-12,
        H=500,
        cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        ω=1e-16,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        zb_forcing_momentum=false,
        zb_forcing_dissipation=false,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=2,
        initial_cond="rest",
        # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
    );
    Slr2 = ShallowWaters.model_setup(Plr2);

    Seuler = deepcopy(Slr);
    current = 1
    for m in (Seuler.Diag.CNNVars.model_Su, Seuler.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-6-8-10-13-15-18-23-28-33-38-41-44-48-51-53-58-63-65-68-73-78-83-86-88daystart_2dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    Sadv = deepcopy(Slr2);
    current = 1
    for m in (Sadv.Diag.CNNVars.model_Su, Sadv.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-6-8-10-13-15-18-23-28-33-38-41-44-48-51-53-58-63-65-68-73-78-83-86-88daystart_2dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    t = 225 * Slr.grid.dtint
    n = 1096

    u_, v_, eta_ = ShallowWaters.add_halo(umulti2[:,:,n], vmulti2[:,:,n], etamulti2[:,:,n], Slr)

    Seuler.Prog.u = u_
    Seuler.Prog.v = v_
    Seuler.Prog.η = eta_

    Sadv.Prog.u = u_
    Sadv.Prog.v = v_
    Sadv.Prog.η = eta_

    du, dv, adv_u, adv_v, fu, fv, detagdx, detagdy, bottomdragu, bottomdragv, viscu, viscv, Fx = compute_mom_budget(Seuler, Sadv, t, n);

    # detagdx2 = zeros(127,128)
    # detagdy2 = zeros(128,127)

    # ShallowWaters.∂x!(detagdx2, (Seuler.constants.g .* etahrcgall[:,:,n])./ Seuler.grid.Δ)
    # ShallowWaters.∂y!(detagdy2, (Seuler.constants.g .* etahrcgall[:,:,n])./ Seuler.grid.Δ)

    # all of the rhs terms are computed with a *scaled* prognostic variable, so to actually 
    # see the momentum budget terms we need to remove that scaling. The division by Delta is because
    # all of the timestepping is done with dt / Delta, and we want the term that gets multiplied by dt
    # I should have corrected for the scalings in the values returned by compute_mom_budget already,
    # but for the stuff I'm pulling after (see below) I need to do it manually

    utendency = ((2 .* du) ./ Seuler.grid.dtint)
    parameterization = Seuler.Diag.CNNVars.S_u ./ (Seuler.grid.Δ * Seuler.constants.scale);
    budgetsum = (adv_u + fv - detagdx + Fx + viscu + bottomdragu + parameterization);

    fig = Figure(fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    utendency .- budgetsum,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Euler tendency - budget sum"),
    # colorrange=(-5e-7, 5e-7)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    fig = Figure(size=(1050, 700), fontsize=15);

    Label(
        fig[0, 3],
        "Online u budget, ensemble 2 day parameterization",
        fontsize = 25,
        tellwidth = false
    )

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    utendency,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Euler tendency"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    budgetsum,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Sum of budget terms"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax2)
    Colorbar(fig[1,4], hm2)

    ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (adv_u),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Nonlinear advection"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax3)
    Colorbar(fig[1,6], hm3)

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    fv,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Coriolis"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    Colorbar(fig[2,2], hm4)
    hidexdecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -detagdx,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Pressure gradient"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    hidedecorations!(ax5)
    Colorbar(fig[2,4], hm5)

    ax6, hm6 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Fx,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Wind stress"),
    # colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax6)
    Colorbar(fig[2,6], hm6)

    ax7, hm7 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    viscu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Viscosity"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    Colorbar(fig[3,2], hm7)

    ax8, hm8 = heatmap(fig[3,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    bottomdragu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Bottom drag"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    hideydecorations!(ax8)
    Colorbar(fig[3,4], hm8)

    ax9, hm9 = heatmap(fig[3,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    parameterization,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Parameterization"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hideydecorations!(ax9)
    Colorbar(fig[3,6], hm9)

end

function momentum_budget_zannabolton()
    Plr = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        zb_forcing_momentum=false,
        zb_forcing_dissipation=true,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=false,
        N=1,
        α=2,
        nx=128,
        Ndays=2,
        initial_cond="rest",
        # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
    );
    Slr = ShallowWaters.model_setup(Plr);


    Plr2 = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=1e-12,
        H=500,
        cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        ω=1e-16,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        zb_forcing_momentum=false,
        zb_forcing_dissipation=true,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=false,
        N=1,
        α=2,
        nx=128,
        Ndays=2,
        initial_cond="rest",
        # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
    );
    Slr2 = ShallowWaters.model_setup(Plr2);

    Seuler = deepcopy(Slr);
    Sadv = deepcopy(Slr2);

    t = 225 * Slr.grid.dtint
    n = 1096

    u_, v_, eta_ = ShallowWaters.add_halo(uzb[:,:,n], vzb[:,:,n], etazb[:,:,n], Slr)

    Seuler.Prog.u = u_
    Seuler.Prog.v = v_
    Seuler.Prog.η = eta_

    Sadv.Prog.u = u_
    Sadv.Prog.v = v_
    Sadv.Prog.η = eta_

    du, dv, adv_u, adv_v, fu, fv, detagdx, detagdy, bottomdragu, bottomdragv, viscu, viscv, Fx = compute_mom_budget(Seuler, Sadv, t, n);

    # detagdx2 = zeros(127,128)
    # detagdy2 = zeros(128,127)

    # ShallowWaters.∂x!(detagdx2, (Seuler.constants.g .* etahrcgall[:,:,n])./ Seuler.grid.Δ)
    # ShallowWaters.∂y!(detagdy2, (Seuler.constants.g .* etahrcgall[:,:,n])./ Seuler.grid.Δ)

    # all of the rhs terms are computed with a *scaled* prognostic variable, so to actually 
    # see the momentum budget terms we need to remove that scaling. The division by Delta is because
    # all of the timestepping is done with dt / Delta, and we want the term that gets multiplied by dt
    # I should have corrected for the scalings in the values returned by compute_mom_budget already,
    # but for the stuff I'm pulling after (see below) I need to do it manually

    eulertendency = ((2 .* du) ./ Seuler.grid.dtint);
    parameterization = Seuler.Diag.ZBVars.S_u ./ (Seuler.grid.Δ * Seuler.constants.scale);
    budgetsum = (adv_u + fv - detagdx + Fx + viscu + bottomdragu + 2 .* parameterization);

    fig = Figure(size=(1050, 700), fontsize=15);

    Label(
        fig[0, 3],
        "Online u budget, ZB20 parameterization",
        fontsize = 20,
        tellwidth = false
    )

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    eulertendency,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Euler tendency"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    budgetsum,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Sum of budget terms"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax2)
    Colorbar(fig[1,4], hm2)

    ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (adv_u),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Advection"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax3)
    Colorbar(fig[1,6], hm3)

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    fv,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Coriolis"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    Colorbar(fig[2,2], hm4)
    hidexdecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -detagdx,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Gravity"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    hidedecorations!(ax5)
    Colorbar(fig[2,4], hm5)

    ax6, hm6 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Fx,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Wind stress"),
    # colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax6)
    Colorbar(fig[2,6], hm6)

    ax7, hm7 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    viscu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Viscosity"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    Colorbar(fig[3,2], hm7)

    ax8, hm8 = heatmap(fig[3,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    bottomdragu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Bottom drag"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    hideydecorations!(ax8)
    Colorbar(fig[3,4], hm8)

    ax9, hm9 = heatmap(fig[3,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    parameterization,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Parameterization"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hideydecorations!(ax9)
    Colorbar(fig[3,6], hm9)

end