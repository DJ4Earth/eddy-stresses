"""
Mostly figure generation, I just wanted to be able to run include("technical_paper.jl")
without all of this also running.
"""

function load_and_create_models()

    T = Float64
    Ndays = 30
    coarse_grained_hrstates = load_object("./offline_files/1024_filtered_downsized_uveta_10days_postspinup_hourlysaves_111925.jld2");
    uhrcg = coarse_grained_hrstates[1];
    vhrcg = coarse_grained_hrstates[2];
    etahrcg = coarse_grained_hrstates[3];

    Pnoparam = ShallowWaters.Parameter(T=T,
        output=true,
        output_dt=24,
        L_ratio=1,
        g=9.81,
        H=500,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
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
        nn_forcing_dissipation=false,
        N=1,
        α=2,
        nx=128,
        Ndays=365
    );

    Snoparam = ShallowWaters.model_setup(Pnoparam);

    u0, v0, eta0, _ = ShallowWaters.add_halo(uhrcg[:,:,1],vhrcg[:,:,1],etahrcg[:,:,1],zeros(128,128),Snoparam);
    initial_cond = [u0, v0, eta0];

    Snoparam.Prog.u .= copy(initial_cond[1]);
    Snoparam.Prog.v .= copy(initial_cond[2]);
    Snoparam.Prog.η .= copy(initial_cond[3]);

    ShallowWaters.time_integration(Snoparam);

    PZB = ShallowWaters.Parameter(T=T,
        output=true,
        output_dt=8,
        L_ratio=1,
        g=9.81,
        H=500,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
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
        Ndays=Ndays
    );

    SZB = ShallowWaters.model_setup(PZB);
    SZB.Prog.u .= copy(initial_cond[1]);
    SZB.Prog.v .= copy(initial_cond[2]);
    SZB.Prog.η .= copy(initial_cond[3]);

    ShallowWaters.time_integration(SZB);

    # the offline problem is unstable even for short integrations, Ndays here is thus capped
    Ndays = 3
    Poffline = ShallowWaters.Parameter(T=T,
        output=true,
        output_dt=8,
        L_ratio=1,
        g=9.81,
        H=500,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
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
        Ndays=Ndays
    );

    Soffline = ShallowWaters.model_setup(Poffline);

    offlineweights = load_object("./tuned_weights/result_offline_150iterations_geluactivation_111925.jld2").solution
    current = 1
    for m in (Soffline.Diag.CNNVars.model_Su, Soffline.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(offlineweights[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    Soffline.Prog.u .= copy(initial_cond[1]);
    Soffline.Prog.v .= copy(initial_cond[2]);
    Soffline.Prog.η .= copy(initial_cond[3]);

    ShallowWaters.time_integration(Soffline);

    # now creating the online version, Ndays can be larger
    Ndays = 30
    Ponline = ShallowWaters.Parameter(T=T,
        output=true,
        output_dt=8,
        L_ratio=1,
        g=9.81,
        H=500,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
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
        Ndays=Ndays
    );

    Sonline = ShallowWaters.model_setup(Ponline);

    onlineweights = load_object("./tuned_weights/result_online_madnlp_states_1dayoptimization_100iterations_112125.jld2").solution
    current = 1
    for m in (Sonline.Diag.CNNVars.model_Su, Sonline.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(onlineweights[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    Sonline.Prog.u .= copy(initial_cond[1]);
    Sonline.Prog.v .= copy(initial_cond[2]);
    Sonline.Prog.η .= copy(initial_cond[3]);

    ShallowWaters.time_integration(Sonline)

    coarse_grained_hrstates = load_object("./offline_files/1024_filtered_downsized_uveta_30days_postspinup_8hoursaves_112125.jld2");
    uhrcg = coarse_grained_hrstates[1];
    vhrcg = coarse_grained_hrstates[2];
    etahrcg = coarse_grained_hrstates[3];

    uhr = ncread("./spinup_files/1024_postspinup_30days_8hoursaves/u.nc", "u");
    vhr = ncread("./spinup_files/1024_postspinup_30days_8hoursaves/v.nc", "v");
    etahr = ncread("./spinup_files/1024_postspinup_30days_8hoursaves/eta.nc", "eta");

    uofflinegelu = ncread("./results/128_offlineparam_postspinup_cginitcond_3days_8hoursaves/u.nc", "u");
    vofflinegelu = ncread("./results/128_offlineparam_postspinup_cginitcond_3days_8hoursaves/v.nc", "v");
    etaofflinegelu = ncread("./results/128_offlineparam_postspinup_cginitcond_3days_8hoursaves/eta.nc", "eta");

    unoparam = ncread("./results/128_noparam_postspinup_cginitcond_30days_8hoursaves/u.nc", "u");
    vnoparam = ncread("./results/128_noparam_postspinup_cginitcond_30days_8hoursaves/v.nc", "v");
    etanoparam = ncread("./results/128_noparam_postspinup_cginitcond_30days_8hoursaves/eta.nc", "eta");

    uonlinegelu = ncread("./results/128_online_geluactivation_offlineinitweights_madnlp_30days_8hoursaves/u.nc", "u");
    vonlinegelu = ncread("./results/128_online_geluactivation_offlineinitweights_madnlp_30days_8hoursaves/v.nc", "v");
    etaonlinegelu = ncread("./results/128_online_geluactivation_offlineinitweights_madnlp_30days_8hoursaves/eta.nc", "eta");

    uzb = ncread("./results/128_ZBparam_postspinup_cginitcond_30days_8hoursaves/u.nc", "u");
    vzb = ncread("./results/128_ZBparam_postspinup_cginitcond_30days_8hoursaves/v.nc", "v");
    etazb = ncread("./results/128_ZBparam_postspinup_cginitcond_30days_8hoursaves/eta.nc", "eta");

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))
    # imfilter(hru[:,:,j], reflect(ker))

end

function plots()

    # energy plots ############################################################

    # nc files
    # u_zb, uhr
    # jld2 files (my save states function)
    # states_noparam, states_nn (untrained), states_trainednn_pd, states_trainednn_pd65, states_trainednn_kespec, states_trainednn_states

    # to get coarse-grained states

    # t is timestep, and I saved every 8 hours up to 30 days
    # this means t can be anything between 1 (the initial condition) and 91 (the final step after 30 days)

    # Prognostic variables #############################################################

    # for showing offline instability
    t = 10
    fig = Figure(size=(950, 275), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etaofflinegelu[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\eta(3 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(etaofflinegelu[:,:,t])),maximum(abs.(etaofflinegelu[:,:,t])))
    );
    Colorbar(fig[1,2], hm1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uofflinegelu[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u(3 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(uofflinegelu[:,:,t])),maximum(abs.(uofflinegelu[:,:,t])))
    );
    Colorbar(fig[1,4], hm2)

    ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    vofflinegelu[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"v(3 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(vofflinegelu[:,:,t])),maximum(abs.(vofflinegelu[:,:,t])))
    );
    Colorbar(fig[1,6], hm3)

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)"], [ga, gb,gc])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # high versus low resolution eta
    t = 1
    fig = Figure(size=(800, 350), fontsize=15);
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etahr[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\mathbf{\eta}(3650 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(etahr[:,:,t])),maximum(abs.(etahr[:,:,t])))
    );
    Colorbar(fig[1,2], hm1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etahrcg[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\overline{\mathbf{\eta}}(3650 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(etahrcg[:,:,t])),maximum(abs.(etahrcg[:,:,t])))
    );
    Colorbar(fig[1,4], hm2)

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)"], [ga, gb])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end


    # just u fields
    t = 46
    fig = Figure(size=(900, 800), fontsize=15);
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uhrcg[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\overline{u}(15 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[1,2], hm1, label="m")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    unoparam[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u(15 \; \text{days}, x, y)\text{, no closure}"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[1,4], hm2, label="m")

    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uzb[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u(15 \; \text{days}, x, y)\text{, ZB20 closure}"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[2,2], hm3, label="m")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uonlinegelu[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u(15 \; \text{days}, x, y)\text{, online closure}"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[2,4], hm4, label="m")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[2, 1] = GridLayout()
    gd = fig[2, 3] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)"], [ga, gb, gc, gd])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # u and v fields
    t = 46
    fig = Figure(size=(1800, 800), fontsize=15);
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uhrcg[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Filtered, coarse-grained u(15 days, x, y)"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[1,2], hm1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    unoparam[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="u(15 days, x, y), no closure"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[1,4], hm2)

    ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uzb[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="u(15 days, x, y), ZB closure"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[1,6], hm3)

    ax4, hm4 = heatmap(fig[1,7], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uonlinegelu[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="u(15 days, x, y), online closure"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[1,8], hm4)

    ax1, hm1 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    vhrcg[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Filtered, coarse-grained v(15 days, x, y)"),
    colorrange=(-maximum(abs.(vhrcg[:,:,t])),maximum(abs.(vhrcg[:,:,t])))
    );
    Colorbar(fig[2,2], hm1)

    ax2, hm2 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    vnoparam[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="v(15 days, x, y), no closure"),
    colorrange=(-maximum(abs.(vhrcg[:,:,t])),maximum(abs.(vhrcg[:,:,t])))
    );
    Colorbar(fig[2,4], hm2)

    ax3, hm3 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    vzb[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="v(15 days, x, y), ZB closure"),
    colorrange=(-maximum(abs.(vhrcg[:,:,t])),maximum(abs.(vhrcg[:,:,t])))
    );
    Colorbar(fig[2,6], hm3)

    ax4, hm4 = heatmap(fig[2,7], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    vonlinegelu[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="v(15 days, x, y), online closure"),
    colorrange=(-maximum(abs.(vhrcg[:,:,t])),maximum(abs.(vhrcg[:,:,t])))
    );
    Colorbar(fig[2,8], hm4)

    # eta fields
    t = 46
    fig = Figure(size=(900, 800), fontsize=15);
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etahrcg[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Filtered, coarse-grained eta(15 days, x, y)"),
    colorrange=(-maximum(abs.(etahrcg[:,:,t])),maximum(abs.(etahrcg[:,:,t])))
    );
    Colorbar(fig[1,2], hm1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etanoparam[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="eta(15 days, x, y), no closure"),
    colorrange=(-maximum(abs.(etahrcg[:,:,t])),maximum(abs.(etahrcg[:,:,t])))
    );
    Colorbar(fig[1,4], hm2)

    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etazb[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="eta(15 days, x, y), ZB closure"),
    colorrange=(-maximum(abs.(etahrcg[:,:,t])),maximum(abs.(etahrcg[:,:,t])))
    );
    Colorbar(fig[2,2], hm3)

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etaonlinegelu[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="eta(15 days, x, y), online closure"),
    colorrange=(-maximum(abs.(etahrcg[:,:,t])),maximum(abs.(etahrcg[:,:,t])))
    );
    Colorbar(fig[2,4], hm4)

    # time-series of prognostic field computed with online parameterization
    t = [4, 10, 46, 91]
    fig = Figure(size=(900, 800), fontsize=15);
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uonlinegelu[:,:,t[1]],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u(1 \text{ day}, x, y)"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[1,2], hm1, label="m/s")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uonlinegelu[:,:,t[2]],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u(3 \text{ days}, x, y)"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[1,4], hm2, label="m/s")

    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uonlinegelu[:,:,t[3]],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u(15 \text{ days}, x, y)"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[2,2], hm3, label="m/s")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uonlinegelu[:,:,t[4]],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u(30 \text{ days}, x, y)"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[2,4], hm4, label="m/s")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[2, 1] = GridLayout()
    gd = fig[2, 3] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)"], [ga, gb, gc, gd])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # Energy ############################################################################

    # high-resolution versus coarse-grained high resolution energy
    t = 31
    fig = Figure(size=(800, 350), fontsize=15);
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uhr[:,1:end-1,t].^2 .+ vhr[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="3.75 km resolution E(10 days, x, y)"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[1,2], hm1)

    ax1, hm1 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Coarse-grained E(10 days, x, y)"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[1,4], hm1)

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)"], [ga, gb])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # coarse-grained versus no parameterization
    t = 31
    fig = Figure(size=(800, 400), fontsize=15);
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Coarse-grained energy"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[1,2], hm1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 1024),
    LinRange(0, 3840, 1024),
    (unoparam[:,1:end-1,t].^2 .+ vnoparam[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30 km, no parameterization"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[1,4], hm1)

    # cg, zb, nn, no param
    t = 31
    fig = Figure(size=(900, 1000), fontsize=15);
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Filtered high-resolution energy"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[1,2], hm1)

    ax1, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (unoparam[:,1:end-1,t].^2 .+ vnoparam[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30km resolution E, no closure"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[1,4], hm1)

    ax1, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uzb[:,1:end-1,t].^2 .+ vzb[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30 km resolution E, ZB closure"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[2,2], hm1)

    ax1, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uonlinegelu[:,1:end-1,t].^2 .+ vonlinegelu[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30 km resolution E, online closure"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[2,4], hm1)

    ax1, hm5 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uofflinegelu[:,1:end-1,10].^2 .+ vofflinegelu[1:end-1,:,10].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30 km resolution E, offline closure after 3 days"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[3,2], hm1)

    # same as the above but without the coarse-grained energy
    # cg, zb, nn, no param
    t = 31
    fig = Figure(size=(900, 800), fontsize=15);

    ax1, hm2 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (unoparam[:,1:end-1,t].^2 .+ vnoparam[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30km resolution E(10 days, x, y), no closure"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[1,2], hm1)

    ax1, hm3 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uzb[:,1:end-1,t].^2 .+ vzb[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30 km resolution E(10 days, x, y), ZB closure"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[1,4], hm1)

    ax1, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uonlinegelu[:,1:end-1,t].^2 .+ vonlinegelu[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30 km resolution E(10 days, x, y), online closure"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[2,2], hm1)

    ax1, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uofflinegelu[:,1:end-1,10].^2 .+ vofflinegelu[1:end-1,:,10].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30 km resolution E(3 days, x, y), offline closure"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[2,4], hm1)

    # spatially averaged energy over integration (30 day integrations on all)

    gelu = []
    noparam = []
    zb = []
    cghr = []
    for j = 1:31
        push!(gelu, sum(uonlinegelu[:,1:end-1,j].^2 .+ vonlinegelu[1:end-1,:,j].^2))
        push!(zb, sum(uzb[:,1:end-1,j].^2 .+ vzb[1:end-1,:,j].^2))
        push!(cghr, sum(uhrcg[:,1:end-1,j].^2 .+ vhrcg[1:end-1,:,j].^2))
        push!(noparam, sum(unoparam[:,1:end-1,j].^2 .+ vnoparam[1:end-1,:,j].^2))
    end

    fig = Figure(size=(1000, 500), fontsize=15);
    lines(fig[1,1], LinRange(0,10, 31), cghr, label="Coarse-grained 3.75km resolution", 
        axis=(
            xlabel="Day",
            ylabel="Energy",
            title="Spatially averaged energy"
        )
    )
    lines!(fig[1,1], LinRange(0,10, 31), noparam, label="30km resolution, no closure")
    lines!(fig[1,1], LinRange(0,10, 31), zb, label="ZB closure")
    lines!(fig[1,1], LinRange(0,10, 31), gelu, label="Online closure")
    axislegend(position = (0,0))


    ###################################################################################

    # KE spectrum #############################################################

    # to get coarse-grained states
    ker = ImageFiltering.Kernel.gaussian((30e3/3750))
    # imfilter(hru[:,:,j], reflect(ker))

    uhr = ncread("./spinup_files/1024_postspinup_30days_8hoursaves/u.nc", "u");
    vhr = ncread("./spinup_files/1024_postspinup_30days_8hoursaves/v.nc", "v");
    etahr = ncread("./spinup_files/1024_postspinup_30days_8hoursaves/eta.nc", "eta");

    coarse_grained_hrstates = load_object("./offline_files/1024_filtered_downsized_uveta_30days_postspinup_8hoursaves_112125.jld2");
    uhrcg = coarse_grained_hrstates[1];
    vhrcg = coarse_grained_hrstates[2];
    etahrcg = coarse_grained_hrstates[3];

    filtered_hrstates = load_object("./offline_files/1024_filtered_uveta_imfilter_30days_postspinup_8hoursaves_112125.jld2");
    uhrfilter= filtered_hrstates[1];
    vhrfilter = filtered_hrstates[2];
    etahrfilter = filtered_hrstates[3];

    totalstates = 91 # saved every 8 hours (this is when the hr cg and low resolution match up)
    up_noparam = zeros(65,totalstates)
    vp_noparam = zeros(65,totalstates)

    up_hr = zeros(513, totalstates)
    vp_hr = zeros(513, totalstates)

    up_zb = zeros(65,totalstates)
    vp_zb = zeros(65,totalstates)

    up_hrfilter = zeros(513,totalstates)
    vp_hrfilter = zeros(513,totalstates)

    up_hrcg = zeros(65,totalstates)
    vp_hrcg = zeros(65,totalstates)

    up_nn = zeros(65,totalstates)
    vp_nn = zeros(65,totalstates)

    for t = 1:totalstates

        up_hr[:,t] = power(periodogram(uhr[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2
        vp_hr[:,t] = power(periodogram(vhr[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2

        up_zb[:,t] = power(periodogram(uzb[:, :, t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_zb[:,t] = power(periodogram(vzb[:, :, t]; radialavg=true, radialsum=false)) ./ 128^2

        up_hrfilter[:,t] = power(periodogram(uhrfilter[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2
        vp_hrfilter[:,t] = power(periodogram(vhrfilter[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2

        up_hrcg[:,t] = power(periodogram(uhrcg[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_hrcg[:,t] = power(periodogram(vhrcg[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_noparam[:,t] = power(periodogram(unoparam[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_noparam[:,t] = power(periodogram(vnoparam[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_nn[:,t] = power(periodogram(uonlinegelu[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_nn[:,t] = power(periodogram(vonlinegelu[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

    end

    lr_wl = (1 ./ freq(periodogram(uzb[:,:,10]; radialavg=true, radialsum=false))) * 30;
    nnu_freq = LinRange(0, 64, 65)
    nnu_freq = nnu_freq ./ 65
    nnu_freq = 1 ./ nnu_freq 
    nnu_freq[1] =  1000
    lr_wl[1] = 1000

    hr_wl = 1 ./ freq(periodogram(uhr[:,:,3]; radialavg=true)) * 3.75;
    hr_wl[1] = 1100

    fig = Figure(size=(1000, 500), fontsize=15);
    t = 91
    lines(fig[1,1], hr_wl[2:65], up_hrfilter[2:65,t] + vp_hrfilter[2:65,t], label="Filtered 3.75km resolution", axis=(
            xscale=log10,
            yscale=log10,
            xlabel="Wavelength (km)",
            ylabel="KE(k)",
            xreversed=true,
            xticks=[700, 100, 30, 10, 2],
            title="KE spectrum")
    )
    lines!(fig[1,1], lr_wl[2:end], up_nn[2:end,t] + vp_nn[2:end,t], label="Online NN closure"
    )
    lines!(fig[1,1], lr_wl[2:end], up_noparam[2:end,t] + vp_noparam[2:end,t], label="30 km resolution, no closure")
    lines!(fig[1,1], lr_wl[2:end], up_zb[2:end,t] + vp_zb[2:end,t], label="ZB closure")
    axislegend(position = (0,0))

    #### comparing time-averaged ke spectra

    up_noparam_avg = zeros(65)
    vp_noparam_avg = zeros(65)

    up_zb_avg = zeros(65)
    vp_zb_avg = zeros(65)

    up_hr_avg = zeros(513)
    vp_hr_avg = zeros(513)

    up_cghr_avg = zeros(65)
    vp_cghr_avg = zeros(65)

    up_filter_avg = zeros(513)
    vp_filter_avg = zeros(513)

    up_nn_avg = zeros(65)
    vp_nn_avg = zeros(65)

    for t = 2:91

        up_noparam_avg += up_noparam[:,t]
        vp_noparam_avg += vp_noparam[:,t]

        up_nn_avg += up_nn[:,t]
        vp_nn_avg += vp_nn[:,t]

        up_zb_avg += up_zb[:,t]
        vp_zb_avg += vp_zb[:,t]

        up_hr_avg += up_hr[:,t]
        vp_hr_avg += vp_hr[:,t]

        up_filter_avg += up_hrfilter[:,t]
        vp_filter_avg += vp_hrfilter[:,t]

        up_cghr_avg += up_cghr[:,t]
        vp_cghr_avg += vp_cghr[:,t]

    end

    fig = Figure(size=(1000, 500), fontsize=15);
    lines(fig[1,1], hr_wl[2:65], (up_cghr_avg[2:65] + vp_cghr_avg[2:65])/31, label="Filtered, coarse-grained 3.75km resolution", axis=(
        xscale=log10,
        yscale=log10,
        xlabel="Wavelength (km)",
        ylabel="KE(k)", xreversed=true,
        xticks=[700, 100, 30, 10, 2],
        title="30 day averaged KE spectrum")
    )
    lines!(fig[1,1], lr_wl[2:end], (up_zb_avg[2:end] + vp_zb_avg[2:end])/31, label="ZB closure")
    lines!(fig[1,1], lr_wl[2:end], (up_noparam_avg[2:end] + vp_noparam_avg[2:end])/totalstates, label="30 km resolution, no closure")
    lines!(fig[1,1], lr_wl[2:end], (up_nn_avg[2:end] + vp_nn_avg[2:end])/totalstates, label="Online NN closure")
    axislegend(position = (0,0))

    #############################################################################################


    # comparing cg energy spectra to energy spectra of the cg

    # up_hr, up_hrfilter, up_hrcg

    true_wl = (1 ./ freq(periodogram(uhr[:,:,1]; radialavg=true))) * 3.75;
    true_wl[1] = 1100
    cg_wl = (1 ./ freq(periodogram(uhrcg[:,:,1]; radialavg=true, radialsum=false))) * 30;
    cg_wl[1] = 1100

    t = 1
    fig = Figure();
    lines(fig[1,1], true_wl[2:end], up_hr[2:end,t] + vp_hr[2:end,t], label="3.75 km resolution", axis=(
            xscale=log10,yscale=log10,xlabel="Wavelength (km)", ylabel="KE(k)", xreversed=true, xticks=[700, 100, 30, 10, 2], title="HR energy spectrum")
    )
    lines!(fig[1,1], cg_wl[2:end], up_hrcg[2:end,t] + vp_hrcg[2:end,t], label="Filtered, coarsened HR spectrum")
    lines!(fig[1,1], true_wl[2:end], up_hrfilter[2:end,t] + vp_hrfilter[2:end,t], label="Filtered HR spectrum")
    axislegend()

    figu = Figure();
    lines(figu[1,1], true_wl[2:65], up_hr[2:65,t], label="HR", axis=(
            xscale=log10,yscale=log10, ylabel="KE(k)", xreversed=true,  title="HR")
    )
    lines!(figu[1,1], cg_wl[2:end], up_hrcg[2:end,t], label="Coarse-grained HR")
    lines!(figu[1,1], true_wl[2:end], up_hrfilter[2:end,t], label="Filtered HR")
    axislegend()

    figv = Figure();
        lines(figv[1,1], true_wl[2:65], vp_hr[2:65,t], label="HR", axis=(
            xscale=log10,yscale=log10, ylabel="KE(k)", xreversed=true,  title="HR")
    )
    lines!(figv[1,1], cg_wl[2:end], vp_hrcg[2:end,t], label="Coarse-grained HR")
    lines!(figv[1,1], true_wl[2:65], vp_hrfilter[2:65,t], label="Filtered HR")
    axislegend()

    figeta = Figure();
        lines(figeta[1,1], true_wl[2:65], etap_hr[2:65,t], label="HR", axis=(
            xscale=log10,yscale=log10, ylabel="KE(k)", xreversed=true,  title="HR")
    )
    lines!(figeta[1,1], cg_wl[2:end], etap_hrcg[2:end,t], label="Coarse-grained HR")
    lines!(figeta[1,1], cg_wl[2:end], etap_hrfilter[2:end,t], label="Filtered HR")
    axislegend()


end

function longer_integration_kespec()

    Ndays = 365

    u_hr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/u.nc", "u")
    v_hr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/v.nc", "v")

    initial_cond = load_object("./coarsegrained_1024_10yearstate_061925.jld2")
    result = load_object("./tuned_weights/weights_aftertraining_kespectrum_twooptimiterations_onetimeseries_10dayintegraton_dailydataduringfinalweek_070925.jld2")

    S_zb = ShallowWaters.model_setup(output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
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
        Ndays=Ndays
    )

    S_zb.Prog.u .= initial_cond[1]
    S_zb.Prog.v .= initial_cond[2]
    S_zb.Prog.η .= initial_cond[3]

    S_before = ShallowWaters.model_setup(output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
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
        Ndays=Ndays
    )

    S_before.Prog.u .= initial_cond[1]
    S_before.Prog.v .= initial_cond[2]
    S_before.Prog.η .= initial_cond[3]

    S_after = ShallowWaters.model_setup(output=true,
        L_ratio=1,
        g=9.81,
        H=500,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
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
        Ndays=Ndays
    )

    S_after.Prog.u .= initial_cond[1]
    S_after.Prog.v .= initial_cond[2]
    S_after.Prog.η .= initial_cond[3]

    current = 1
    for model in (S_after.Diag.NNVars.model_diag, S_after.Diag.NNVars.model_offdiag)
        for layers in model[1]
            for array in layers
                sz = prod(size(array))
                array .= reshape(result.minimizer[current:(current + sz - 1)], size(array)...)
                current += sz
            end
        end
    end

    ShallowWaters.time_integration(S_after)
    ShallowWaters.time_integration(S_before)
    ShallowWaters.time_integration(S_zb)

    temp_before = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(
    S_before.Prog.u,
    S_before.Prog.v,
    S_before.Prog.η,
    S_before.Prog.sst,
    S_before
    )...)
    temp_after = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(
    S_after.Prog.u,
    S_after.Prog.v,
    S_after.Prog.η,
    S_after.Prog.sst,
    S_after
    )...)
    temp_zb = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(
    S_zb.Prog.u,
    S_zb.Prog.v,
    S_zb.Prog.η,
    S_zb.Prog.sst,
    S_zb
    )...)

    up_before = zeros(65)
    vp_before = zeros(65)

    up_nn = zeros(65)
    vp_nn = zeros(65)

    up_zb = zeros(65)
    vp_zb = zeros(65)

    up_true = zeros(513)
    vp_true = zeros(513)

    up_before[:] = power(periodogram(temp_before.u; radialavg=true, radialsum=false)) ./ 128^2
    vp_before[:] = power(periodogram(temp_before.v; radialavg=true, radialsum=false)) ./ 128^2

    up_nn[:] = power(periodogram(temp_after.u; radialavg=true, radialsum=false)) ./ 128^2
    vp_nn[:] = power(periodogram(temp_after.v; radialavg=true, radialsum=false)) ./ 128^2

    up_zb[:] = power(periodogram(temp_zb.u; radialavg=true, radialsum=false)) ./ 128^2
    vp_zb[:] = power(periodogram(temp_zb.v; radialavg=true, radialsum=false)) ./ 128^2

    t = 365
    up_true[:] = power(periodogram(u_hr[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2
    vp_true[:] = power(periodogram(v_hr[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2

    nn_wl = (1 ./ freq(periodogram(temp_zb.u; radialavg=true, radialsum=false))) * 30;
    nnu_freq = LinRange(0, 64, 65)
    nnu_freq = nnu_freq ./ 65
    nnu_freq = 1 ./ nnu_freq 
    nnu_freq[1] =  1000
    nn_wl[1] = 1000

    true_wl = 1 ./ freq(periodogram(u_hr[:,:,3]; radialavg=true)) * 3.75;
    true_wl[1] = 1100

    fig2 = Figure(size=(800, 500));
    lines(fig2[1,1], nn_wl[2:end], up_before[2:end] + vp_before[2:end], label="NN before training", axis=(
            xscale=log10,yscale=log10,xlabel="Wavelength (km)", ylabel="KE(k)", xreversed=true, title="KE Spectrum after training")
    )
    lines!(fig2[1,1], nn_wl[2:end], up_nn[2:end] + vp_nn[2:end], label="NN after training")
    lines!(fig2[1,1], nn_wl[2:end], up_zb[2:end] + vp_zb[2:end], label="ZB")
    lines!(fig2[1,1], true_wl[2:end], up_true[2:end] + vp_true[2:end], label="HR")
    axislegend()

end

