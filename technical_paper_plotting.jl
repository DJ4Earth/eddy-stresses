"""
Mostly figure generation, I just wanted to be able to run include("technical_paper.jl")
without all of this also running. The function deserialize opens any saved checkpoints/deprecated now
"""




function deserialize(x)
    s = IOBuffer(x)
    Serialization.deserialize(s)
end

function stuff()

    f = Figure()
    ax = Axis(f[1,1])
    h = heatmap!(ax, 0:128:3840,
    0:127:3840,
    temp.u,
    colormap=:balance,
    colorrange=(-maximum(temp.u), maximum(temp.u))
    )
    Colorbar(f[1,2], h)
    
    # investigating derivatives 

    # norm of the prognostic variables derivatives
    detanorm = load_object("./technicalpaper_datafiles/finaldatafiles/tech")
    dunorm = load_object("./technicalpaper_timeaveragedobjective_starting3monthsin_dunorm_dividedbynxny_500m_1year_111424.jld2")
    dvnorm = load_object("./technicalpaper_timeaveragedobjective_starting3monthsin_dvnorm_dividedbynxny_500m_1year_111424.jld2")

    timestep1 = 1:0.99726:364
    f = Figure(size = (1000, 500));
    ax1 = Axis(f[1, 1],
        title = "Norm of adjoint derivative in 100km model, 12-month integration",
        xlabel = "t (days)",
        ylabel = L"||\partial J / \partial u||",
    )
    timestep2 = 1:0.99726:364
    ax2 = Axis(f[2,1],
    # title = "Norm of adjoint derivative in 100km, 12-month integration with double viscosity, spinup also with double viscosity",
    xlabel = "t (days)",
    # yscale=log10,
    ylabel = L"||\partial J / \partial v||",
    )
    timestep3 = 1:0.99726:360
    ax3 = Axis(f[3, 1],
    title = "Norm of adjoint derivative, 12-month integration with double viscosity, spinup also with double viscosity",
    xlabel = "t (days)",
    # yscale = log,
    ylabel = L"||\partial J / \partial x||")

    lines!(ax1, timestep1, dunorm12dv, label = L"||\partial J / \partial u(t)||")
    lines!(ax1, timestep1, dvnorm2, label = L"\partial J / \partial v(t)")
    axislegend(ax1, position = :rt)
    lines!(ax2, timestep2, dvnorm12dv, label = L"||\partial J / \partial v(t)||")
    lines!(ax2, timestep2, dvnorm4, label = L"\partial J / \partial v(t)")
    axislegend(ax2, position = :rt)
    lines!(ax3, timestep3, dunorm12, label = L"\partial J / \partial u(t)")
    lines!(ax3, timestep3, dvnorm12, label = L"\partial J / \partial v(t)")
    axislegend(ax3, position = :rt)

    save("technicalpaper_normprog_divnxny_1000m_114224.png", f)

    # norm of the forcing and parameter derivatives

    dcDnorm = load_object("./technicalpaper_timeaveragedobjective_dcDnorm_dividedbynxny_1000m_1year_111424.jld2");
    dFxnorm = load_object("./technicalpaper_timeaveragedobjective_dFxnorm_dividedbynxny_1000m_1year_1114224.jld2");

    timestep = 365:-0.99726:1
    f = Figure();
    ax1 = Axis(f[1, 1],
        title = "Norm of adjoint derivative w.r.t. wind-stress, time-averaged objective, H = 1000m",
        xlabel = "t (days)",
        ylabel = L"||\partial J / \partial F_x||"
    );
    lines!(ax1, timestep, dFxnorm, label = L"\partial J / \partial F_x")
    axislegend(ax1, position = :rt);

    ax2 = Axis(f[2, 1],
    title = "Norm of adjoint derivative w.r.t. bottom drag, time-averaged objective, H = 1000m",
    xlabel = "t (days)",
    ylabel = L"||\partial J / \partial c_D||");
    lines!(ax2, timestep, dcDnorm, label = L"\partial J / \partial c_D");
    axislegend(ax2, position = :rt);


    save("technicalpaper_parametersensitivities_normcDFx_divnxny_1000mdepth_111424.png", f)

    # computing time averaged velocity/variance/mean flow

    u = ncread("./30km_1year_output/u.nc", "u")
    v = ncread("./30km_1year_output/v.nc", "v")

    avg_u = zeros(127,128)
    avg_v = zeros(128,127)
    for t = 1:366
        avg_u += u[:, :, t]
        avg_v += v[:, :, t]
    end

    f = Figure()
    ax1 = Axis(f[1, 1],
    title = "Average x-velocity",
    xlabel = "x",
    ylabel = "y"
    )
    heatmap!(ax1, 0:127:3840,
        0:128:3840,
        avg_u,
        colorrange = (-150, 150),
        colormap=:balance
    )

    f = Figure()
    ax2 = Axis(f[1, 1],
    title = "Time averaged y-velocity",
    xlabel = "x",
    ylabel = "y"
    )
    h = heatmap!(ax2, 0:128:3840,
    0:127:3840,
    avg_v,
    colorrange = (-600, 600),
    colormap=:balance
    )
    Colorbar(f[1,2], h)
    save("technicalpaper_avgv_111324.png")

    # loss function is final spatially averaged energy
    # initial condition sensitivity

    primal500 = h5open("primal_technicalpaper_check_110524.h5.h5", "r")
    blob = read(primal500["1"])
    prim = deserialize(blob)

    adj50012 = h5open("./technicalpaper_datafiles/finaldatafiles/adjoint_technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_500m_12months_float32_112024.h5", "r")
    day = read(adj50012["225"])
    month = read(adj50012["6721"])
    sixmonth = read(adj50012["40321"])
    twelvemonth = read(adj50012["80641"])

    adj1 = deserialize(day)
    adj30 = deserialize(month)
    adj180 = deserialize(sixmonth)
    adj360 = deserialize(twelvemonth)

    fig = Figure(size=(900, 900));

    dS360 = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(adj360.Prog.u,
    adj360.Prog.v,
    adj360.Prog.η, 
    adj360.Prog.sst,adj360)...)
    ax1 = Axis(fig[1,1], 
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Objective function sensitivity to u -0 month",
    );
    hm1 = CairoMakie.heatmap!(ax1, 0:127:3840,
        0:128:3840,
        dS360.u,
        colormap=:balance,
        colorrange = (-maximum(dS360.u), maximum(dS360.u))
    );
    Colorbar(fig[2, 1], 
        hm1,
        vertical=false
    )

    dS180 = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(adj180.Prog.u,
    adj180.Prog.v,
    adj180.Prog.η, 
    adj180.Prog.sst,adj180)...)
    ax2 = Axis(fig[1,2],
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Objective function sensitivity to u -6 month",
    );
    hm2 = CairoMakie.heatmap!(ax2, 0:128:3840,
        0:128:3840,
        dS180.u,
        colormap=:balance,
        colorrange = (-maximum(dS180.u),maximum(dS180.u))
    );
    Colorbar(fig[2, 2], hm2, vertical=false)

    dS30 = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(adj30.Prog.u,
    adj30.Prog.v,
    adj30.Prog.η, 
    adj30.Prog.sst,adj30)...)
    ax3 = Axis(fig[3,1], 
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Objective function sensitivity to u -11 month",
    );
    hm1 = CairoMakie.heatmap!(ax3, 0:127:3840,
        0:128:3840,
        dS30.u,
        colormap=:balance,
        colorrange = (-maximum(dS30.u),maximum(dS30.u))
    );
    Colorbar(fig[4, 1], 
        hm1,
        vertical=false
    )

    dS1 = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(adj1.Prog.u,
    adj1.Prog.v,
    adj1.Prog.η, 
    adj1.Prog.sst,adj1)...)
    ax4 = Axis(fig[3,2],
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Objective function sensitivity to u -12 month",
    );
    hm2 = CairoMakie.heatmap!(ax4, 0:128:3840,
        0:128:3840,
        dS1.u,
        colormap=:balance,
        colorrange = (-maximum(dS1.u), maximum(dS1.u))
    );
    Colorbar(fig[4, 2], hm2, vertical=false)

    fig

    adj50012 = load_object("./technicalpaper_datafiles/finaldatafiles/technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_startingfromspinup_finaladjoint_struct_500mdepth_12months_float32start_112024.jld2")
    dS50012 = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(adj50012.Prog.u,
    adj50012.Prog.v,
    adj50012.Prog.η, 
    adj50012.Prog.sst,adj50012)...)

    # Initial condition

    fig = Figure(size=(900, 450));

    ax1 = Axis(fig[1,1], 
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = L"\partial J / \partial u(t_0, x, y)",
    );
    hm1 = CairoMakie.heatmap!(ax1, 0:127:3840,
        0:128:3840,
        dS50012.u,
        colormap=:balance,
        colorrange = (-maximum(dS50012.u), maximum(dS50012.u))
    );
    Colorbar(fig[2, 1], 
        hm1,
        vertical=false
    )

    ax2 = Axis(fig[1,2],
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = L"\partial J / \partial v(t_0, x, y)",
    );
    hm2 = CairoMakie.heatmap!(ax2, 0:128:3840,
        0:128:3840,
        dS50012.v,
        colormap=:balance,
        colorrange = (-maximum(dS50012.v),maximum(dS50012.v))
    );
    Colorbar(fig[2, 2], hm2, vertical=false)
    
    ax3 = Axis(fig[1,1], 
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = L"\partial J / \partial F_x",
    );
    hm1 = CairoMakie.heatmap!(ax3, 0:127:3840,
        0:128:3840,
        adj50012.forcing.Fx,
        colormap=:balance,
        colorrange = (-maximum(adj50012.forcing.Fx), maximum(adj50012.forcing.Fx))
    );
    Colorbar(fig[1, 2], 
        hm1,
        vertical=true
    )

    ax4 = Axis(fig[1,1],
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = L"\partial J / \partial c_D",
    );
    hm2 = CairoMakie.heatmap!(ax4, 0:128:3840,
        0:128:3840,
        adj50012.constants.cDfield[2:end-1,2:end-1],
        colormap=:balance,
        colorrange = (-maximum(adj50012.constants.cDfield), maximum(adj50012.constants.cDfield))
    );
    Colorbar(fig[1, 2], hm2, vertical=true)

    fig

    # wind stress sensitivity
    # Makie plot
    fig = Figure(size = (1000,500));
    ax1 = Axis(fig[1,1], 
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Wind-stress sensitivity, H = 500m, time-averaged objective"
    );
    hm1 = CairoMakie.heatmap!(ax1, 0:128:3840,
        0:128:3840,
        adj3.forcing.Fx,
        colormap=:balance,
        colorrange=(-20,20)
    );
    Colorbar(fig[:, end+1], hm1);

    ax2 = Axis(fig[1,3], 
    xlabel = "x (km)",
    ylabel = "y (km)",
    title = "Wind stress sensitivity, H = 500m, non time-averaged objective"
    );
    hm2 = CairoMakie.heatmap!(ax2, 0:128:3840,
    0:128:3840,
    adj2.forcing.Fx,
    colorrange = (-3e8, 3e8),
    colormap=:balance,
    );
    Colorbar(fig[1, 4], hm2);
    fig


    # bottom drag coefficient sensitivity
    fig = Figure(size=(1000,500));
    ax1 = Axis(fig[1,1], 
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Bottom drag sensitivity, H = 500m, time-averaged objective"
    );
    hm1 = CairoMakie.heatmap!(ax1, 0:128:3840,
        0:128:3840,
        adj3.constants.cDfield[2:end-1,2:end-1],
        colorrange = (-600, 600),
        colormap=:balance,
    );
    Colorbar(fig[1, 2], hm1);

    ax2 = Axis(fig[1,3], 
    xlabel = "x (km)",
    ylabel = "y (km)",
    title = "Bottom drag sensitivity, H = 500m, non time-averaged objective"
    );
    hm2 = CairoMakie.heatmap!(ax2, 0:128:3840,
    0:128:3840,
    adj2.constants.cDfield[2:end-1,2:end-1],
    colorrange = (-4e9, 4e9),
    colormap=:balance,
    );
    Colorbar(fig[1, 4], hm2);
    fig

    #### Energy plot 

    t = LinRange(0, 11, 888686)
    t1 = LinRange(0, 1, 80640)
    t2 = LinRange(10, 11, 80641)
    fig = Figure(size=(1000,700));
    ax1 = Axis(fig[1,1], title="Energy during 10-year spinup + 1 year",
    ylabel="Energy")
    ax2 = Axis(fig[2,1], title="Energy during first year of spinup", ylabel="Energy")
    ax3 = Axis(fig[3,1], title="Energy 1-year beyond spinup", ylabel="Energy", xlabel="Years")
    lines!(ax1, t, energy)
    lines!(ax1, t1, energy[1:360*224], label="First year")
    lines!(ax1, t2, energy[808046:end], label="One year beyond spinup")
    axislegend(ax1, position=:rb)
    lines!(ax2, t1, energy[1:360*224])
    lines!(ax3,t2, energy[808046:end])


    fig = Figure();
    ax = Axis(fig[1,1])
    lines!(ax, energy)
    lines!(ax, energy2)
    lines!(ax, S3.parameters.data[738328:end-1])

end

function create_adjoint_gif()

    # primal_fid = h5open("technicalpaper.h5")
    adj_fid = h5open("./technicalpaper_datafiles/finaldatafiles/50km/adjoint_technicalpaper_50km_timeavgobj_onlyfinalmonth_everytimestep_startingfromspinup_everytimestep_12months_float32_120324.h5")
    primal_fid = h5open("./primal_technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_doubleviscosity_startingfromrest_everytimestep_500m_12months_float32_120324.h5")
    # states = ncread("../data_files_gamma0.3/1024_spinup/eta.nc", "eta")

    # unorm_anim = Animation()
    # vnorm_anim = Animation()
    # etanorm_anim = Animation()

    dunorm = []
    dvnorm = []
    detanorm = []

    dcDnorm = []
    dFxnorm = []

    J = []

    # for j = 1:224:224*30*12
    for j = 1:224:47927 # 50 km run
    # for j = 1:66:23964*3 # 100km run
    # for j = 1:87:31546 # 75km run


        blob = read(adj_fid[string(j)])
        adj_chkp = deserialize(blob)

        # push!(J, primal_chkp.parameters.J)


        temp = ShallowWaters.PrognosticVars{Float32}(
            ShallowWaters.remove_halo(adj_chkp.Prog.u,
            adj_chkp.Prog.v,
            adj_chkp.Prog.η,
            adj_chkp.Prog.sst,adj_chkp)...
        )

        push!(dcDnorm, sum(adj_chkp.constants.cDfield.^2) / 128^2)
        push!(dFxnorm, sum(adj_chkp.forcing.Fx.^2) / (127 * 128))
        push!(dunorm, sum(temp.u.^2) / (127*128))
        push!(dvnorm, sum(temp.v.^2) / (127*128))
        push!(detanorm, sum(temp.η.^2) / (128 * 128))


        # frame(eta_anim, heatmap(temp.η',
        #     # title=L"\partial \mathcal{E}(t_f)/\partial u(%$j)",
        #     title=L"\partial \mathcal{E} / \partial \eta(%$j)",
        #     xlabel=L"x",
        #     ylabel=L"y",
        #     c=:balance,
        #     dpi=300,
        #     colorbar_title="         ",
        #     xlabelfontsize=14,
        #     ylabelfontsize=14,
        #     clim=(-15000,15000),
        #     colorbar_titlefontsize=14,
        #     colorbar_tickfontsize=8,
        #     xtickfontsize=8,
        #     ytickfontsize=8)
        # )

        # frame(v_anim, heatmap(temp.v',
        #     title=L"\partial \mathcal{E}(t_f)/\partial v(%$j)",
        #     clim=(-1e9, 1e9),
        #     legend=:none,
        #     xlabel=L"x",
        #     ylabel=L"y",
        #     c=:balance,
        #     dpi=300,
        #     colorbar_title="         ",
        #     xlabelfontsize=14,
        #     ylabelfontsize=14,
        #     colorbar_titlefontsize=14,
        #     colorbar_tickfontsize=8,
        #     xtickfontsize=8,
        #     ytickfontsize=8)
        # )

        # frame(eta_anim, heatmap(temp.η',
        #     title=L"\partial \mathcal{E}(t_f)/\partial \eta(%$j)",
        #     clim=(-50, 50),
        #     legend=:none,
        #     xlabel=L"x",
        #     ylabel=L"y",
        #     c=:balance,
        #     dpi=300,
        #     colorbar_title="         ",
        #     xlabelfontsize=14,
        #     ylabelfontsize=14,
        #     colorbar_titlefontsize=14,
        #     colorbar_tickfontsize=8,
        #     xtickfontsize=8,
        #     ytickfontsize=8)
        # )

    end

    # gif(eta_anim, "eta_integration2_110424.gif", fps = 8)
    # gif(u_anim, "du_integration_365_energy_withclosure_fps7_031424.png", fps = 7)
    # gif(v_anim, "dv_integration_365_energy_withclosure_fps7_031424.png", fps = 7)

    @save "technicalpaper_timeaveragedobjective_finalmonthobj_dcDnormrest_500m_12months_111824.jld2" dcDnorm12rest
    @save "technicalpaper_timeaveragedobjective_finalmonthobj_dFxnormrest_500m_12months_111824.jld2" dFxnorm12rest
    @save "technicalpaper_timeaveragedobjective_finalmonthobj_dunormrest_500m_12months_111824.jld2" dunorm12rest
    @save "technicalpaper_timeaveragedobjective_finalmonthobj_dvnormrest_500m_12months_111824.jld2" dvnorm12rest
    @save "technicalpaper_timeaveragedobjective_finalmonthobj_detanormrest_500m_12months_111824.jld2" detanorm12rest

end

function energy_plot()

    f = Figure();
    ax = Axis(f[1,1], xlabel="Timestep", ylabel="Energy")
    for j = 1:5

        S = ShallowWaters.model_setup(output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        α=2,
        nx=128,
        Ndays=12*30,
        initial_cond="ncfile",
        initpath="./data_files_gamma0.3/10yearspinup_128_noslipbc_fromrest_float32params"
        )

        snaps = Int(floor(sqrt(S.grid.nt)))
        revolve = Revolve{ShallowWaters.ModelSetup}(S.grid.nt, snaps;
            verbose=1,
            gc=true,
            write_checkpoints=false
        )

        S.Prog.u = S.Prog.u + 0.001 .* randn(Float32, 131, 132)
        S.Prog.v = S.Prog.v + 0.001 .* randn(Float32, 132, 131)
        S.Prog.η = S.Prog.η + 0.001 .* randn(Float32, 130, 130)

        _ = checkpointed_integration(S, revolve)

        lines!(ax, S.parameters.data[1:end-1])
        hlines!(ax, sum(S.parameters.data) / S.grid.nt)

    end

end

function fd_plots()

    # first generating and saving the unperturbed objective function over time
    S_unperturbed = ShallowWaters.model_setup(
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
        α=2,
        nx=128,
        Ndays=2*30
        # initial_cond="ncfile",
        # initpath="./run_0001/"
    )

    snaps = Int(floor(sqrt(S_unperturbed.grid.nt)))
    revolve = Revolve{ShallowWaters.ModelSetup}(S_unperturbed.grid.nt,
        snaps;
        verbose=1,
        gc=true,
        write_checkpoints=false
    )

    J_unperturbed = checkpointed_integration(S_unperturbed, revolve)

    # running and storing with a perturbed initial condition, still not using 
    # autodiff or checkpointing

    S_perturbed = ShallowWaters.model_setup(
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
        α=2,
        nx=128,
        Ndays=2*30
        # initial_cond="ncfile",
        # initpath="./run_0001/"
    )

    # adjust the value of the perturbation here
    perturbation = 1e-4
    S_perturbed.Prog.u[62,20] += perturbation

    revolve_perturbed = Revolve{ShallowWaters.ModelSetup}(S_perturbed.grid.nt,
    snaps;
    verbose=1,
    gc=true,
    write_checkpoints=true)

    J_perturbed = checkpointed_integration(S_perturbed, revolve_perturbed)

    t = S_unperturbed.grid.nt-224*30:S_unperturbed.grid.nt
    adj_fid =  h5open("./technicalpaper_datafiles/finaldatafiles/30km/running_from_rest/adjoint_technicalpaper_timeavgobj_onlyfinalmonth_startingfromrest_everytimestep_500m_12months_float32_112024.h5")
    rhs = []
    lhs = []
    for j = (S_unperturbed.grid.nt-224*30):S_unperturbed.grid.nt

        push!(lhs, (S_perturbed.parameters.data[j] - S_unperturbed.parameters.data[j]) / perturbation)

    end

    final_adj = load_object("./technicalpaper_datafiles/finaldatafiles/30km/running_from_rest/technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_startingfromrest_finaladjoint_struct_500mdepth_2months_float32start_112024.jld2")
    final_adj_nohalo = ShallowWaters.PrognosticVars{Float32}(
        ShallowWaters.remove_halo(final_adj.Prog.u,
        final_adj.Prog.v,
        final_adj.Prog.η,
        final_adj.Prog.sst,final_adj)...
    )

    for k = 1:224:S_unperturbed.grid.nt

        blob = read(adj_fid[string(k)])
        adj_chkp = deserialize(blob)

        temp = ShallowWaters.PrognosticVars{Float32}(
            ShallowWaters.remove_halo(adj_chkp.Prog.u,
            adj_chkp.Prog.v,
            adj_chkp.Prog.η,
            adj_chkp.Prog.sst,adj_chkp)...
        )
        push!(rhs, temp.u[62,20])#*perturbation)

    end

    steps = [1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]

    diffs = []

    for s in steps

        S_inner = ShallowWaters.model_setup(
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
            α=2,
            nx=128,
            Ndays=2*30
            # initial_cond="ncfile",
            # initpath="./run_0001/"
        )

        S_inner.Prog.u[62,20] += s

        J_inner = checkpointed_integration(S_inner, revolve_perturbed)

        push!(diffs, (J_inner - J_unperturbed) / s)

    end

end

