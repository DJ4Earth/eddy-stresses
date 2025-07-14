"""
Mostly figure generation, I just wanted to be able to run include("technical_paper.jl")
without all of this also running. The function deserialize opens any saved checkpoints/deprecated now
"""
function deserialize(x)
    s = IOBuffer(x)
    Serialization.deserialize(s)
end

function plots()

    Ndays = 30

    u0 = load_object("coarsegrained_1024_10yearstate_061925.jld2")[1]
    v0 = load_object("coarsegrained_1024_10yearstate_061925.jld2")[2]
    eta0 = load_object("coarsegrained_1024_10yearstate_061925.jld2")[3]

    initial_cond = [u0, v0, eta0]

    Snoparam = ShallowWaters.model_setup(output=false,
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
        Ndays=Ndays
    )
    Snparam.Prog.u .= initial_cond[1]
    Snoparam.Prog.v .= initial_cond[2]
    Snoparam.Prog.η .= initial_cond[3]

    Snn = ShallowWaters.model_setup(output=false,
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
    Snn.Prog.u .= initial_cond[1]
    Snn.Prog.v .= initial_cond[2]
    Snn.Prog.η .= initial_cond[3]

    param_guess = result.minimizer
    current = 1
    for model in (S.Diag.NNVars.model_diag, S.Diag.NNVars.model_offdiag)
        for layers in model[1]
            for array in layers
                sz = prod(size(array))
                array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
                current += sz
            end
        end
    end

    Szb = ShallowWaters.model_setup(output=false,
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
    Szb.Prog.u .= initial_cond[1]
    Szb.Prog.v .= initial_cond[2]
    Szb.Prog.η .= initial_cond[3]

    states_zb = save_states(Szb)
    states_nn = save_states(Snn)

    # energy plots

    t = 6733
    fig1 = Figure(size=(600, 500));
    ax1, hm1 = heatmap(fig1[1,1], (true_states[t].u[:, 1:end-1].^2 .+ true_states[t].v[1:end-1, :].^2),
    colormap=:amp,
    axis=(xlabel=L"x", ylabel=L"y", title=L"\mathcal{E}"),
    colorrange=(0,
    maximum(true_states[t].u[:, 1:end-1].^2 .+ true_states[t].v[1:end-1, :].^2))
    );
    Colorbar(fig1[1,2], hm1)

    t = 224
    fig1 = Figure(size=(800, 700));
    ax1, hm1 = heatmap(fig1[1,1], states_zb[t].u[:, 1:end-1].^2 .+ states_zb[t].v[1:end-1, :].^2,
    colormap=:amp,
    axis=(xlabel=L"x", ylabel=L"y", title=L"\tilde{\mathcal{E}}(+)")
    # colorrange=(0,
    # maximum(true_states[t].u[:, 1:end-1].^2 .+ true_states[t].v[1:end-1, :].^2)),
    );
    Colorbar(fig1[1,2], hm1)

    ax3, hm3 = heatmap(fig1[2, 1], states_nn[t].u[:, 1:end-1].^2 .+ states_nn[t].v[1:end-1, :].^2,
    colormap=:amp,
    axis=(xlabel=L"x", ylabel=L"y", title=L"\tilde{\mathcal{E}}")
    # colorrange=(0,
    # maximum(true_states[t].u[:, 1:end-1].^2 .+ true_states[t].v[1:end-1, :].^2)),
    );
    Colorbar(fig1[2,2], hm3)


    # pre training plots
    u_hr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/u.nc", "u")
    v_hr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/v.nc", "v")

    u_zb = ncread("./spinup_files/128_zbforcingmomentum_withcginitcond_30days_070925/u.nc", "u")
    v_zb = ncread("./spinup_files/128_zbforcingmomentum_withcginitcond_30days_070925/v.nc", "v")

    u_nn = ncread("./spinup_files/128_postspinup_nnforcing_coarsegrainedinitcond_nottrained_10days/u.nc", "u")
    v_nn = ncread("./spinup_files/128_postspinup_nnforcing_coarsegrainedinitcond_nottrained_10days/v.nc", "v")

    totalstates = 224
    up_nn = zeros(65, totalstates)
    vp_nn = zeros(65, totalstates)

    up_zb = zeros(65, totalstates)
    vp_zb = zeros(65, totalstates)

    up_true = zeros(513, totalstates)
    vp_true = zeros(513, totalstates)

    for t = 1:10
        up_nn[:,t] = power(periodogram(u_nn[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_nn[:,t] = power(periodogram(v_nn[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_zb[:,t] = power(periodogram(u_zb[:, :, t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_zb[:,t] = power(periodogram(v_zb[:, :, t]; radialavg=true, radialsum=false)) ./ 128^2

        up_true[:,t] = power(periodogram(u_hr[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2
        vp_true[:,t] = power(periodogram(v_hr[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2
    end

    t = 10
    nn_wl = (1 ./ freq(periodogram(u_zb[:,:,t]; radialavg=true, radialsum=false))) * 30;
    nnu_freq = LinRange(0, 64, 65)
    nnu_freq = nnu_freq ./ 65
    nnu_freq = 1 ./ nnu_freq 
    nnu_freq[1] =  1000
    nn_wl[1] = 1000

    true_wl = 1 ./ freq(periodogram(u_hr[:,:,3]; radialavg=true)) * 3.75;
    true_wl[1] = 1100

    fig2 = Figure(size=(800, 500));
    t = 10
    lines(fig2[1,1], nn_wl[2:end], up_nn[2:end,t] + vp_nn[2:end,t], label="NN", axis=(
            xscale=log10,yscale=log10,xlabel="Wavelength (km)", ylabel="KE(k)", xreversed=true, xticks=[700, 100, 30, 10, 2], title="KE Spectrum before training")
    )
    lines!(fig2[1,1], nn_wl[2:end], up_zb[2:end,t] + vp_zb[2:end,t], label="ZB")
    lines!(fig2[1,1], true_wl[2:end], up_true[2:end,t] + vp_true[2:end,t], label="HR")
    axislegend()


    # after training

    u_hr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/u.nc", "u")
    v_hr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/v.nc", "v")

    u_zb = ncread("./spinup_files/128_zbforcingmomentum_withcginitcond_30days_070925/u.nc", "u")
    v_zb = ncread("./spinup_files/128_zbforcingmomentum_withcginitcond_30days_070925/v.nc", "v")

    u0 = load_object("coarsegrained_1024_10yearstate_061925.jld2")[1]
    v0 = load_object("coarsegrained_1024_10yearstate_061925.jld2")[2]
    eta0 = load_object("coarsegrained_1024_10yearstate_061925.jld2")[3]

    Ndays = 10
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

    S2 = deepcopy(S)

    initial_cond = [u0, v0, eta0]

    S.Prog.u .= initial_cond[1]
    S.Prog.v .= initial_cond[2]
    S.Prog.η .= initial_cond[3]

    S2.Prog.u .= initial_cond[1]
    S2.Prog.v .= initial_cond[2]
    S2.Prog.η .= initial_cond[3]

    current = 1
    for model in (S.Diag.NNVars.model_diag, S.Diag.NNVars.model_offdiag)
        for layers in model[1]
            for array in layers
                sz = prod(size(array))
                array .= reshape(result.minimizer[current:(current + sz - 1)], size(array)...)
                current += sz
            end
        end
    end

    weights = load_object("./tuned_weights/minimizer_134days_statelossfunction_3iterationsLBFGS_dailydata_071125.jld2")
    current = 1
    for model in (S2.Diag.NNVars.model_diag, S2.Diag.NNVars.model_offdiag)
        for layers in model[1]
            for array in layers
                sz = prod(size(array))
                array .= reshape(weights[current:(current + sz - 1)], size(array)...)
                current += sz
            end
        end
    end

    ShallowWaters.time_integration(S2)
    temp2 = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(
    S2.Prog.u,
    S2.Prog.v,
    S2.Prog.η,
    S2.Prog.sst,
    S2
    )...)

    ShallowWaters.time_integration(S)
    temp = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(
    S.Prog.u,
    S.Prog.v,
    S.Prog.η,
    S.Prog.sst,
    S
    )...)

    up_nn = zeros(65)
    vp_nn = zeros(65)

    up2_nn = zeros(65)
    vp2_nn = zeros(65)

    up_zb = zeros(65)
    vp_zb = zeros(65)

    up_true = zeros(513)
    vp_true = zeros(513)
    t = 10
    up_nn[:] = power(periodogram(temp.u; radialavg=true, radialsum=false)) ./ 128^2
    vp_nn[:] = power(periodogram(temp.v; radialavg=true, radialsum=false)) ./ 128^2

    up2_nn[:] = power(periodogram(temp2.u; radialavg=true, radialsum=false)) ./ 128^2
    vp2_nn[:] = power(periodogram(temp2.u; radialavg=true, radialsum=false)) ./ 128^2

    up_zb[:] = power(periodogram(u_zb[:, :, t]; radialavg=true, radialsum=false)) ./ 128^2
    vp_zb[:] = power(periodogram(v_zb[:, :, t]; radialavg=true, radialsum=false)) ./ 128^2

    up_true[:] = power(periodogram(u_hr[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2
    vp_true[:] = power(periodogram(v_hr[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2

    nn_wl = (1 ./ freq(periodogram(u_zb[:,:,t]; radialavg=true, radialsum=false))) * 30;
    nnu_freq = LinRange(0, 64, 65)
    nnu_freq = nnu_freq ./ 65
    nnu_freq = 1 ./ nnu_freq 
    nnu_freq[1] =  1000
    nn_wl[1] = 1000

    true_wl = 1 ./ freq(periodogram(u_hr[:,:,3]; radialavg=true)) * 3.75;
    true_wl[1] = 1100

    fig2 = Figure(size=(800, 500));
    t = 10
    lines(fig2[1,1], nn_wl[2:end], up_nn[2:end] + vp_nn[2:end], label="NN", axis=(
            xscale=log10,yscale=log10,xlabel="Wavelength (km)", ylabel="KE(k)", xreversed=true, xticks=[700, 100, 30, 10, 2], title="KE Spectrum after training")
    )
    lines!(fig2[1,1], nn_wl[2:end], up_zb[2:end] + vp_zb[2:end], label="ZB")
    lines!(fig2[1,1], true_wl[2:end], up_true[2:end] + vp_true[2:end], label="HR")
    lines!(fig2[1,1], nn_wl[2:end], up2_nn[2:end] + vp2_nn[2:end], label="NN states")
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

