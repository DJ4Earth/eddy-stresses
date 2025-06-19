"""
Mostly figure generation, I just wanted to be able to run include("technical_paper.jl")
without all of this also running. The function deserialize opens any saved checkpoints/deprecated now
"""
function deserialize(x)
    s = IOBuffer(x)
    Serialization.deserialize(s)
end

function plots()

    Ndays = 10

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
        handwritten=false,
        N=1,
        α=2,
        nx=128,
        Ndays=Ndays,
        initpath="./data_files_gamma0.3/128_spinup_wforcing_dissipation_wfilter_1pass_noslipbc"
    )

    Snn_kespec = deepcopy(S)
    Snn_energy = deepcopy(S)

    result_kespec = 
    Snn_kespec.Diag.NNVars.model_diag[1][1] .= reshape(result.minimizer[1:34], 2, 17)
    Snn_kespec.Diag.NNVars.model_offdiag[1][1] .= reshape(result.minimizer[35:end], 1, 22)

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
    handwritten=false,
    N=1,
    α=2,
    nx=128,
    Ndays=Ndays,
    initpath="./data_files_gamma0.3/128_spinup_wforcing_dissipation_wfilter_1pass_noslipbc"
    )

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

    t = 673
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

    u_zb = ncread("./spinup_files/128_postspinup_zbforcing_dissipation_10days/u.nc", "u")
    v_zb = ncread("./spinup_files/128_postspinup_zbforcing_dissipation_10days/v.nc", "v")

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
            xscale=log10,yscale=log10,xlabel="Wavelength (km)", ylabel="KE(k)", xreversed=true, xticks=[100, 30, 10, 2], title="KE Spectrum before training")
    )
    lines!(fig2[1,1], nn_wl[2:end], up_zb[2:end,t] + vp_zb[2:end,t], label="ZB")
    lines!(fig2[1,1], true_wl[2:end], up_true[2:end,t] + vp_true[2:end,t], label="HR")
    axislegend()


    # after training

    u_hr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/u.nc", "u")
    v_hr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/v.nc", "v")

    u_zb = ncread("./spinup_files/128_postspinup_zbforcing_dissipation_10days/u.nc", "u")
    v_zb = ncread("./spinup_files/128_postspinup_zbforcing_dissipation_10days/v.nc", "v")

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
        handwritten=false,
        N=1,
        α=2,
        nx=128,
        Ndays=Ndays
    )
    S.Diag.NNVars.model_diag[1][1] .= reshape(result.minimizer[1:34], 2, 17)
    S.Diag.NNVars.model_offdiag[1][1] .= reshape(result.minimizer[35:56], 1, 22)
    S.Diag.NNVars.model_diag[1][2] .= reshape(result.minimizer[57:58], 2, 1)
    S.Diag.NNVars.model_offdiag[1][2] .= result.minimizer[end]

    S.Diag.NNVars.model_diag[1][1] .= reshape(param_guess[1:34], 2, 17)
    S.Diag.NNVars.model_offdiag[1][1] .= reshape(param_guess[35:56], 1, 22)
    S.Diag.NNVars.model_diag[1][2] .= reshape(param_guess[57:58], 2, 1)
    S.Diag.NNVars.model_offdiag[1][2] .= param_guess[end]

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

    up_zb = zeros(65)
    vp_zb = zeros(65)

    up_true = zeros(513)
    vp_true = zeros(513)
    t = 10
    up_nn[:] = power(periodogram(temp.u; radialavg=true, radialsum=false)) ./ 128^2
    vp_nn[:] = power(periodogram(temp.v; radialavg=true, radialsum=false)) ./ 128^2

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
            xscale=log10,yscale=log10,xlabel="Wavelength (km)", ylabel="KE(k)", xreversed=true, title="KE Spectrum after training")
    )
    lines!(fig2[1,1], nn_wl[2:end], up_zb[2:end] + vp_zb[2:end], label="ZB")
    lines!(fig2[1,1], true_wl[2:end], up_true[2:end] + vp_true[2:end], label="HR")
    axislegend()

end

