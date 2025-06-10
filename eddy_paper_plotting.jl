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

    Snn.Diag.NNVars.model_diag[1][1] .= reshape(result.minimizer[1:34], 2, 17)
    Snn.Diag.NNVars.model_offdiag[1][1] .= reshape(result.minimizer[35:end], 1, 22)

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

    up_nn = zeros(65, 673)
    vp_nn = zeros(65, 673)

    up_zb = zeros(65, 673)
    vp_zb = zeros(65, 673)

    up_true = zeros(65, 673)
    vp_true = zeros(65, 673)

    for t = 1:673
        up_nn[:,t] = power(periodogram(states_nn[t].u; radialavg=true))
        vp_nn[:,t] = power(periodogram(states_nn[t].v; radialavg=true))

        up_zb[:,t] = power(periodogram(states_zb[t].u; radialavg=true))
        vp_zb[:,t] = power(periodogram(states_zb[t].v; radialavg=true))

        # up_true[:,t] = power(periodogram(true_states[t].u; radialavg=true))
        # vp_true[:,t] = power(periodogram(true_states[t].v; radialavg=true))
    end

    fftu_nn = fft(up_nn,[2])
    fftv_nn = fft(vp_nn,[2])
    nn_wl = 1 ./ freq(periodogram(states_nn[3].u; radialavg=true));
    nnu_freq = LinRange(0, 672, 673)
    nnu_freq = nnu_freq ./ 673
    nnu_freq = 1 ./ nnu_freq 
    nnu_freq[1] =  1000
    nn_wl[1] = 1000

    fig2 = Figure(size=(800, 500));
    t = 673
    lines(fig2[1,1], nn_wl[2:end], up_nn[2:end,t] + vp_nn[2:end,t], label="NN", axis=(xscale=log10,yscale=log10,xlabel="Wavelength (km)", ylabel="KE(k)", xreversed=true, xticks=[100, 30, 10, 2]))
    lines!(fig2[1,1], nn_wl[2:end], up_zb[2:end,t] + vp_zb[2:end,t], label="ZB")
    lines!(fig2[1,1], true_wl[2:end], up_true[2:end,t] + vp_true[2:end,t], label="Truth")
    axislegend()

end

