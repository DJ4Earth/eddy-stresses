
function spectrum_plots()

    # KE spectrum #############################################################

    # to get coarse-grained states
    ker = ImageFiltering.Kernel.gaussian((30e3/3750))
    # imfilter(hru[:,:,j], reflect(ker))

    # for first three years
    totalstates = 1096

    # for all ten years
    totalstates = 1461
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

    up_gelu1day = zeros(65,totalstates)
    vp_gelu1day = zeros(65,totalstates)

    up_gelu5day = zeros(65,totalstates)
    vp_gelu5day = zeros(65,totalstates)

    up_10day = zeros(65, totalstates)
    vp_10day = zeros(65, totalstates)

    up_20day = zeros(65, totalstates)
    vp_20day = zeros(65, totalstates)

    up_30day = zeros(65, totalstates)
    vp_30day = zeros(65, totalstates)

    up_multi10 = zeros(65,totalstates)
    vp_multi10 = zeros(65,totalstates)

    up_multi20 = zeros(65,totalstates)
    vp_multi20 = zeros(65,totalstates)

    up_multi3 = zeros(65,totalstates)
    vp_multi3 = zeros(65,totalstates)

    up_multi3more = zeros(65,totalstates)
    vp_multi3more = zeros(65,totalstates)

    up_multi1 = zeros(65,totalstates)
    vp_multi1 = zeros(65,totalstates)

    up_multi2 = zeros(65,totalstates)
    vp_multi2 = zeros(65,totalstates)

    up_geluKEspecpd = zeros(65,totalstates)
    vp_geluKEspecpd = zeros(65,totalstates)

    up_geluKEspec = zeros(65,totalstates)
    vp_geluKEspec = zeros(65,totalstates)

    up_geluhybrid = zeros(65,totalstates)
    vp_geluhybrid = zeros(65,totalstates)

    up_gelufourier = zeros(65,totalstates)
    vp_gelufourier = zeros(65,totalstates)

    up_5dayeta = zeros(65, 366)
    vp_5dayeta = zeros(65, 366)

    up_10dayeta = zeros(65, 366)
    vp_10dayeta = zeros(65, 366)

    # up_relu1day = zeros(65,totalstates)
    # vp_relu1day = zeros(65,totalstates)

    # up_relu5day = zeros(65,totalstates)
    # vp_relu5day = zeros(65,totalstates)

    # up_reluKEspecpd = zeros(65,totalstates)
    # vp_reluKEspecpd = zeros(65,totalstates)

    @views for t = 1:totalstates

        # up_hr[:,t] = power(periodogram(uhr[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2
        # vp_hr[:,t] = power(periodogram(vhr[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2

        up_zb[:,t] = power(periodogram(uzball[:, :, t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_zb[:,t] = power(periodogram(vzball[:, :, t]; radialavg=true, radialsum=false)) ./ 128^2

        up_hrcg[:,t] = power(periodogram(uhrcgall[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_hrcg[:,t] = power(periodogram(vhrcgall[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_20day[:,t] = power(periodogram(u20sall[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_20day[:,t] = power(periodogram(v20sall[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_30day[:,t] = power(periodogram(u30sall[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_30day[:,t] = power(periodogram(v30sall[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_noparam[:,t] = power(periodogram(unoparamall[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_noparam[:,t] = power(periodogram(vnoparamall[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        # up_gelu5day[:,t] = power(periodogram(u5s[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        # vp_gelu5day[:,t] = power(periodogram(v5s[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        # up_multi1[:,t] = power(periodogram(umulti1[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        # vp_multi1[:,t] = power(periodogram(vmulti1[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_multi2[:,t] = power(periodogram(umulti2all[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_multi2[:,t] = power(periodogram(vmulti2all[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_multi10[:,t] = power(periodogram(umulti10all[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_multi10[:,t] = power(periodogram(vmulti10all[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_multi20[:,t] = power(periodogram(umulti20all[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_multi20[:,t] = power(periodogram(vmulti20all[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_multi3[:,t] = power(periodogram(umulti3all[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_multi3[:,t] = power(periodogram(vmulti3all[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        # up_multi3more[:,t] = power(periodogram(umulti3more[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        # vp_multi3more[:,t] = power(periodogram(vmulti3more[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        # up_gelu1day[:,t] = power(periodogram(u1daystategelu[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        # vp_gelu1day[:,t] = power(periodogram(v1daystategelu[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        # up_geluKEspecpd[:,t] = power(periodogram(ukespecpd[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        # vp_geluKEspecpd[:,t] = power(periodogram(vkespecpd[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        # up_geluKEspec[:,t] = power(periodogram(ukespec[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        # vp_geluKEspec[:,t] = power(periodogram(vkespec[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        # up_geluhybrid[:,t] = power(periodogram(uhybrid[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        # vp_geluhybrid[:,t] = power(periodogram(vhybrid[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        # up_gelufourier[:,t] = power(periodogram(ufourier[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        # vp_gelufourier[:,t] = power(periodogram(vfourier[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        # up_relu1day[:,t] = power(periodogram(u1daystaterelu[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        # vp_relu1day[:,t] = power(periodogram(v1daystaterelu[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        # up_relu5day[:,t] = power(periodogram(u5daystaterelu[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        # vp_relu5day[:,t] = power(periodogram(v5daystaterelu[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        # up_reluKEspecpd[:,t] = power(periodogram(ukespecpdrelu[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        # vp_reluKEspecpd[:,t] = power(periodogram(vkespecpdrelu[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

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
    t = 360
    ax = Axis(
        xscale=log10,
        yscale=log10,
        xlabel="Wavelength (km)",
        ylabel="KE(k)",
        xreversed=true,
        xticks=[700, 100, 30, 10, 2],
        title="KE spectrum"
    )
    lines!(ax, lr_wl[2:65], up_hrcg[2:65,t] + vp_hrcg[2:65,t], label="Filtered, coarse-grained 3.75km resolution")
    # lines!(fig[1,1], lr_wl[2:end], up_gelu1day[2:end,t] + vp_gelu1day[2:end,t], label="Online NN closure, 1 day")
    lines!(ax, lr_wl[2:end], up_gelu5day[2:end,t] + vp_gelu5day[2:end,t], label="Online NN closure, 5 day")
    lines!(ax, lr_wl[2:end], up_20day[2:end,t] + vp_20day[2:end,t], label="Online NN closure, 20 day")
    # lines!(fig[1,1], lr_wl[2:end], up_geluKEspecpd[2:end,t] + up_geluKEspecpd[2:end,t], label="Online NN closure, KE spectrum PD")
    lines!(ax, lr_wl[2:end], up_noparam[2:end,t] + vp_noparam[2:end,t], label="30 km resolution, no closure")
    lines!(ax, lr_wl[2:end], up_zb[2:end,t] + vp_zb[2:end,t], label="ZB closure")
    lines!(ax, lr_wl[2:end], up_multi1[2:end,t] + vp_multi1[2:end,t], label="Online NN closure, batched 1 day")
    lines!(ax, lr_wl[2:end], up_multi2[2:end,t] + vp_multi2[2:end,t], label="Online NN closure, batched 2 day")
    lines!(ax, lr_wl[2:end], up_multi10[2:end,t] + vp_multi10[2:end,t], label="Online NN closure, multi 10 day")
    lines!(ax, lr_wl[2:end], up_multi20[2:end,t] + vp_multi20[2:end,t], label="Online NN closure, multi 20 day")
    lines!(ax, lr_wl[2:end], up_multi3_old[2:end,t] + vp_multi3_old[2:end,t], label="Online NN closure, multi 3 day old")
    lines!(ax, lr_wl[2:end], up_multi3_new[2:end,t] + vp_multi3_new[2:end,t], label="Online NN closure, multi 3 day new")

    axislegend(position = (0,0))

    #### comparing time-averaged ke spectra

    up_noparam_avg = zeros(65)
    vp_noparam_avg = zeros(65)

    up_zb_avg = zeros(65)
    vp_zb_avg = zeros(65)

    up_10day_avg = zeros(65)
    vp_10day_avg = zeros(65)

    up_20day_avg = zeros(65)
    vp_20day_avg = zeros(65)

    up_30day_avg = zeros(65)
    vp_30day_avg = zeros(65)

    up_hr_avg = zeros(513)
    vp_hr_avg = zeros(513)

    up_cghr_avg = zeros(65)
    vp_cghr_avg = zeros(65)

    up_filter_avg = zeros(513)
    vp_filter_avg = zeros(513)

    up_gelu1day_avg = zeros(65)
    vp_gelu1day_avg = zeros(65)

    up_5day_avg = zeros(65)
    vp_5day_avg = zeros(65)

    up_multi10_avg = zeros(65)
    vp_multi10_avg = zeros(65)

    up_multi1_avg = zeros(65)
    vp_multi1_avg = zeros(65)

    up_multi20_avg = zeros(65)
    vp_multi20_avg = zeros(65)

    up_multi2_avg = zeros(65)
    vp_multi2_avg = zeros(65)

    up_multi3_avg = zeros(65)
    vp_multi3_avg = zeros(65)

    up_multi3more_avg = zeros(65)
    vp_multi3more_avg = zeros(65)

    # up_geluKEspec_avg = zeros(65)
    # vp_geluKEspec_avg = zeros(65)

    # up_geluKEspecpd_avg = zeros(65)
    # vp_geluKEspecpd_avg = zeros(65)

    # up_geluhybrid_avg = zeros(65)
    # vp_geluhybrid_avg = zeros(65)

    # up_gelufourier_avg = zeros(65)
    # vp_gelufourier_avg = zeros(65)

    # up_relu1day_avg = zeros(65)
    # vp_relu1day_avg = zeros(65)

    # up_relu5day_avg = zeros(65)
    # vp_relu5day_avg = zeros(65)

    # up_reluKEspec_avg = zeros(65)
    # vp_reluKEspec_avg = zeros(65)

    totalstates = 1096
    for t = 1:totalstates

        up_noparam_avg += up_noparam[:,t]
        vp_noparam_avg += vp_noparam[:,t]

        up_zb_avg += up_zb[:,t]
        vp_zb_avg += vp_zb[:,t]

        up_20day_avg += up_20day[:, t]
        vp_20day_avg += vp_20day[:, t]

        up_30day_avg += up_30day[:, t]
        vp_30day_avg += vp_30day[:, t]

        # up_hr_avg += up_hr[:,t]
        # vp_hr_avg += vp_hr[:,t]

        up_filter_avg += up_hrfilter[:,t]
        vp_filter_avg += vp_hrfilter[:,t]

        up_cghr_avg += up_hrcg[:,t]
        vp_cghr_avg += vp_hrcg[:,t]

        up_gelu1day_avg += up_gelu1day[:,t]
        vp_gelu1day_avg += vp_gelu1day[:,t]

        up_5day_avg += up_gelu5day[:,t]
        vp_5day_avg += vp_gelu5day[:,t]

        up_10day_avg += up_10day[:, t]
        vp_10day_avg += vp_10day[:, t]

        up_multi1_avg += up_multi1[:, t]
        vp_multi1_avg += vp_multi1[:, t]

        up_multi2_avg += up_multi2[:, t]
        vp_multi2_avg += vp_multi2[:, t]

        up_multi10_avg += up_multi10[:, t]
        vp_multi10_avg += vp_multi10[:, t]

        up_multi20_avg += up_multi20[:, t]
        vp_multi20_avg += vp_multi20[:, t]

        up_multi3_avg += up_multi3[:, t]
        vp_multi3_avg += vp_multi3[:, t]

        up_multi3more_avg += up_multi3[:, t]
        vp_multi3more_avg += vp_multi3[:, t]

        # up_geluKEspec_avg += up_geluKEspec[:,t]
        # vp_geluKEspec_avg += vp_geluKEspec[:,t]

        # up_geluKEspecpd_avg += up_geluKEspecpd[:,t]
        # vp_geluKEspecpd_avg += vp_geluKEspecpd[:,t]

        # up_geluhybrid_avg += up_geluhybrid[:,t]
        # vp_geluhybrid_avg += vp_geluhybrid[:,t]

        # up_gelufourier_avg += up_gelufourier[:,t]
        # vp_gelufourier_avg += vp_gelufourier[:,t]

    end

    fig = Figure(size=(900, 400), fontsize=15);
    ax = Axis(fig[1,1],
        xscale=log10,
        yscale=log10,
        xlabel="Wavelength (km)",
        ylabel=L"\text{KE(k) } (m^3/s^2)", xreversed=true,
        xticks=[700, 400, 100, 30, 10, 2],
        title="3-year averaged KE spectrum"
    )
    lines!(ax, lr_wl[2:end], (up_cghr_avg[2:end] + vp_cghr_avg[2:end])/totalstates, label="Filtered, coarse-grained 3.75 km", color=:black)
    lines!(ax, lr_wl[2:end], (up_zb_avg[2:end] + vp_zb_avg[2:end])/totalstates, label="ZB20", color=:red)
    lines!(ax, lr_wl[2:end], (up_noparam_avg[2:end] + vp_noparam_avg[2:end])/totalstates, label="No closure, 30 km", color=:gray)

    # lines!(ax, lr_wl[2:end], (up_multi1_avg[2:end] + vp_multi1_avg[2:end])/1096, label="Online NN closure, batched 1 day")#, linestyle=:dash)
    lines!(ax, lr_wl[2:end], (up_multi2_avg[2:end] + vp_multi2_avg[2:end])/totalstates, label="Ensemble 2 day")#, linestyle=:dash)
    lines!(ax, lr_wl[2:end], (up_multi3more_avg[2:end] + vp_multi3more_avg[2:end])/totalstates, label="Ensemble 3 day")#, linestyle=:dashdot)
    lines!(ax, lr_wl[2:end], (up_multi10_avg[2:end] + vp_multi10_avg[2:end])/totalstates, label="Ensemble 10 day")#,linestyle=:dot)
    axislegend(position = :lb)


    ax2 = Axis(fig[1,1],
        xscale=log10,
        yscale=log10,
        xlabel="Wavelength (km)",
        ylabel=L"\text{KE(k) } (m^3/s^2)", xreversed=true,
        xticks=[700, 400, 100, 30, 10, 2],
        title="3-year averaged kinetic energy spectrum"
    )
    lines!(ax2, lr_wl[2:end], (up_cghr_avg[2:end] + vp_cghr_avg[2:end])/totalstates, label="Filtered, coarse-grained 3.75 km", color=:black)
    lines!(ax2, lr_wl[2:end], (up_zb_avg[2:end] + vp_zb_avg[2:end])/totalstates, label="ZB20", color=:red)
    lines!(ax2, lr_wl[2:end], (up_noparam_avg[2:end] + vp_noparam_avg[2:end])/totalstates, label="No closure, 30 km", color=:gray)
    lines!(ax2, lr_wl[2:end], (up_5day_avg[2:end] + vp_5day_avg[2:end])/1096, label="5 day")
    lines!(ax2, lr_wl[2:end], (up_10day_avg[2:end] + vp_10day_avg[2:end])/1096, label="10 day")
    lines!(ax2, lr_wl[2:end], (up_20day_avg[2:end] + vp_20day_avg[2:end])/1096, label="20 day")
    lines!(ax2, lr_wl[2:end], (up_30day_avg[2:end] + vp_30day_avg[2:end])/1096, label="30 day")#, linestyle=:dashdot)
    lines!(ax, lr_wl[2:end], (up_multi20_avg[2:end] + vp_multi20_avg[2:end])/totalstates, label="Ensemble 20 day")#, color=:red)

    # Legend(fig[1, 2], ax2, orientation = :vertical)


    # lines!(fig[1,1], lr_wl[2:end], (up_5dayeta_avg[2:end] + vp_5dayeta_avg[2:end])/1098, label="Online NN closure, 5 day with eta")
    # lines!(fig[1,1], lr_wl[2:end], (up_geluKEspec_avg[2:end] + vp_geluKEspec_avg[2:end])/31, label="Online NN closure, KE spec")
    # lines!(fig[1,1], lr_wl[2:end], (up_geluKEspecpd_avg[2:end] + vp_geluKEspecpd_avg[2:end])/31, label="Online NN closure, KE spec percent-diff")
    # lines!(fig[1,1], lr_wl[2:end], (up_geluhybrid_avg[2:end] + vp_geluhybrid_avg[2:end])/31, label="Online NN closure, Hybrid")
    # lines!(fig[1,1], lr_wl[2:end], (up_gelufourier_avg[2:end] + vp_gelufourier_avg[2:end])/31, label="Online NN closure, Fourier")
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
