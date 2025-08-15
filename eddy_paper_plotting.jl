"""
Mostly figure generation, I just wanted to be able to run include("technical_paper.jl")
without all of this also running. The function deserialize opens any saved checkpoints/deprecated now
"""
function deserialize(x)
    s = IOBuffer(x)
    Serialization.deserialize(s)
end

function load_and_create_models()

    offlineweights

    # load once
    u_zb = ncread("./spinup_files/128_zbforcingdissipation_cginitcond_postspinup_365days_071825/u.nc", "u")
    v_zb = ncread("./spinup_files/128_zbforcingdissipation_cginitcond_postspinup_365days_071825/v.nc", "v")
    eta_zb = ncread("./spinup_files/128_zbforcingdissipation_cginitcond_postspinup_365days_071825/eta.nc", "eta")

    uhr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/u.nc", "u")
    vhr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/v.nc", "v")
    etahr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/eta.nc", "eta")

    u0 = load_object("./spinup_files/coarsegrained_1024_10yearstate_061925.jld2")[1]
    v0 = load_object("./spinup_files/coarsegrained_1024_10yearstate_061925.jld2")[2]
    eta0 = load_object("./spinup_files/coarsegrained_1024_10yearstate_061925.jld2")[3]

    Ndays = 365
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
    Snoparam.Prog.u .= copy(initial_cond[1])
    Snoparam.Prog.v .= copy(initial_cond[2])
    Snoparam.Prog.η .= copy(initial_cond[3])

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
    Snn.Prog.u .= copy(initial_cond[1])
    Snn.Prog.v .= copy(initial_cond[2])
    Snn.Prog.η .= copy(initial_cond[3])

    param_guess = result.minimizer
    # param_guess = load_object("./initialweights_standarddeviation1_justrandomnumbers.jld2")
    current = 1
    for model in (Snn.Diag.NNVars.model_diag, Snn.Diag.NNVars.model_offdiag)
        for layers in model[1]
            for array in layers
                sz = prod(size(array))
                array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
                current += sz
            end
        end
    end

    Strainednn_states = ShallowWaters.model_setup(output=false,
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
    Strainednn_states.Prog.u .= copy(initial_cond[1])
    Strainednn_states.Prog.v .= copy(initial_cond[2])
    Strainednn_states.Prog.η .= copy(initial_cond[3])

    # param_guess = result.minimizer
    param_guess = load_object("./tuned_weights/multistate_dailydata_1:2:10daysintegration_result_071725.jld2").minimizer
    current = 1
    for model in (Strainednn_states.Diag.NNVars.model_diag, Strainednn_states.Diag.NNVars.model_offdiag)
        for layers in model[1]
            for array in layers
                sz = prod(size(array))
                array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
                current += sz
            end
        end
    end

    Strainednn_kespec = ShallowWaters.model_setup(output=false,
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
    Strainednn_kespec.Prog.u .= copy(initial_cond[1])
    Strainednn_kespec.Prog.v .= copy(initial_cond[2])
    Strainednn.Prog.η .= copy(initial_cond[3])

    # param_guess = result.minimizer
    param_guess = load_object("./tuned_weights/result_kespec_tendays_dailydata_imagefiltering_5iterationsLBFGS_072225.jld2").minimizer
    current = 1
    for model in (Strainednn_kespec.Diag.NNVars.model_diag, Strainednn_kespec.Diag.NNVars.model_offdiag)
        for layers in model[1]
            for array in layers
                sz = prod(size(array))
                array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
                current += sz
            end
        end
    end

    Strainednn_pd = ShallowWaters.model_setup(output=false,
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
    Strainednn_pd.Prog.u .= copy(initial_cond[1])
    Strainednn_pd.Prog.v .= copy(initial_cond[2])
    Strainednn_pd.Prog.η .= copy(initial_cond[3])

    # param_guess = result.minimizer
    param_guess = load_object("./tuned_weights/result_percentdiff_tendays_dailydata_imagefiltering_10iterationsLBFGS_072225.jld2").minimizer
    current = 1
    for model in (Strainednn_pd.Diag.NNVars.model_diag, Strainednn_pd.Diag.NNVars.model_offdiag)
        for layers in model[1]
            for array in layers
                sz = prod(size(array))
                array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
                current += sz
            end
        end
    end

    Strainednn_pd65 = ShallowWaters.model_setup(output=false,
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
    Strainednn_pd65.Prog.u .= copy(initial_cond[1])
    Strainednn_pd65.Prog.v .= copy(initial_cond[2])
    Strainednn_pd65.Prog.η .= copy(initial_cond[3])

    # param_guess = result.minimizer
    param_guess = load_object("./tuned_weights/result_percentdiff_imagefiltering_dailydata_truncatingat65_10days_10iterationsLBFGS_072225.jld2").minimizer
    current = 1
    for model in (Strainednn_pd65.Diag.NNVars.model_diag, Strainednn_pd65.Diag.NNVars.model_offdiag)
        for layers in model[1]
            for array in layers
                sz = prod(size(array))
                array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
                current += sz
            end
        end
    end

    # states_nn = save_states(Snn) # untrained
    # states_noparam = save_states(Snoparam)
    # states_trainednn_kespec = save_states(Strainednn_kespec)
    # states_trainednn_states = save_states(Strainednn_states)
    # states_trainednn_pd = save_states(Strainednn_pd)
    # states_trainednn_pd65 = save_states(Strainednn_pd65)

    states_nn = load_object("./spinup_files/nnresult_nottrained_oneyearintegration_dailysaves_072325.jld2")
    states_noparam = load_object("./spinup_files/output_noparam_oneyearintegration_dailysaves_072325.jld2")
    states_trainednn_kespec = load_object("./spinup_files/nnresult_trained__kespec_oneyearintegration_dailysaves_072325.jld2")
    states_trainednn_states = load_object("./spinup_files/nnresult_trained_stateloss_oneyearintegration_dailysaves_072325.jld2")
    states_trainednn_pd = load_object("./spinup_files/nnresult_trained_kespecpd_oneyearintegration_dailysaves_072325.jld2")

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
    ker = ImageFiltering.Kernel.gaussian((30e3/3750))
    # imfilter(hru[:,:,j], reflect(ker))

    t = 366
    fig = Figure(size=(800, 400), fontsize=15);
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 1024),
    LinRange(0, 3840, 1024),
    (uhr[:,1:end-1,t].^2 .+ vhr[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="E"),
    colorrange=(0,
    maximum((uhr[:,1:end-1,t].^2 .+ vhr[1:end-1,:,t].^2))
    ));
    Colorbar(fig[1,2], hm1)

    uhrcg = imfilter(uhr[:,:,t], reflect(ker))
    vhrcg = imfilter(vhr[:,:,t], reflect(ker))
    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 1024),
    LinRange(0, 3840, 1024),
    (uhrcg[:,1:end-1].^2 .+ vhrcg[1:end-1,:].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Coarse-grained E"),
    colorrange=(0,
    maximum((uhr[:,1:end-1,t].^2 .+ vhr[1:end-1,:,t].^2))
    ));
    Colorbar(fig[1,4], hm1)

    # cg, zb, nn, no param
    fig = Figure(size=(900, 1000), fontsize=15);
    t = 31
    uhrcg = imfilter(uhr[:,:,t], reflect(ker))
    vhrcg = imfilter(vhr[:,:,t], reflect(ker))
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 1024),
    LinRange(0, 3840, 1024),
    (uhrcg[:,1:end-1].^2 .+ vhrcg[1:end-1,:].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Coarse-grained E"),
    colorrange=(0,
    maximum((uhrcg[:,1:end-1].^2 .+ vhrcg[1:end-1,:].^2))
    ));
    Colorbar(fig[1,2], hm1)

    ax1, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (states_noparam[30].u[:,1:end-1].^2 .+ states_noparam[30].v[1:end-1,:].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30km resolution E, no closure"),
    colorrange=(0,
    maximum((uhrcg[:,1:end-1].^2 .+ vhrcg[1:end-1,:].^2))
    ));
    Colorbar(fig[1,4], hm1)

    ax1, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (u_zb[:,1:end-1,t].^2 .+ v_zb[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30km resolution E with ZB closure"),
    colorrange=(0,
    maximum((uhrcg[:,1:end-1].^2 .+ vhrcg[1:end-1,:].^2))
    ));
    Colorbar(fig[2,2], hm1)

    ax1, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (states_nn[t-1].u[:,1:end-1].^2 .+ states_nn[t-1].v[1:end-1,:].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30km resolution E with untrained NN closure"),
    colorrange=(0,
    maximum((uhrcg[:,1:end-1].^2 .+ vhrcg[1:end-1,:].^2))
    ));
    Colorbar(fig[2,4], hm1)

    ax1, hm5 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (states_offline[t-1].u[:,1:end-1].^2 .+ states_offline[t-1].v[1:end-1,:].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30km resolution E with ''offline`` NN closure"),
    colorrange=(0,
    maximum((uhrcg[:,1:end-1].^2 .+ vhrcg[1:end-1,:].^2))
    ));
    Colorbar(fig[3,2], hm1)

    ax1, hm6 = heatmap(fig[3,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (states_trainednn_states[t-1].u[:,1:end-1].^2 .+ states_trainednn_states[t-1].v[1:end-1,:].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30km resolution E with state tuned NN"),
    colorrange=(0,
    maximum((uhrcg[:,1:end-1].^2 .+ vhrcg[1:end-1,:].^2))
    ));
    Colorbar(fig[3,4], hm1)

    # comparing trained NN results

    ###################################################################################

    ## Different NN results, energy ##########################################################

    # nc files
    # u_zb, uhr
    # jld2 files (my save states function)
    # states_noparam, states_nn (untrained), states_trainednn_pd, states_trainednn_pd65, states_trainednn_kespec, states_trainednn_states

    # to get coarse-grained states
    ker = ImageFiltering.Kernel.gaussian((30e3/3750))
    # imfilter(hru[:,:,j], reflect(ker))

    fig = Figure(size=(900, 900), fontsize=15);
    t = 31
    ax1, hm = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (states_trainednn_states[t].u[:,1:end-1].^2 .+ states_trainednn_states[t].v[1:end-1,:].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Trained closure, states"),
    colorrange=(0,
    maximum((uhrcg[:,1:end-1].^2 .+ vhrcg[1:end-1,:].^2))
    ));
    Colorbar(fig[1,2], hm1)

    ax1, hm = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (states_trainednn_kespec[t].u[:,1:end-1].^2 .+ states_trainednn_kespec[t].v[1:end-1,:].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Trained closure, KE spectrum"),
    colorrange=(0,
    maximum((uhrcg[:,1:end-1].^2 .+ vhrcg[1:end-1,:].^2))
    ));
    Colorbar(fig[1,4], hm1)

    ax1, hm = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (states_trainednn_pd[t].u[:,1:end-1].^2 .+ states_trainednn_pd[t].v[1:end-1,:].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Trained closure, KE spectrum percent-diff"),
    colorrange=(0,
    maximum((uhrcg[:,1:end-1].^2 .+ vhrcg[1:end-1,:].^2))
    ));
    Colorbar(fig[2,2], hm1)

    ax1, hm = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (states_trainednn_pd65[t].u[:,1:end-1].^2 .+ states_trainednn_pd65[t].v[1:end-1,:].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Trained closure, KE spectrum percent-diff 65"),
    colorrange=(0,
    maximum((uhrcg[:,1:end-1].^2 .+ vhrcg[1:end-1,:].^2))
    ));
    Colorbar(fig[2,4], hm1)

    ##############################################################################################

    # comparing results, KE spectrum #############################################################

    # nc files
    # u_zb, uhr
    # jld2 files (my save states function)
    # states_noparam, states_nn (untrained), states_trainednn_pd, states_trainednn_pd65, states_trainednn_kespec, states_trainednn_states

    # to get coarse-grained states
    ker = ImageFiltering.Kernel.gaussian((30e3/3750))
    # imfilter(hru[:,:,j], reflect(ker))

    totalstates = 364 # saved hourly
    up_noparam = zeros(65,totalstates)
    vp_noparam = zeros(65,totalstates)

    up_zb = zeros(65,366)
    vp_zb = zeros(65,366)

    up_hr = zeros(513,366)
    vp_hr = zeros(513,366)

    up_cghr = zeros(513,366)
    vp_cghr = zeros(513,366)

    up_nn = zeros(65,totalstates)
    vp_nn = zeros(65,totalstates)

    up_nnkesp = zeros(65,totalstates)
    vp_nnkesp = zeros(65,totalstates)

    up_nnstates = zeros(65,totalstates)
    vp_nnstates = zeros(65,totalstates)

    up_nnpd = zeros(65,totalstates)
    vp_nnpd = zeros(65,totalstates)

    up_offline = zeros(65, totalstates)
    vp_offline = zeros(65, totalstates)

    for t = 1:366

        up_zb[:,t] = power(periodogram(u_zb[:, :, t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_zb[:,t] = power(periodogram(v_zb[:, :, t]; radialavg=true, radialsum=false)) ./ 128^2

        up_hr[:,t] = power(periodogram(uhr[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2
        vp_hr[:,t] = power(periodogram(vhr[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2

        up_cghr[:,t] = power(periodogram(imfilter(uhr[:,:,t], reflect(ker)); radialavg=true, radialsum=false)) ./ 1024^2
        vp_cghr[:,t] = power(periodogram(imfilter(vhr[:,:,t], reflect(ker)); radialavg=true, radialsum=false)) ./ 1024^2

    end

    for t = 1:totalstates

        up_noparam[:,t] = power(periodogram(states_noparam[t].u; radialavg=true, radialsum=false)) ./ 128^2
        vp_noparam[:,t] = power(periodogram(states_noparam[t].v; radialavg=true, radialsum=false)) ./ 128^2

        #untrained
        up_nn[:,t] = power(periodogram(states_nn[t].u; radialavg=true, radialsum=false)) ./ 128^2
        vp_nn[:,t] = power(periodogram(states_nn[t].v; radialavg=true, radialsum=false)) ./ 128^2

        # all trained
        up_nnstates[:,t] = power(periodogram(states_trainednn_states[t].u; radialavg=true, radialsum=false)) ./ 128^2
        vp_nnstates[:,t] = power(periodogram(states_trainednn_states[t].v; radialavg=true, radialsum=false)) ./ 128^2

        up_nnkesp[:,t] = power(periodogram(states_trainednn_kespec[t].u; radialavg=true, radialsum=false)) ./ 128^2
        vp_nnkesp[:,t] = power(periodogram(states_trainednn_kespec[t].v; radialavg=true, radialsum=false)) ./ 128^2

        up_nnpd[:,t] = power(periodogram(states_trainednn_pd[t].u; radialavg=true, radialsum=false)) ./ 128^2
        vp_nnpd[:,t] = power(periodogram(states_trainednn_pd[t].v; radialavg=true, radialsum=false)) ./ 128^2

        up_offline[:,t] = power(periodogram(states_offline[t].u; radialavg=true, radialsum=false)) ./ 128^2
        vp_offline[:,t] = power(periodogram(states_offline[t].v; radialavg=true, radialsum=false)) ./ 128^2

        # up_nnpd65[:,t] = power(periodogram(states_trainednn_pd65[t].u; radialavg=true, radialsum=false)) ./ 128^2
        # vp_nnpd65[:,t] = power(periodogram(states_trainednn_pd65[t].v; radialavg=true, radialsum=false)) ./ 128^2
    end

    lr_wl = (1 ./ freq(periodogram(u_zb[:,:,10]; radialavg=true, radialsum=false))) * 30;
    nnu_freq = LinRange(0, 64, 65)
    nnu_freq = nnu_freq ./ 65
    nnu_freq = 1 ./ nnu_freq 
    nnu_freq[1] =  1000
    lr_wl[1] = 1000

    hr_wl = 1 ./ freq(periodogram(uhr[:,:,3]; radialavg=true)) * 3.75;
    hr_wl[1] = 1100

    fig = Figure(size=(1000, 500), fontsize=15);
    t = 90
    lines(fig[1,1], hr_wl[2:65], up_cghr[2:65,t] + vp_cghr[2:65,t], label="Coarse-grained 3.75km resolution", axis=(
            xscale=log10,
            yscale=log10,
            xlabel="Wavelength (km)",
            ylabel="KE(k)",
            xreversed=true,
            xticks=[700, 100, 30, 10, 2],
            title="KE spectrum, 30 day integration")
    )
    lines!(fig[1,1], lr_wl[2:end], up_nn[2:end,t] + vp_nn[2:end,t], label="Untrained NN closure"
    )
    lines!(fig[1,1], lr_wl[2:end], up_noparam[2:end,t] + vp_noparam[2:end,t], label="30km resolution, no closure")
    lines!(fig[1,1], lr_wl[2:end], up_zb[2:end,t] + vp_zb[2:end,t], label="ZB closure")
    lines!(fig[1,1], lr_wl[2:end], up_nnstates[2:end,t] + vp_nnstates[2:end,t], label="NN closure, state loss")
    lines!(fig[1,1], lr_wl[2:end], up_nnkesp[2:end,t] + vp_nnkesp[2:end,t], label="NN closure, KE spectra loss")
    lines!(fig[1,1], lr_wl[2:end], up_nnpd[2:end,t] + vp_nnpd[2:end,t], label="NN closure, KE spectra percent diff loss")
    lines!(fig[1,1], lr_wl[2:end], up_offline[2:end,t] + vp_offline[2:end,t], label="NN closure, offline trained")
    axislegend(position = (0,0))

    #### comparing time-averaged ke spectra

    up_noparam_avg = zeros(65)
    vp_noparam_avg = zeros(65)

    up_zb_avg = zeros(65)
    vp_zb_avg = zeros(65)

    up_hr_avg = zeros(513)
    vp_hr_avg = zeros(513)

    up_cghr_avg = zeros(513)
    vp_cghr_avg = zeros(513)

    up_nn_avg = zeros(65)
    vp_nn_avg = zeros(65)

    up_nnkesp_avg = zeros(65)
    vp_nnkesp_avg = zeros(65)

    up_nnstates_avg = zeros(65)
    vp_nnstates_avg = zeros(65)

    up_nnpd_avg = zeros(65)
    vp_nnpd_avg = zeros(65)

    up_offline_avg = zeros(65)
    vp_offline_avg = zeros(65)

    daystoaverage = 364
    for t = 1:daystoaverage

        up_noparam_avg += up_noparam[:,t]
        vp_noparam_avg += vp_noparam[:,t]

        # untrained
        up_nn_avg += up_nn[:,t]
        vp_nn_avg += vp_nn[:,t]

        up_nnkesp_avg += up_nnkesp[:,t]
        vp_nnkesp_avg += vp_nnkesp[:,t]

        up_nnstates_avg += up_nnstates[:,t]
        vp_nnstates_avg += vp_nnstates[:,t]

        up_nnpd_avg += up_nnpd[:,t]
        vp_nnpd_avg += vp_nnpd[:,t]

        up_offline_avg += up_offline[:,t]
        vp_offline_avg += vp_offline[:,t]

    end

    for t = 1:(daystoaverage+1)

        up_zb_avg += up_zb[:,t]
        vp_zb_avg += vp_zb[:,t]

        up_hr_avg += up_hr[:,t]
        vp_hr_avg += vp_hr[:,t]

        up_cghr_avg += up_cghr[:,t]
        vp_cghr_avg += vp_cghr[:,t]

    end

    fig = Figure(size=(1000, 500), fontsize=15);
    lines(fig[1,1], hr_wl[2:65], (up_cghr_avg[2:65] + vp_cghr_avg[2:65])/31, label="Coarse-grained 3.75km resolution", axis=(
        xscale=log10,
        yscale=log10,
        xlabel="Wavelength (km)",
        ylabel="KE(k)", xreversed=true,
        xticks=[700, 100, 30, 10, 2],
        title="One-year averaged KE spectrum")
    )
    lines!(fig[1,1], lr_wl[2:end], (up_zb_avg[2:end] + vp_zb_avg[2:end])/31, label="ZB closure")
    lines!(fig[1,1], lr_wl[2:end], (up_noparam_avg[2:end] + vp_noparam_avg[2:end])/totalstates, label="30km resolution, no closure")
    lines!(fig[1,1], lr_wl[2:end], (up_nn_avg[2:end] + vp_nn_avg[2:end])/totalstates, label="Untrained NN closure")
    # lines!(fig[1,1], lr_wl[2:end], (up_nnstates_avg[2:end] + vp_nnstates_avg[2:end])/totalstates, label="NN closure, state loss")
    # lines!(fig[1,1], lr_wl[2:end], (up_nnkesp_avg[2:end] + vp_nnkesp_avg[2:end])/totalstates, label="NN closure, KE spectra loss")
    # lines!(fig[1,1], lr_wl[2:end], (up_nnpd_avg[2:end] + vp_nnpd_avg[2:end])/totalstates, label="NN closure, KE spectra percent diff loss")
    lines!(fig[1,1], lr_wl[2:end], (up_offline_avg[2:end] + vp_offline_avg[2:end])/totalstates, label="NN closure, offline loss")
    axislegend(position = (0,0))

    #############################################################################################


    # comparing cg energy spectra to energy spectra of the cg

    u_hr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/u.nc", "u")
    v_hr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/v.nc", "v")
    eta_hr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/eta.nc", "eta")

    ucg2_hr = load_object("./spinup_files/coarsegrained_1024_10yearstate_061925.jld2")[1]
    vcg2_hr = load_object("./spinup_files/coarsegrained_1024_10yearstate_061925.jld2")[2]
    etacg2_hr = load_object("./spinup_files/coarsegrained_1024_10yearstate_061925.jld2")[3]

    Slr = ShallowWaters.model_setup(output=false,
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
        Ndays=1
    )

    temp = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(
            ucg2_hr,
            vcg2_hr,
            etacg2_hr,
            Slr.Prog.sst,
            Slr
        )...)

    up_hr= zeros(513)
    ufp_hr= zeros(513)

    vp_hr = zeros(513)
    etap_hr = zeros(513)
    upcg2 = zeros(65)
    vpcg2 = zeros(65)
    etapcg2 = zeros(65)
    upcg = zeros(65)
    vpcg = zeros(65)
    etapcg = zeros(65)

    up_hr[:] = power(periodogram(u_hr[:,:,2]; radialavg=true, radialsum=false)) ./ 1024^2
    vp_hr[:] = power(periodogram(v_hr[:,:,1]; radialavg=true, radialsum=false)) ./ 1024^2
    etap_hr[:] = power(periodogram(eta_hr[:,:,1]; radialavg=true, radialsum=false)) ./ 1024^2

    ufp_hr[:] = power(periodogram(newu; radialavg=true, radialsum=false)) ./ 1024^2

    upcg[:] = power(periodogram(ucg; radialavg=true, radialsum=false)) ./ 128^2
    vpcg[:] = power(periodogram(vcg; radialavg=true, radialsum=false)) ./ 128^2
    etapcg[:] = power(periodogram(etacg; radialavg = true, radialsum=false)) ./128^2

    upcg2[:] = power(periodogram(temp.u; radialavg=true, radialsum=false)) ./ 128^2
    vpcg2[:] = power(periodogram(temp.v; radialavg=true, radialsum=false)) ./ 128^2
    etapcg2[:] = power(periodogram(temp.η; radialavg = true, radialsum=false)) ./128^2

    true_wl = 1 ./ freq(periodogram(u[:,:,1]; radialavg=true)) * 3.75;
    true_wl[1] = 1100
    cg_wl = (1 ./ freq(periodogram(temp.u; radialavg=true, radialsum=false))) * 30;
    cg_wl[1] = 1100

    fig = Figure();
    lines(fig[1,1], true_wl[2:end], up_hr[2:end] + vp_hr[2:end], label="HR", axis=(
            xscale=log10,yscale=log10,xlabel="Wavelength (km)", ylabel="KE(k)", xreversed=true, xticks=[700, 100, 30, 10, 2], title="HR energy spectrum")
    )
    lines!(fig[1,1], cg_wl[2:end], upcg[2:end] + vpcg[2:end], label="Coarse-grained HR spectrum")
    axislegend()

    figu = Figure();
    lines(figu[1,1], true_wl[2:65], up_hr[2:65], label="HR", axis=(
            xscale=log10,yscale=log10, ylabel="KE(k)", xreversed=true,  title="HR")
    )
    lines!(figu[1,1], cg_wl[2:end], upcg[2:end], label="Coarse-grained HR recomputed")
    lines!(figu[1,1], cg_wl[2:end], upcg2[2:end], label="Coarse-grained HR saved")
    axislegend()

    figv = Figure();
        lines(figv[1,1], true_wl[2:65], vp_hr[2:65], label="HR", axis=(
            xscale=log10,yscale=log10, ylabel="KE(k)", xreversed=true,  title="HR")
    )
    lines!(figv[1,1], cg_wl[2:end], vpcg[2:end], label="Coarse-grained HR")
    lines!(figv[1,1], cg_wl[2:end], vpcg2[2:end], label="Coarse-grained HR")
    axislegend()

    figeta = Figure();
        lines(figeta[1,1], true_wl[2:65], etap_hr[2:65], label="HR", axis=(
            xscale=log10,yscale=log10, ylabel="KE(k)", xreversed=true,  title="HR")
    )
    lines!(figeta[1,1], cg_wl[2:end], etapcg[2:end], label="Coarse-grained HR")
    lines!(figeta[1,1], cg_wl[2:end], etapcg2[2:end], label="Coarse-grained HR")
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

