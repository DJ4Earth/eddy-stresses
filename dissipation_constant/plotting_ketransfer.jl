
# using the definition 
#       T(k_x, k_y) = Re( F(u)^* F(S_x) + F(v)^* F(S_y) )
# where ^* is the complex conjugate. It's not clear to me if S should come from the same
# timestep or the prior, because the prior is what when into computing u and v
function computing_ketransfer_fromS()

    # computing the total subgrid forcing

    # duhrcg = load_object("./dissipation_constant/computing_trueS/hrcgtendencies_first3years_dailysaves_dudvdeta.jld2")[1];
    # dvhrcg = load_object("./dissipation_constant/computing_trueS/hrcgtendencies_first3years_dailysaves_dudvdeta.jld2")[2];

    # duhr = load_object("./dissipation_constant/computing_trueS/hrtendencies_first3years_dailysaves_dudvdeta.jld2")[1];
    # dvhr = load_object("./dissipation_constant/computing_trueS/hrtendencies_first3years_dailysaves_dudvdeta.jld2")[2];

    # advu_hr = load_object("./nonlinear_advec_fromhrstates_advu_advv.jld2")[1];
    # advv_hr = load_object("./nonlinear_advec_fromhrstates_advu_advv.jld2")[2];

    # advu_cg = load_object("./nonlinear_advec_fromcghrstates_advu_advv.jld2")[1];
    # advv_cg = load_object("./nonlinear_advec_fromcghrstates_advu_advv.jld2")[2];

    # duhrcg = load_object("./dissipation_constant/ke_transfers/tendencies_dwindowuhrcg_dwindowvhrcg_dwindowetahrcg_overline(window(state)).jld2")[1];
    # dvhrcg = load_object("./dissipation_constant/ke_transfers/tendencies_dwindowuhrcg_dwindowvhrcg_dwindowetahrcg_overline(window(state)).jld2")[2];

    # duhr = load_object("./dissipation_constant/ke_transfers/tendencies_dwindowuhr_dwindowvhr_dwindowetahr.jld2")[1];
    # dvhr = load_object("./dissipation_constant/ke_transfers/tendencies_dwindowuhr_dwindowvhr_dwindowetahr.jld2")[2];

    duhrcg = load_object("./dissipation_constant/ke_transfers/tendencies_dwindowuhrcg_dwindowvhrcg_nodeta_didnotapplywindowtoeta_overline(window(state)).jld2")[1];
    dvhrcg = load_object("./dissipation_constant/ke_transfers/tendencies_dwindowuhrcg_dwindowvhrcg_nodeta_didnotapplywindowtoeta_overline(window(state)).jld2")[2];

    duhr = load_object("./dissipation_constant/ke_transfers/tendencies_dwindowuhr_dwindowvhr_nodeta_didnotapplywindowtoeta.jld2")[1];
    dvhr = load_object("./dissipation_constant/ke_transfers/tendencies_dwindowuhr_dwindowvhr_nodeta_didnotapplywindowtoeta.jld2")[2];

    Su_true_noetawindow = zeros(127, 128, 1096)
    Sv_true_noetawindow = zeros(128, 127, 1096)

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    for j = 1:1096

        dufiltered = imfilter(@view(duhr[:,:,j]), reflect(ker))
        dvfiltered = imfilter(@view(dvhr[:,:,j]), reflect(ker))
    
        @views duhrdownsized = (dufiltered[8:8:end, 4:8:end] .+ dufiltered[8:8:end, 5:8:end]) .* 0.5
        @views dvhrdownsized = (dvfiltered[4:8:end, 8:8:end] .+ dvfiltered[5:8:end, 8:8:end]) .* 0.5

        # advufiltered = imfilter(@view(advu_hr[:,:,j]), reflect(ker))
        # advvfiltered = imfilter(@view(advv_hr[:,:,j]), reflect(ker))

        # advuhrdownsized = (advufiltered[8:8:end, 4:8:end] .+ advufiltered[8:8:end, 5:8:end]) .* 0.5
        # advvhrdownsized = (advvfiltered[4:8:end, 8:8:end] .+ advvfiltered[5:8:end, 8:8:end]) .* 0.5

        @views Su_true_noetawindow[:,:,j] .= -(duhrdownsized[:,:,1]./48) + (duhrcg[:,:,j]./384)
        @views Sv_true_noetawindow[:,:,j] .= -(dvhrdownsized[:,:,1]./48) + (dvhrcg[:,:,j]./384)

    end

    T = Float64
    Ponline = ShallowWaters.Parameter(T=T,
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
    );

    T = Float64
    PZB = ShallowWaters.Parameter(T=T,
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
        zb_forcing_dissipation=true,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=false,
        N=1,
        α=2,
        nx=128,
        Ndays=1
    );

    Shrcg = ShallowWaters.model_setup(Ponline);
    SZB = ShallowWaters.model_setup(PZB);

    S10 = ShallowWaters.model_setup(Ponline);
    current = 1
    for m in (S10.Diag.CNNVars.model_Su, S10.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_states_10dayoptimization_startfrom5daystate_noeta_gelu_30iterations.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    S20 = ShallowWaters.model_setup(Ponline);
    current = 1
    for m in (S20.Diag.CNNVars.model_Su, S20.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_states_20dayoptimization_startfrom10daystate_noeta_10iterations.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    S5 = ShallowWaters.model_setup(Ponline);

    S30 = ShallowWaters.model_setup(Ponline);
    current = 1
    for m in (S30.Diag.CNNVars.model_Su, S30.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_state_30dayoptimzation_startfrom20day_constantdissipation_6iterations_8hourdata_200maxhistory_fixedcfl.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    Smulti2 = ShallowWaters.model_setup(Ponline);
    current = 1
    for m in (Smulti2.Diag.CNNVars.model_Su, Smulti2.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-6-8-10-13-15-18-23-28-33-38-41-44-48-51-53-58-63-65-68-73-78-83-86-88daystart_2dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    Smulti3 = ShallowWaters.model_setup(Ponline);
    current = 1
    for m in (Smulti3.Diag.CNNVars.model_Su, Smulti3.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86daystart_3dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end


    Smulti10 = ShallowWaters.model_setup(Ponline);
    current = 1
    for m in (Smulti10.Diag.CNNVars.model_Su, Smulti10.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/result_multistate_5-20-35-50-65-75daystart_10dayoptimization_initialweights20daystate_fixedcfl_15iterations_constdissipation.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    Smulti20 = ShallowWaters.model_setup(Ponline);
    current = 1
    for m in (Smulti20.Diag.CNNVars.model_Su, Smulti20.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/result_multistate_5-25-45-65daystart_20dayoptimization_initialweights20daystate_fixedcfl_15iterations_constdissipation.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    true_S = load_object("./dissipation_constant/computing_trueS/trueS_fromtendencies_withRK4_SuSv_first3years_dailysaves_032526.jld2");
    Suhr = true_S[1];
    Svhr = true_S[2];

    approx_S = load_object("./dissipation_constant/computing_trueS/approxS_nonlinearadvec_nottendencies_hasextraDelta_SuSv_first3years_040726.jld2");
    Suapprox = approx_S[1];
    Svapprox = approx_S[2];

    # for computing a KE transfer from some/all of the domain
    xvals = :
    # xvals = Int((150/30)):Int((3840-150)/30)
    # xvals = Int((300/30)):Int((3840-300)/30)
    yvals = xvals

    lr_freq = 1/30 .* freq(periodogram(umulti2all[xvals,yvals,10]; radialavg=true, radialsum=false));
    nfft = nextfastfft(size(uhrcgall[xvals,yvals,1]))

    totalu_hrcg = zeros(length(lr_freq))
    totalv_hrcg = zeros(length(lr_freq))

    totalu_approx = zeros(length(lr_freq))
    totalv_approx = zeros(length(lr_freq))

    totalu_hrcg2 = zeros(length(lr_freq))
    totalv_hrcg2 = zeros(length(lr_freq))

    totalu_5 = zeros(length(lr_freq))
    totalv_5 = zeros(length(lr_freq))

    totalu_10 = zeros(length(lr_freq))
    totalv_10 = zeros(length(lr_freq))

    totalu_20 = zeros(length(lr_freq))
    totalv_20 = zeros(length(lr_freq))
    
    totalu_30 = zeros(length(lr_freq))
    totalv_30 = zeros(length(lr_freq))

    totalu_multi1 = zeros(length(lr_freq))
    totalv_multi1 = zeros(length(lr_freq))

    totalu_multi1more = zeros(length(lr_freq))
    totalv_multi1more = zeros(length(lr_freq))

    totalu_multi2 = zeros(length(lr_freq))
    totalv_multi2 = zeros(length(lr_freq))

    totalu_multi3 = zeros(length(lr_freq))
    totalv_multi3 = zeros(length(lr_freq))

    totalu_multi5 = zeros(length(lr_freq))
    totalv_multi5 = zeros(length(lr_freq))

    totalu_multi10 = zeros(length(lr_freq))
    totalv_multi10 = zeros(length(lr_freq))

    totalu_multi20 = zeros(length(lr_freq))
    totalv_multi20 = zeros(length(lr_freq))

    totalu_ZB = zeros(length(lr_freq))
    totalv_ZB = zeros(length(lr_freq))

    totalu_3 = zeros(length(lr_freq))
    totalv_3 = zeros(length(lr_freq))

    totalstates = 1096
    for t = 1:totalstates

        # from approximation to nonlinear advection
        outu_approx, inputu_approx, inputSu_approx = paddingu(uhrcgall[xvals, yvals, t], Suapprox[xvals,yvals,t], nfft[1])
        fft2pow2radial!(outu_approx, rfft(inputu_approx), rfft(inputSu_approx), nfft...)
        outv_approx, inputv_approx, inputSv_approx = paddingv(vhrcgall[xvals, yvals, t], Svapprox[xvals,yvals,t], nfft[1])
        fft2pow2radial!(outv_approx, rfft(inputv_approx), rfft(inputSv_approx), nfft...)

        totalu_approx += outu_approx
        totalv_approx += outv_approx

        # from total tendencies
        outu_hrcg, inputu_hrcg, inputSu_hrcg = paddingu(uhrcgall[xvals, yvals, t], Suhr[xvals,yvals,t], nfft[1])
        fft2pow2radial!(outu_hrcg, rfft(inputu_hrcg), rfft(inputSu_hrcg), nfft...)
        outv_hrcg, inputv_hrcg, inputSv_hrcg = paddingv(vhrcgall[xvals, yvals, t], Svhr[xvals,yvals,t], nfft[1])
        fft2pow2radial!(outv_hrcg, rfft(inputv_hrcg), rfft(inputSv_hrcg), nfft...)

        totalu_hrcg += outu_hrcg
        totalv_hrcg += outv_hrcg

        # ZB20 ############################

        uzb_, vzb_, _ = ShallowWaters.add_halo(Float64.(uzb[:,:,t]), Float64.(vzb[:,:,t]), Float64.(etazb[:,:,t]), zeros(128,128), SZB);
        ShallowWaters.ZB_momentum(uzb_, vzb_, SZB, SZB.Diag);

        outu_ZB, inputu_ZB, inputSu_ZB = paddingu(uzb[xvals, yvals, t], SZB.Diag.ZBVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_ZB, rfft(inputu_ZB), rfft(inputSu_ZB), nfft...)
        outv_ZB, inputv_ZB, inputSv_ZB = paddingv(vzb[xvals, yvals, t], SZB.Diag.ZBVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_ZB, rfft(inputv_ZB), rfft(inputSv_ZB), nfft...)

        totalu_ZB += outu_ZB
        totalv_ZB += outv_ZB

        # 5 day optimization ################

        # u5, v5, eta5 = ShallowWaters.add_halo(Float64.(u5day[:,:,t]), Float64.(v5day[:,:,t]), Float64.(eta5day[:,:,t]), zeros(128,128), S5);
        # ShallowWaters.CNN_momentum(u5, v5, S5);
        # outu_5, inputu_5, inputSu_5 = paddingu(u5day[xvals, yvals, t], S5.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_5, rfft(inputu_5), rfft(inputSu_5), nfft...)
        # outv_5, inputv_5, inputSv_5 = paddingv(v5day[xvals, yvals, t], S5.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_5, rfft(inputv_5), rfft(inputSv_5), nfft...)

        # totalu_5 += outu_5
        # totalv_5 += outv_5

        # 20 day optimization ##############

        u20, v20, _ = ShallowWaters.add_halo(Float64.(u20s[:,:,t]), Float64.(v20s[:,:,t]), Float64.(eta20s[:,:,t]), zeros(128,128), S20);
        ShallowWaters.CNN_momentum(u20, v20, S20)
        outu_20, inputu_20, inputSu_20 = paddingu(u20s[xvals, yvals, t], S20.Diag.CNNVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_20, rfft(inputu_20), rfft(inputSu_20), nfft...)
        outv_20, inputv_20, inputSv_20 = paddingv(v20s[xvals, yvals, t], S20.Diag.CNNVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_20, rfft(inputv_20), rfft(inputSv_20), nfft...)

        totalu_20 += outu_20
        totalv_20 += outv_20

        # 30 day optimization ###############

        u30, v30, _ = ShallowWaters.add_halo(Float64.(u30s[:,:,t]), Float64.(v30s[:,:,t]), Float64.(eta30s[:,:,t]), zeros(128,128), S30);
        ShallowWaters.CNN_momentum(u30, v30, S30)
        outu_30, inputu_30, inputSu_30 = paddingu(u30s[xvals, yvals, t], S30.Diag.CNNVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_30, rfft(inputu_30), rfft(inputSu_30), nfft...)
        outv_30, inputv_30, inputSv_30 = paddingv(v30s[xvals, yvals, t], S30.Diag.CNNVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_30, rfft(inputv_30), rfft(inputSv_30), nfft...)

        totalu_30 += outu_30
        totalv_30 += outv_30

        # batched 2 day ##################################

        umulti2_, vmulti2_, _ = ShallowWaters.add_halo(Float64.(umulti2[:,:,t]), Float64.(vmulti2[:,:,t]), Float64.(etamulti2[:,:,t]), zeros(128,128), Smulti2);
        ShallowWaters.CNN_momentum(umulti2_, vmulti2_, Smulti2)

        outu_multi2, inputu_multi2, inputSu_multi2 = paddingu(umulti2[xvals, yvals, t], Smulti2.Diag.CNNVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_multi2, rfft(inputu_multi2), rfft(inputSu_multi2), nfft...)
        outv_multi2, inputv_multi2, inputSv_multi2 = paddingv(vmulti2[xvals, yvals, t], Smulti2.Diag.CNNVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_multi2, rfft(inputv_multi2), rfft(inputSv_multi2), nfft...)

        totalu_multi2 += outu_multi2
        totalv_multi2 += outv_multi2

        # batched 3 day ###########################

        umulti3_, vmulti3_, _ = ShallowWaters.add_halo(Float64.(umulti3[:,:,t]), Float64.(vmulti3[:,:,t]), Float64.(etamulti3[:,:,t]), zeros(128,128), Smulti3);
        ShallowWaters.CNN_momentum(umulti3_, vmulti3_, Smulti3)
        outu_multi3, inputu_multi3, inputSu_multi3 = paddingu(umulti3[xvals, yvals, t], Smulti3.Diag.CNNVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_multi3, rfft(inputu_multi3), rfft(inputSu_multi3), nfft...)
        outv_multi3, inputv_multi3, inputSv_multi3 = paddingv(vmulti3[xvals, yvals, t], Smulti3.Diag.CNNVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_multi3, rfft(inputv_multi3), rfft(inputSv_multi3), nfft...)

        totalu_multi3 += outu_multi3
        totalv_multi3 += outv_multi3

        # batched 5 day ######################

        # umulti5_, vmulti5_, _ = ShallowWaters.add_halo(Float64.(umulti5[:,:,t]), Float64.(vmulti5[:,:,t]), Float64.(etamulti5[:,:,t]), zeros(128,128), Smulti5);
        # ShallowWaters.CNN_momentum(umulti5_, vmulti5_, Smulti5)
        # outu_multi5, inputu_multi5, inputSu_multi5 = paddingu(umulti5[xvals, yvals, t], Smulti5.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_multi5, rfft(inputu_multi5), rfft(inputSu_multi5), nfft...)
        # outv_multi5, inputv_multi5, inputSv_multi5 = paddingv(vmulti5[xvals, yvals, t], Smulti5.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_multi5, rfft(inputv_multi5), rfft(inputSv_multi5), nfft...)

        # totalu_multi5 += outu_multi5
        # totalv_multi5 += outv_multi5

        # batched 10 day #####################

        umulti10_, vmulti10_, _ = ShallowWaters.add_halo(Float64.(umulti10[:,:,t]), Float64.(vmulti10[:,:,t]), Float64.(etamulti10[:,:,t]), zeros(128,128), Smulti10);
        ShallowWaters.CNN_momentum(umulti10_, vmulti10_, Smulti10)
        outu_multi10, inputu_multi10, inputSu_multi10 = paddingu(umulti10[xvals, yvals, t], Smulti10.Diag.CNNVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_multi10, rfft(inputu_multi10), rfft(inputSu_multi10), nfft...)
        outv_multi10, inputv_multi10, inputSv_multi10 = paddingv(vmulti10[xvals, yvals, t], Smulti10.Diag.CNNVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_multi10, rfft(inputv_multi10), rfft(inputSv_multi10), nfft...)

        totalu_multi10 += outu_multi10
        totalv_multi10 += outv_multi10

        # batched 20 day #####################

        umulti20_, vmulti20_, _ = ShallowWaters.add_halo(Float64.(umulti20[:,:,t]), Float64.(vmulti20[:,:,t]), Float64.(etamulti20[:,:,t]), zeros(128,128), Smulti20);
        ShallowWaters.CNN_momentum(umulti20_, vmulti20_, Smulti20)
        outu_multi20, inputu_multi20, inputSu_multi20 = paddingu(umulti20[xvals, yvals, t], Smulti20.Diag.CNNVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_multi20, rfft(inputu_multi20), rfft(inputSu_multi20), nfft...)
        outv_multi20, inputv_multi20, inputSv_multi20 = paddingv(vmulti20[xvals, yvals, t], Smulti20.Diag.CNNVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_multi20, rfft(inputv_multi20), rfft(inputSv_multi20), nfft...)

        totalu_multi20 += outu_multi20
        totalv_multi20 += outv_multi20

    end

end

function computing_ketransfer_fromviscosity()
    # same as above but with the viscosity term
    visc_multi2 = load_object("./dissipation_constant/alternate_S_files/viscosity_ensemble2day_3years_dailysaves_MuMv.jld2");
    visc_multi3 = load_object("./dissipation_constant/alternate_S_files/viscosity_ensemble3day_3years_dailysaves_MuMv.jld2");
    visc_multi10 = load_object("./dissipation_constant/alternate_S_files/viscosity_ensemble10day_3years_dailysaves_MuMv.jld2");
    visc_ZB20 = load_object("./dissipation_constant/alternate_S_files/viscosity_ZB_3years_dailysaves_MuMv.jld2");

    visc_hrcg = load_object("./dissipation_constant/alternate_S_files/viscosity_cghr_3years_dailysaves_MuMv.jld2");

    T = Float64
    Ponline = ShallowWaters.Parameter(T=T,
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
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=1
    );

    T = Float64
    PZB = ShallowWaters.Parameter(T=T,
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
        zb_forcing_dissipation=true,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=false,
        N=1,
        α=2,
        nx=128,
        Ndays=1
    );

    Shrcg = ShallowWaters.model_setup(Ponline);
    SZB = ShallowWaters.model_setup(PZB);

    Smulti2 = ShallowWaters.model_setup(Ponline);
    Smulti3 = ShallowWaters.model_setup(Ponline);

    Smulti10 = ShallowWaters.model_setup(Ponline);

    lr_freq = 1/30 .* freq(periodogram(umulti1[:,:,10]; radialavg=true, radialsum=false));
    nfft = nextfastfft(size(uhrcgall[:,:,1]))

    viscu_hrcg = zeros(65)
    viscv_hrcg = zeros(65)

    viscu_hrcg2 = zeros(65)
    viscv_hrcg2 = zeros(65)

    viscu_5 = zeros(65)
    viscv_5 = zeros(65)

    viscu_10 = zeros(65)
    viscv_10 = zeros(65)

    viscu_20 = zeros(65)
    viscv_20 = zeros(65)
    
    viscu_30 = zeros(65)
    viscv_30 = zeros(65)

    viscu_multi1 = zeros(65)
    viscv_multi1 = zeros(65)

    viscu_multi1more = zeros(65)
    viscv_multi1more = zeros(65)

    viscu_multi2 = zeros(65)
    viscv_multi2 = zeros(65)

    viscu_multi3 = zeros(65)
    viscv_multi3 = zeros(65)

    viscu_multi5 = zeros(65)
    viscv_multi5 = zeros(65)

    viscu_multi10 = zeros(65)
    viscv_multi10 = zeros(65)

    viscu_multi20 = zeros(65)
    viscv_multi20 = zeros(65)

    viscu_ZB = zeros(65)
    viscv_ZB = zeros(65)

    viscu_3 = zeros(65)
    viscv_3 = zeros(65)

    totalstates = 1096

    for t = 1:totalstates

        # from total viscosity, computing as \overline{visc(u)} - visc(\overline{u})
        outu_hrcg, inputu_hrcg, inputSu_hrcg = paddingu(uhrcgall[:, :, t], visc_hrcg[1][xvals, yvals, t], nfft[1])
        fft2pow2radial!(outu_hrcg, rfft(inputu_hrcg), rfft(inputSu_hrcg), nfft...)
        outv_hrcg, inputv_hrcg, inputSv_hrcg = paddingv(vhrcgall[:, :, t], visc_hrcg[2][xvals, yvals, t], nfft[1])
        fft2pow2radial!(outv_hrcg, rfft(inputv_hrcg), rfft(inputSv_hrcg), nfft...)

        viscu_hrcg += outu_hrcg
        viscv_hrcg += outv_hrcg

        # ZB20 ############################

        outu_ZB, inputu_ZB, inputSu_ZB = paddingu(uzb[:, :, t], visc_ZB20[1][xvals, yvals, t], nfft[1])
        fft2pow2radial!(outu_ZB, rfft(inputu_ZB), rfft(inputSu_ZB), nfft...)
        outv_ZB, inputv_ZB, inputSv_ZB = paddingv(vzb[:, :, t], visc_ZB20[2][xvals, yvals, t], nfft[1])
        fft2pow2radial!(outv_ZB, rfft(inputv_ZB), rfft(inputSv_ZB), nfft...)

        viscu_ZB += outu_ZB
        viscv_ZB += outv_ZB

        # 30 day optimization ###############

        # u30, v30, _ = ShallowWaters.add_halo(Float64.(u30s[xvals, yvals, t]), Float64.(v30s[xvals, yvals, t]), Float64.(eta30s[xvals, yvals, t]), zeros(128,128), S30);
        # ShallowWaters.CNN_momentum(u30, v30, S30)
        # outu_30, inputu_30, inputSu_30 = paddingu(u30s[:, :, t], S30.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_30, rfft(inputu_30), rfft(inputSu_30), nfft...)
        # outv_30, inputv_30, inputSv_30 = paddingv(v30s[:, :, t], S30.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_30, rfft(inputv_30), rfft(inputSv_30), nfft...)

        # viscu_30 += outu_30
        # viscv_30 += outv_30

        # batched 2 day ##################################

        outu_multi2, inputu_multi2, inputSu_multi2 = paddingu(umulti2[:, :, t], visc_multi2[1][xvals, yvals, t], nfft[1])
        fft2pow2radial!(outu_multi2, rfft(inputu_multi2), rfft(inputSu_multi2), nfft...)
        outv_multi2, inputv_multi2, inputSv_multi2 = paddingv(vmulti2[:, :, t], visc_multi2[2][xvals, yvals, t], nfft[1])
        fft2pow2radial!(outv_multi2, rfft(inputv_multi2), rfft(inputSv_multi2), nfft...)

        viscu_multi2 += outu_multi2
        viscv_multi2 += outv_multi2

        # batched 3 day ###########################

        outu_multi3, inputu_multi3, inputSu_multi3 = paddingu(umulti3[:, :, t], visc_multi3[1][xvals, yvals, t], nfft[1])
        fft2pow2radial!(outu_multi3, rfft(inputu_multi3), rfft(inputSu_multi3), nfft...)
        outv_multi3, inputv_multi3, inputSv_multi3 = paddingv(vmulti3[:, :, t], visc_multi3[2][xvals, yvals, t], nfft[1])
        fft2pow2radial!(outv_multi3, rfft(inputv_multi3), rfft(inputSv_multi3), nfft...)

        viscu_multi3 += outu_multi3
        viscv_multi3 += outv_multi3

        # batched 10 day #####################

        outu_multi10, inputu_multi10, inputSu_multi10 = paddingu(umulti10[:, :, t], visc_multi10[1][xvals, yvals, t], nfft[1])
        fft2pow2radial!(outu_multi10, rfft(inputu_multi10), rfft(inputSu_multi10), nfft...)
        outv_multi10, inputv_multi10, inputSv_multi10 = paddingv(vmulti10[:, :, t], visc_multi10[2][:,:,t], nfft[1])
        fft2pow2radial!(outv_multi10, rfft(inputv_multi10), rfft(inputSv_multi10), nfft...)

        viscu_multi10 += outu_multi10
        viscv_multi10 += outv_multi10

    end

end

function computing_ketransfer_frombottomdrag()

    # Lastly, doing this with the bottom drag term
    bd_multi2 = load_object("./dissipation_constant/alternate_S_files/bottomdrag_ensemble2day_3years_dailysaves_BuBv.jld2");
    bd_multi3 = load_object("./dissipation_constant/alternate_S_files/bottomdrag_ensemble3day_3years_dailysaves_BuBv.jld2");
    bd_multi10 = load_object("./dissipation_constant/alternate_S_files/bottomdrag_ensemble10day_3years_dailysaves_BuBv.jld2");
    bd_ZB20 = load_object("./dissipation_constant/alternate_S_files/bottomdrag_ZB_3years_dailysaves_BuBv.jld2");

    bd_hrcg = load_object("./dissipation_constant/alternate_S_files/bottomdrag_cghr_3years_dailysaves_BuBv.jld2");

    T = Float64
    Ponline = ShallowWaters.Parameter(T=T,
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
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=1
    );

    T = Float64
    PZB = ShallowWaters.Parameter(T=T,
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
        zb_forcing_dissipation=true,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=false,
        N=1,
        α=2,
        nx=128,
        Ndays=1
    );

    Shrcg = ShallowWaters.model_setup(Ponline);
    SZB = ShallowWaters.model_setup(PZB);

    Smulti2 = ShallowWaters.model_setup(Ponline);
    Smulti3 = ShallowWaters.model_setup(Ponline);

    Smulti10 = ShallowWaters.model_setup(Ponline);

    lr_freq = 1/30 .* freq(periodogram(umulti2[:,:,10]; radialavg=true, radialsum=false));
    nfft = nextfastfft(size(uhrcgall[:,:,1]))

    bdu_hrcg = zeros(65)
    bdv_hrcg = zeros(65)

    bdu_multi2 = zeros(65)
    bdv_multi2 = zeros(65)

    bdu_multi3 = zeros(65)
    bdv_multi3 = zeros(65)

    bdu_multi10 = zeros(65)
    bdv_multi10 = zeros(65)

    bdu_ZB = zeros(65)
    bdv_ZB = zeros(65)

    totalstates = 1096

    for t = 1:totalstates

        # from coarse-grained high-resolution
        outu_hrcg, inputu_hrcg, inputSu_hrcg = paddingu(uhrcgall[:, :, t], bd_hrcg[1][:,:,t], nfft[1])
        fft2pow2radial!(outu_hrcg, rfft(inputu_hrcg), rfft(inputSu_hrcg), nfft...)
        outv_hrcg, inputv_hrcg, inputSv_hrcg = paddingv(vhrcgall[:, :, t], bd_hrcg[2][:,:,t], nfft[1])
        fft2pow2radial!(outv_hrcg, rfft(inputv_hrcg), rfft(inputSv_hrcg), nfft...)

        bdu_hrcg += outu_hrcg
        bdv_hrcg += outv_hrcg

        # ZB20 ############################

        outu_ZB, inputu_ZB, inputSu_ZB = paddingu(uzb[:, :, t], bd_ZB20[1][:,:,t], nfft[1])
        fft2pow2radial!(outu_ZB, rfft(inputu_ZB), rfft(inputSu_ZB), nfft...)
        outv_ZB, inputv_ZB, inputSv_ZB = paddingv(vzb[:, :, t], bd_ZB20[2][:,:,t], nfft[1])
        fft2pow2radial!(outv_ZB, rfft(inputv_ZB), rfft(inputSv_ZB), nfft...)

        bdu_ZB += outu_ZB
        bdv_ZB += outv_ZB

        # 30 day optimization ###############

        # u30, v30, _ = ShallowWaters.add_halo(Float64.(u30s[:,:,t]), Float64.(v30s[:,:,t]), Float64.(eta30s[:,:,t]), zeros(128,128), S30);
        # ShallowWaters.CNN_momentum(u30, v30, S30)
        # outu_30, inputu_30, inputSu_30 = paddingu(u30s[:, :, t], S30.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_30, rfft(inputu_30), rfft(inputSu_30), nfft...)
        # outv_30, inputv_30, inputSv_30 = paddingv(v30s[:, :, t], S30.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_30, rfft(inputv_30), rfft(inputSv_30), nfft...)

        # bdu_30 += outu_30
        # bdv_30 += outv_30

        # batched 2 day ##################################

        outu_multi2, inputu_multi2, inputSu_multi2 = paddingu(umulti2[:, :, t], bd_multi2[1][:,:,t], nfft[1])
        fft2pow2radial!(outu_multi2, rfft(inputu_multi2), rfft(inputSu_multi2), nfft...)
        outv_multi2, inputv_multi2, inputSv_multi2 = paddingv(vmulti2[:, :, t], bd_multi2[2][:,:,t], nfft[1])
        fft2pow2radial!(outv_multi2, rfft(inputv_multi2), rfft(inputSv_multi2), nfft...)

        bdu_multi2 += outu_multi2
        bdv_multi2 += outv_multi2

        # batched 3 day ###########################

        outu_multi3, inputu_multi3, inputSu_multi3 = paddingu(umulti3[:, :, t], bd_multi3[1][:,:,t], nfft[1])
        fft2pow2radial!(outu_multi3, rfft(inputu_multi3), rfft(inputSu_multi3), nfft...)
        outv_multi3, inputv_multi3, inputSv_multi3 = paddingv(vmulti3[:, :, t], bd_multi3[2][:,:,t], nfft[1])
        fft2pow2radial!(outv_multi3, rfft(inputv_multi3), rfft(inputSv_multi3), nfft...)

        bdu_multi3 += outu_multi3
        bdv_multi3 += outv_multi3

        # batched 10 day #####################

        outu_multi10, inputu_multi10, inputSu_multi10 = paddingu(umulti10[:, :, t], bd_multi10[1][:,:,t], nfft[1])
        fft2pow2radial!(outu_multi10, rfft(inputu_multi10), rfft(inputSu_multi10), nfft...)
        outv_multi10, inputv_multi10, inputSv_multi10 = paddingv(vmulti10[:, :, t], bd_multi10[2][:,:,t], nfft[1])
        fft2pow2radial!(outv_multi10, rfft(inputv_multi10), rfft(inputSv_multi10), nfft...)

        bdu_multi10 += outu_multi10
        bdv_multi10 += outv_multi10

    end

end

function ketransfer_check()

    T = Float64
    PZB = ShallowWaters.Parameter(T=T,
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
        zb_forcing_dissipation=true,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=false,
        N=1,
        α=2,
        nx=128,
        Ndays=1
    );

    SZB = ShallowWaters.model_setup(PZB);

    lr_freq = 1/30 .* freq(periodogram(umulti1[:,:,10]; radialavg=true, radialsum=false));
    nfft = nextfastfft(size(uhrcgall[:,:,1]))

    totalu_ZB_check = zeros(65)
    totalv_ZB_check = zeros(65)

    totalstates = 1096
    for t = 1:totalstates

        # ZB20 ############################

        uzb_, vzb_, _ = ShallowWaters.add_halo(Float64.(uzb[:,:,t]), Float64.(vzb[:,:,t]), Float64.(etazb[:,:,t]), zeros(128,128), SZB);
        ShallowWaters.ZB_momentum(uzb_, vzb_, SZB, SZB.Diag);

        outu_ZB, inputu_ZB, inputSu_ZB = paddingu(uzb[:, :, t], SZB.Diag.ZBVars.S_u + Mu_all[:,:,t], nfft[1])
        fft2pow2radial!(outu_ZB, rfft(inputu_ZB), rfft(inputSu_ZB), nfft...)
        outv_ZB, inputv_ZB, inputSv_ZB = paddingv(vzb[:, :, t], SZB.Diag.ZBVars.S_v + Mv_all[:,:,t], nfft[1])
        fft2pow2radial!(outv_ZB, rfft(inputv_ZB), rfft(inputSv_ZB), nfft...)

        totalu_ZB_check += outu_ZB
        totalv_ZB_check += outv_ZB

    end

end

function ketransfer_plots()

    colors = Makie.wong_colors()

    xvals = Int((300/30)):Int((3840-300)/30)
    yvals = Int((300/30)):Int((3840-300)/30)
    # xvals = :
    # yvals = :

    lr_freq = 1/30 .* freq(periodogram(umulti1[xvals,yvals,10]; radialavg=true, radialsum=false));
    nfft = nextfastfft(size(uhrcgall[xvals,yvals,1]))

    # this scaling s is determined as 
    #   (1) a division of the total states we average over (1096), 
    #   (2) a division by \Delta, because I want to compare the S values that are multiplied by dt,
    #       and in Milan's code the forcing that is multiplied by dt has a division by \Delta, which is the
    #       \Delta that I divide by here. The final division is 64 = model.grid.scale, and this is divided because the 
    #       parameterizations are added to a haloed grid, which has a scale factor, and we don't want that scale factor when we compute
    #       the KE transfer
    # Confusingly, the hrcg and approx ones have different scalings. Approximation to the true SGS forcing
    # has no division by s because I never added the scaling factor, and the true SGS forcing doesn't because I divded by the scale
    # when I computed it from my single step function. 
    # Lastly, the true SGS forcing has a multiplication by -1 because when you math
    # out the nonlinear advection term, the nonlinear advection (and thus all the parameterizations), should be approximating -1 
    # time the total tendency term that I computed and saved
    s = 1096 * 30000 * 64
    fig = Figure(size=(800, 350), fontsize=15);
    ax = Axis(fig[1,1],
        xscale = log10,
        xlabel="Wavenumber (1/km)",
        ylabel="KE(k)",
        title="Kinetic energy transfer"
    )
    # lines!(ax, -lr_freq.*(totalu_hrcg + totalv_hrcg)./(1096), label="Total SGS forcing", color=:black)
    # I accidentally divided by \Delta^2 when I computed the SGS forcing from the nonlinear advection approximation, so that's why
    # this one has a multiplication by \Delta
    lines!(ax, (lr_freq * 30000).*(totalu_approx + totalv_approx) ./ (1096), label="Approximate SGS forcing", color=:red)
    lines!(ax, lr_freq.*(totalu_ZB + totalv_ZB) ./ s, label="ZB20", color=colors[1])
    lines!(ax, lr_freq.*(totalu_multi2 + totalv_multi2) ./ s, label="Ensemble 2 day", color=colors[2])

    # lines!(ax, lr_freq.*(totalu_30 + totalv_30) / s, label="30 day",color=colors[6])
    # lines!(ax, lr_freq.*(totalu_20 + totalv_20) / s, label="20 day", color=colors[7])
    lines!(ax, lr_freq.*(totalu_multi3 + totalv_multi3) ./ s, label="Ensemble 3 day", color=colors[3])
    lines!(ax, lr_freq.*(totalu_multi10 + totalv_multi10) ./ s, label="Ensemble 10 day", color=colors[4])
    # lines!(ax, lr_freq.*(totalu_multi20 + totalv_multi20) / s, label="Ensemble 20 day")
    axislegend(ax, position=:rt)
    # Legend(fig[1,2], ax)

    fig = Figure(size=(800,350),fontsize=15);
    s = 1096 * 30000 * 64
    ax = Axis(fig[1,1],
        xscale = log10,
        xlabel="Wavenumber (1/km)",
        ylabel="KE(k)",
        title="Kinetic energy transfer of viscosity"
    )
    lines!(ax, lr_freq.*(viscu_hrcg + viscv_hrcg) ./ s, label="Total SGS forcing", color=:black)
    lines!(ax, lr_freq.*(viscu_ZB + viscv_ZB) ./ s, label="ZB20", color=colors[1])
    lines!(ax, lr_freq.*(viscu_multi2 + viscv_multi2) ./ s, label="Ensemble 2 day",color=colors[2])
    lines!(ax, lr_freq.*(viscu_multi3 + viscv_multi3) ./ s, label="Ensemble 3 day",color=colors[3])
    # lines!(ax, lr_freq.*(totalu_multi5 + totalv_multi5) / 1096, label="Multi 5 day state optimization")
    lines!(ax, lr_freq.*(viscu_multi10 + viscv_multi10) ./ s, label="Ensemble 10 day",color=colors[4])
    # lines!(ax, lr_freq.*(totalu_multi20 + totalv_multi20) / s, label="Ensemble 20 day")
    axislegend(ax, position=:rb)
    # Legend(fig[2,2], ax)

    fig = Figure(size=(800,350),fontsize=15);
    s = 1096 * 30000 * 64
    ax = Axis(fig[1,1],
        xscale = log10,
        xlabel="Wavenumber (1/km)",
        ylabel="KE(k)",
        title="Kinetic energy transfer of bottom drag"
    )
    lines!(ax, lr_freq.*(bdu_hrcg + bdv_hrcg) ./ s, label="Total SGS forcing", color=:black)
    lines!(ax, lr_freq.*(bdu_ZB + bdv_ZB) ./ s, label="ZB20", color=colors[1])
    lines!(ax, lr_freq.*(bdu_multi2 + bdv_multi2) ./ s, label="Ensemble 2 day",color=colors[2])
    lines!(ax, lr_freq.*(bdu_multi3 + bdv_multi3) ./ s, label="Ensemble 3 day",color=colors[3])
    # lines!(ax, lr_freq.*(totalu_multi5 + totalv_multi5) / 1096, label="Multi 5 day state optimization")
    lines!(ax, lr_freq.*(bdu_multi10 + bdv_multi10) ./ s, label="Ensemble 10 day",color=colors[4])
    # lines!(ax, lr_freq.*(totalu_multi20 + totalv_multi20) / s, label="Ensemble 20 day")
    axislegend(ax, position=:rb)
    # Legend(fig[2,2], ax)

    s = 1096 * 30000 * 64
    fig = Figure(size=(800, 300), fontsize=15);
    ax = Axis(fig[1,1],
        xscale = log10,
        xlabel="Wavenumber (1/km)",
        ylabel="KE(k)",
        title="Kinetic energy transfer of S + Viscosity"
    )
    lines!(ax, lr_freq.*((-totalu_hrcg - totalv_hrcg)/1096 + (viscu_hrcg + viscv_hrcg)./s), label="Subgrid forcing",color=:black)
    lines!(ax, lr_freq.*(totalu_ZB + totalv_ZB + viscu_ZB + viscv_ZB) / s, label="ZB20",color=colors[1])
    lines!(ax, lr_freq.*(totalu_multi2 + totalv_multi2 + viscu_multi2 + viscv_multi2) / s, label="Ensemble 2 day",color=colors[2])

    # lines!(ax, lr_freq.*(totalu_30 + totalv_30) / s, label="30 day")
    # lines!(ax, lr_freq.*(totalu_20 + totalv_20) / s, label="20 day")
    lines!(ax, lr_freq.*(totalu_multi3 + totalv_multi3 + viscu_multi3 + viscv_multi3) / s, label="Ensemble 3 day",color=colors[3])
    # lines!(ax, lr_freq.*(totalu_multi5 + totalv_multi5) / 1096, label="Multi 5 day state optimization")
    lines!(ax, lr_freq.*(totalu_multi10 + totalv_multi10 + viscu_multi10 + viscv_multi10) / s, label="Ensemble 10 day",color=colors[4])
    # lines!(ax, lr_freq.*(totalu_multi20 + totalv_multi20) / s, label="Ensemble 20 day")
    axislegend(ax, position=:rt)

    s = 1096 * 30000 * 64
    fig = Figure(size=(800, 300), fontsize=15);
    ax = Axis(fig[1,1],
        xscale = log10,
        xlabel="Wavenumber (1/km)",
        ylabel="KE(k)",
        title="Kinetic energy transfer of S + Viscosity + Bottom drag"
    )
    lines!(ax, lr_freq.*((-totalu_hrcg - totalv_hrcg)/1096 + (viscu_hrcg + viscv_hrcg + bdu_hrcg + bdv_hrcg)./s), label="Subgrid forcing",color=:black)
    lines!(ax, lr_freq.*(totalu_ZB + totalv_ZB + viscu_ZB + viscv_ZB + bdu_ZB + bdv_ZB) / s, label="ZB20",color=colors[1])
    lines!(ax, lr_freq.*(totalu_multi2 + totalv_multi2 + viscu_multi2 + viscv_multi2 +bdu_multi2 + bdv_multi2) / s, label="Ensemble 2 day",color=colors[2])

    # lines!(ax, lr_freq.*(totalu_30 + totalv_30) / s, label="30 day")
    # lines!(ax, lr_freq.*(totalu_20 + totalv_20) / s, label="20 day")
    lines!(ax, lr_freq.*(totalu_multi3 + totalv_multi3 + viscu_multi3 + viscv_multi3 + bdu_multi3 + bdv_multi3) / s, label="Ensemble 3 day",color=colors[3])
    # lines!(ax, lr_freq.*(totalu_multi5 + totalv_multi5) / 1096, label="Multi 5 day state optimization")
    lines!(ax, lr_freq.*(totalu_multi10 + totalv_multi10 + viscu_multi10 + viscv_multi10 + bdu_multi10 + bdv_multi10) / s, label="Ensemble 10 day",color=colors[4])
    # lines!(ax, lr_freq.*(totalu_multi20 + totalv_multi20) / s, label="Ensemble 20 day")
    axislegend(ax, position=:rt)

    s = 1096 * 30000 * 64
    fig = Figure(size=(800, 600), fontsize=15);

    lr_freq = 1/30 .* freq(periodogram(umulti1[:,:,10]; radialavg=true, radialsum=false));
    ax = Axis(fig[1,1],
        xscale = log10,
        xlabel="Wavenumber (1/km)",
        ylabel="KE(k)",
        title="Kinetic energy transfer on whole domain, 3840 km by 3840 km"
    )
    lines!(ax, -lr_freq.*(wholedomain_hrcg)/(1096), label="Total SGS forcing", color=:black)
    # I accidentally divided by \Delta^2 when I computed the SGS forcing from the nonlinear advection approximation, so that's why
    # this one has a multiplication by \Delta
    lines!(ax, (lr_freq * 30000).*(wholedomain_approx) ./ (1096), label="Approximate SGS forcing", color=:red)
    lines!(ax, lr_freq.*(wholedomain_ZB) ./ s, label="ZB20", color=colors[1])
    lines!(ax, lr_freq.*(wholedomain_multi2) ./ s, label="Ensemble 2 day", color=colors[2])
    lines!(ax, lr_freq.*(wholedomain_multi3) ./ s, label="Ensemble 3 day", color=colors[3])
    lines!(ax, lr_freq.*(wholedomain_multi10) ./ s, label="Ensemble 10 day", color=colors[4])

    xvals = Int((150/30)):Int((3840-150)/30)
    lr_freq = 1/30 .* freq(periodogram(umulti2[xvals,xvals,10]; radialavg=true, radialsum=false));
    ax2 = Axis(fig[2,1],
        xscale = log10,
        xlabel="Wavenumber (1/km)",
        ylabel="KE(k)",
        title="Kinetic energy transfer on subdomain of 3690 km by 3690 km"
    )
    lines!(ax2, -lr_freq.*(minus150_hrcg)/(1096), label="Total SGS forcing", color=:black)
    # I accidentally divided by \Delta^2 when I computed the SGS forcing from the nonlinear advection approximation, so that's why
    # this one has a multiplication by \Delta
    lines!(ax2, (lr_freq * 30000).*(minus150_approx) ./ (1096), label="Approximate SGS forcing", color=:red)
    lines!(ax2, lr_freq.*(minus150_ZB) ./ s, label="ZB20", color=colors[1])
    lines!(ax2, lr_freq.*(minus150_multi2) ./ s, label="Ensemble 2 day", color=colors[2])
    lines!(ax2, lr_freq.*(minus150_multi3) ./ s, label="Ensemble 3 day", color=colors[3])
    lines!(ax2, lr_freq.*(minus150_multi10) ./ s, label="Ensemble 10 day", color=colors[4])

    xvals = Int((300/30)):Int((3840-300)/30)
    lr_freq = 1/30 .* freq(periodogram(umulti2[xvals,xvals,10]; radialavg=true, radialsum=false));
    ax3 = Axis(fig[3,1],
        xscale = log10,
        xlabel="Wavenumber (1/km)",
        ylabel="KE(k)",
        title="Kinetic energy transfer on subdomain of 3540 km by 3540 km"
    )
    lines!(ax3, -lr_freq.*(minus300_hrcg)/(1096), label="Total SGS forcing", color=:black)
    # I accidentally divided by \Delta^2 when I computed the SGS forcing from the nonlinear advection approximation, so that's why
    # this one has a multiplication by \Delta
    lines!(ax3, (lr_freq * 30000).*(minus300_approx) ./ (1096), label="Approximate SGS forcing", color=:red)
    lines!(ax3, lr_freq.*(minus300_ZB) ./ s, label="ZB20", color=colors[1])
    lines!(ax3, lr_freq.*(minus300_multi2) ./ s, label="Ensemble 2 day", color=colors[2])
    lines!(ax3, lr_freq.*(minus300_multi3) ./ s, label="Ensemble 3 day", color=colors[3])
    lines!(ax3, lr_freq.*(minus300_multi10) ./ s, label="Ensemble 10 day", color=colors[4])

    Legend(fig[1:3, 2], ax3, orientation = :vertical)

    ga = fig[1, 1] = GridLayout()
    gb = fig[2, 1] = GridLayout()
    gc = fig[3, 1] = GridLayout()

    for (label, layout) in zip(["(a)", "(b)", "(c)"], [ga, gb, gc])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # comparing windowed results
    s = 1096 * 30000 * 64
    fig = Figure(size=(800, 600), fontsize=15);

    lr_freq = 1/30 .* freq(periodogram(umulti1[:,:,10]; radialavg=true, radialsum=false));
    ax = Axis(fig[1,1],
        xscale = log10,
        xlabel="Wavenumber (1/km)",
        ylabel="KE(k)",
        title="KE transfer, no window"
    )
    lines!(ax, -lr_freq.*(totalu_hrcg + totalv_hrcg)./(1096), label="Total SGS forcing", color=:black)
    # I accidentally divided by \Delta^2 when I computed the SGS forcing from the nonlinear advection approximation, so that's why
    # this one has a multiplication by \Delta
    lines!(ax, (lr_freq * 30000).*(totalu_approx + totalv_approx) ./ (1096), label="Approximate SGS forcing", color=:red)
    lines!(ax, lr_freq.*(totalu_ZB + totalv_ZB) ./ s, label="ZB20", color=colors[1])
    lines!(ax, lr_freq.*(totalu_multi2 + totalv_multi2) ./ s, label="Ensemble 2 day", color=colors[2])
    lines!(ax, lr_freq.*(totalu_multi3 + totalv_multi3) ./ s, label="Ensemble 3 day", color=colors[3])
    lines!(ax, lr_freq.*(totalu_multi10 + totalv_multi10) ./ s, label="Ensemble 10 day", color=colors[4])
    # lines!(ax, lr_freq.*(totalu_multi20 + totalv_multi20) / s, label="Ensemble 20 day")

    ax2 = Axis(fig[2,1],
        xscale = log10,
        xlabel="Wavenumber (1/km)",
        ylabel="KE(k)",
        title="KE transfer, window after"
    )
    lines!(ax2, -lr_freq.*(totalu_hrcg_windowafter + totalv_hrcg_windowafter)./(1096), label="Total SGS forcing", color=:black)
    # I accidentally divided by \Delta^2 when I computed the SGS forcing from the nonlinear advection approximation, so that's why
    # this one has a multiplication by \Delta
    lines!(ax2, (lr_freq * 30000).*(totalu_approx_windowafter + totalv_approx_windowafter) ./ (1096), label="Approximate SGS forcing", color=:red)
    lines!(ax2, lr_freq.*(totalu_ZB_windowafter + totalv_ZB_windowafter) ./ s, label="ZB20", color=colors[1])
    lines!(ax2, lr_freq.*(totalu_multi2_windowafter + totalv_multi2_windowafter) ./ s, label="Ensemble 2 day", color=colors[2])
    lines!(ax2, lr_freq.*(totalu_multi3_windowafter + totalv_multi3_windowafter) ./ s, label="Ensemble 3 day", color=colors[3])
    lines!(ax2, lr_freq.*(totalu_multi10_windowafter + totalv_multi10_windowafter) ./ s, label="Ensemble 10 day", color=colors[4])
    # lines!(ax2, lr_freq.*(totalu_multi20_windowafter + totalv_multi20_windowafter) / s, label="Ensemble 20 day")

    ax3 = Axis(fig[3,1],
        xscale = log10,
        xlabel="Wavenumber (1/km)",
        ylabel="KE(k)",
        title="KE transfer, window before"
    )
    # lines!(ax3, -lr_freq.*(totalu_hrcg_windowbefore + totalv_hrcg_windowbefore)./(1096), label="Total SGS forcing", color=:black)
    # I accidentally divided by \Delta^2 when I computed the SGS forcing from the nonlinear advection approximation, so that's why
    # this one has a multiplication by \Delta
    lines!(ax3, (lr_freq * 30000).*(totalu_approx_windowbefore + totalv_approx_windowbefore) ./ (1096), label="Approximate SGS forcing", color=:red)
    lines!(ax3, lr_freq.*(totalu_ZB_windowbefore + totalv_ZB_windowbefore) ./ s, label="ZB20", color=colors[1])
    lines!(ax3, lr_freq.*(totalu_multi2_windowbefore + totalv_multi2_windowbefore) ./ s, label="Ensemble 2 day", color=colors[2])
    lines!(ax3, lr_freq.*(totalu_multi3_windowbefore + totalv_multi3_windowbefore) ./ s, label="Ensemble 3 day", color=colors[3])
    lines!(ax3, lr_freq.*(totalu_multi10_windowbefore + totalv_multi10_windowbefore) ./ s, label="Ensemble 10 day", color=colors[4])
    # lines!(ax3, lr_freq.*(totalu_multi20_windowbefore + totalv_multi20_windowbefore) / s, label="Ensemble 20 day")
    Legend(fig[1:3, 2], ax, orientation = :vertical)


    fig = Figure(size=(700, 500));
    ax = Axis(fig[1,1],
        xscale = log10,
        xlabel="Wavenumber (1/km)",
        ylabel="KE(k)",
        title="KE transfer, computed three ways"
    )
    lines!(ax, (-lr_freq ).*(totalu_hrcg + totalv_hrcg)./1096, label="Total SGS forcing, no window", color=colors[1])
    lines!(ax, (-lr_freq ).*(totalu_hrcg_windowafter + totalv_hrcg_windowafter)./1096, label="Total SGS forcing, window after", color=colors[2])
    lines!(ax, (-lr_freq ).*(totalu_hrcg_windowbefore + totalv_hrcg_windowbefore)./1096, label="Total SGS forcing, window before", color=colors[3])
    lines!(ax, (-lr_freq ).*(totalu_hrcg_windowbefore_noeta + totalv_hrcg_windowbefore_noeta)./1096, label="Total SGS forcing, window before", color=colors[4])

    axislegend(ax, position=:rt)

end
