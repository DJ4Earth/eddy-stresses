function computing_ketransfer_windowafter()

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

    totalu_hrcg_windowafter = zeros(length(lr_freq))
    totalv_hrcg_windowafter = zeros(length(lr_freq))

    totalu_approx_windowafter = zeros(length(lr_freq))
    totalv_approx_windowafter = zeros(length(lr_freq))

    totalu_multi2_windowafter = zeros(length(lr_freq))
    totalv_multi2_windowafter = zeros(length(lr_freq))

    totalu_multi3_windowafter = zeros(length(lr_freq))
    totalv_multi3_windowafter = zeros(length(lr_freq))

    totalu_multi10_windowafter = zeros(length(lr_freq))
    totalv_multi10_windowafter = zeros(length(lr_freq))

    totalu_multi20_windowafter = zeros(length(lr_freq))
    totalv_multi20_windowafter = zeros(length(lr_freq))

    totalu_ZB_windowafter = zeros(length(lr_freq))
    totalv_ZB_windowafter = zeros(length(lr_freq))

    alpha = 0.1
    winu = tukey((127, 128), alpha)
    winv = tukey((128,127), alpha)
    totalstates = 1096
    for t = 1:totalstates

        # from approximation to nonlinear advection
        outu_approx, inputu_approx, inputSu_approx = paddingu(winu.*uhrcgall[xvals, yvals, t], winu.*Suapprox[xvals,yvals,t], nfft[1])
        fft2pow2radial!(outu_approx, rfft(inputu_approx), rfft(inputSu_approx), nfft...)
        outv_approx, inputv_approx, inputSv_approx = paddingv(winv.*vhrcgall[xvals, yvals, t], winv.*Svapprox[xvals,yvals,t], nfft[1])
        fft2pow2radial!(outv_approx, rfft(inputv_approx), rfft(inputSv_approx), nfft...)

        totalu_approx_windowafter += outu_approx
        totalv_approx_windowafter += outv_approx

        # from total tendencies
        outu_hrcg, inputu_hrcg, inputSu_hrcg = paddingu(winu.*uhrcgall[xvals, yvals, t], winu.*Suhr[xvals,yvals,t], nfft[1])
        fft2pow2radial!(outu_hrcg, rfft(inputu_hrcg), rfft(inputSu_hrcg), nfft...)
        outv_hrcg, inputv_hrcg, inputSv_hrcg = paddingv(winv.*vhrcgall[xvals, yvals, t], winv.*Svhr[xvals,yvals,t], nfft[1])
        fft2pow2radial!(outv_hrcg, rfft(inputv_hrcg), rfft(inputSv_hrcg), nfft...)

        totalu_hrcg_windowafter += outu_hrcg
        totalv_hrcg_windowafter += outv_hrcg

        # ZB20 ############################

        uzb_, vzb_, _ = ShallowWaters.add_halo(Float64.(uzb[:,:,t]), 
            Float64.(vzb[:,:,t]),
            Float64.(etazb[:,:,t]),
            zeros(128,128),
            SZB);
        ShallowWaters.ZB_momentum(uzb_, vzb_, SZB, SZB.Diag);

        outu_ZB, inputu_ZB, inputSu_ZB = paddingu(winu.*uzb[xvals, yvals, t], winu.*SZB.Diag.ZBVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_ZB, rfft(inputu_ZB), rfft(inputSu_ZB), nfft...)
        outv_ZB, inputv_ZB, inputSv_ZB = paddingv(winv.*vzb[xvals, yvals, t], winv.*SZB.Diag.ZBVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_ZB, rfft(inputv_ZB), rfft(inputSv_ZB), nfft...)

        totalu_ZB_windowafter += outu_ZB
        totalv_ZB_windowafter += outv_ZB

        # batched 2 day ##################################

        umulti2_, vmulti2_, _ = ShallowWaters.add_halo(Float64.(umulti2[:,:,t]), Float64.(vmulti2[:,:,t]), Float64.(etamulti2[:,:,t]), zeros(128,128), Smulti2);
        ShallowWaters.CNN_momentum(umulti2_, vmulti2_, Smulti2)

        outu_multi2, inputu_multi2, inputSu_multi2 = paddingu(winu.*umulti2[xvals, yvals, t], winu.*Smulti2.Diag.CNNVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_multi2, rfft(inputu_multi2), rfft(inputSu_multi2), nfft...)
        outv_multi2, inputv_multi2, inputSv_multi2 = paddingv(winv.*vmulti2[xvals, yvals, t], winv.*Smulti2.Diag.CNNVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_multi2, rfft(inputv_multi2), rfft(inputSv_multi2), nfft...)

        totalu_multi2_windowafter += outu_multi2
        totalv_multi2_windowafter += outv_multi2

        # batched 3 day ###########################

        umulti3_, vmulti3_, _ = ShallowWaters.add_halo(Float64.(umulti3[:,:,t]), Float64.(vmulti3[:,:,t]), Float64.(etamulti3[:,:,t]), zeros(128,128), Smulti3);
        ShallowWaters.CNN_momentum(umulti3_, vmulti3_, Smulti3)

        outu_multi3, inputu_multi3, inputSu_multi3 = paddingu(winu.*umulti3[xvals, yvals, t], winu.*Smulti3.Diag.CNNVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_multi3, rfft(inputu_multi3), rfft(inputSu_multi3), nfft...)
        outv_multi3, inputv_multi3, inputSv_multi3 = paddingv(winv.*vmulti3[xvals, yvals, t], winv.*Smulti3.Diag.CNNVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_multi3, rfft(inputv_multi3), rfft(inputSv_multi3), nfft...)

        totalu_multi3_windowafter += outu_multi3
        totalv_multi3_windowafter += outv_multi3

        # batched 10 day #####################

        umulti10_, vmulti10_, _ = ShallowWaters.add_halo(Float64.(umulti10[:,:,t]), Float64.(vmulti10[:,:,t]), Float64.(etamulti10[:,:,t]), zeros(128,128), Smulti10);
        ShallowWaters.CNN_momentum(umulti10_, vmulti10_, Smulti10)
        outu_multi10, inputu_multi10, inputSu_multi10 = paddingu(winu.*umulti10[xvals, yvals, t], winu.*Smulti10.Diag.CNNVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_multi10, rfft(inputu_multi10), rfft(inputSu_multi10), nfft...)
        outv_multi10, inputv_multi10, inputSv_multi10 = paddingv(winv.*vmulti10[xvals, yvals, t], winv.*Smulti10.Diag.CNNVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_multi10, rfft(inputv_multi10), rfft(inputSv_multi10), nfft...)

        totalu_multi10_windowafter += outu_multi10
        totalv_multi10_windowafter += outv_multi10

        # batched 20 day #####################

        umulti20_, vmulti20_, _ = ShallowWaters.add_halo(Float64.(umulti20[:,:,t]), Float64.(vmulti20[:,:,t]), Float64.(etamulti20[:,:,t]), zeros(128,128), Smulti20);
        ShallowWaters.CNN_momentum(umulti20_, vmulti20_, Smulti20)
        outu_multi20, inputu_multi20, inputSu_multi20 = paddingu(winu.*umulti20[xvals, yvals, t], winu.*Smulti20.Diag.CNNVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_multi20, rfft(inputu_multi20), rfft(inputSu_multi20), nfft...)
        outv_multi20, inputv_multi20, inputSv_multi20 = paddingv(winv.*vmulti20[xvals, yvals, t], winv.*Smulti20.Diag.CNNVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_multi20, rfft(inputv_multi20), rfft(inputSv_multi20), nfft...)

        totalu_multi20_windowafter += outu_multi20
        totalv_multi20_windowafter += outv_multi20

    end

end

function computing_ketransfer_windowbefore()

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

    # true_S = load_object("./dissipation_constant/ke_transfers/trueS_fromwindowedtendencies_withRK4_SuSv_first3years_dailysaves.jld2");
    # Suhr = true_S[1];
    # Svhr = true_S[2];
    true_S = load_object("./dissipation_constant/ke_transfers/trueS_fromwindowedtendencies_withRK4_SuSv_first3years_dailysaves_S(overline(window(u))).jld2");
    Suhr = true_S[1];
    Svhr = true_S[2];

    approx_S = load_object("./dissipation_constant/ke_transfers/approxS_nonlinearadvec_nottendencies_fromwindowedstates_hasextraDelta_SuSv_first3years.jld2");
    Suapprox = approx_S[1];
    Svapprox = approx_S[2];

    # for computing a KE transfer from some/all of the domain
    xvals = :
    # xvals = Int((150/30)):Int((3840-150)/30)
    # xvals = Int((300/30)):Int((3840-300)/30)
    yvals = xvals

    lr_freq = 1/30 .* freq(periodogram(umulti2[xvals,yvals,10]; radialavg=true, radialsum=false));
    nfft = nextfastfft(size(uhrcgall[xvals,yvals,1]))

    totalu_hrcg_windowbefore = zeros(length(lr_freq))
    totalv_hrcg_windowbefore = zeros(length(lr_freq))

    totalu_approx_windowbefore = zeros(length(lr_freq))
    totalv_approx_windowbefore = zeros(length(lr_freq))

    totalu_multi2_windowbefore = zeros(length(lr_freq))
    totalv_multi2_windowbefore = zeros(length(lr_freq))

    totalu_multi3_windowbefore = zeros(length(lr_freq))
    totalv_multi3_windowbefore = zeros(length(lr_freq))

    totalu_multi10_windowbefore = zeros(length(lr_freq))
    totalv_multi10_windowbefore = zeros(length(lr_freq))

    totalu_multi20_windowbefore = zeros(length(lr_freq))
    totalv_multi20_windowbefore = zeros(length(lr_freq))

    totalu_ZB_windowbefore = zeros(length(lr_freq))
    totalv_ZB_windowbefore = zeros(length(lr_freq))

    alpha = 0.1
    winu = tukey((SZB.grid.nux, SZB.grid.nuy), alpha)
    winv = tukey((SZB.grid.nvx, SZB.grid.nvy), alpha)
    totalstates = 1096
    for t = 1:totalstates

        # from approximation to nonlinear advection
        outu_approx, inputu_approx, inputSu_approx = paddingu(winu.*uhrcgall[xvals, yvals, t], Suapprox[xvals,yvals,t], nfft[1])
        fft2pow2radial!(outu_approx, rfft(inputu_approx), rfft(inputSu_approx), nfft...)
        outv_approx, inputv_approx, inputSv_approx = paddingv(winv.*vhrcgall[xvals, yvals, t], Svapprox[xvals,yvals,t], nfft[1])
        fft2pow2radial!(outv_approx, rfft(inputv_approx), rfft(inputSv_approx), nfft...)

        totalu_approx_windowbefore += outu_approx
        totalv_approx_windowbefore += outv_approx

        # from total tendencies
        outu_hrcg, inputu_hrcg, inputSu_hrcg = paddingu(winu.*uhrcgall[xvals, yvals, t], Subefore[xvals,yvals,t], nfft[1])
        fft2pow2radial!(outu_hrcg, rfft(inputu_hrcg), rfft(inputSu_hrcg), nfft...)
        outv_hrcg, inputv_hrcg, inputSv_hrcg = paddingv(winv.*vhrcgall[xvals, yvals, t], Svbefore[xvals,yvals,t], nfft[1])
        fft2pow2radial!(outv_hrcg, rfft(inputv_hrcg), rfft(inputSv_hrcg), nfft...)

        totalu_hrcg_windowbefore += outu_hrcg
        totalv_hrcg_windowbefore += outv_hrcg

        # ZB20 ############################

        uzb_, vzb_, _ = ShallowWaters.add_halo(Float64.(winu.*uzb[:,:,t]), 
            Float64.(winv.*vzb[:,:,t]),
            Float64.(etazb[:,:,t]),
            zeros(128,128),
            SZB);
        ShallowWaters.ZB_momentum(uzb_, vzb_, SZB, SZB.Diag);

        outu_ZB, inputu_ZB, inputSu_ZB = paddingu(winu.*uzb[xvals, yvals, t], SZB.Diag.ZBVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_ZB, rfft(inputu_ZB), rfft(inputSu_ZB), nfft...)
        outv_ZB, inputv_ZB, inputSv_ZB = paddingv(winv.*vzb[xvals, yvals, t], SZB.Diag.ZBVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_ZB, rfft(inputv_ZB), rfft(inputSv_ZB), nfft...)

        totalu_ZB_windowbefore += outu_ZB
        totalv_ZB_windowbefore += outv_ZB


        # batched 2 day ##################################

        umulti2_, vmulti2_, _ = ShallowWaters.add_halo(Float64.(winu.*umulti2[:,:,t]), Float64.(winv.*vmulti2[:,:,t]), Float64.(etamulti2[:,:,t]), zeros(128,128), Smulti2);
        ShallowWaters.CNN_momentum(umulti2_, vmulti2_, Smulti2)

        outu_multi2, inputu_multi2, inputSu_multi2 = paddingu(winu.*umulti2[xvals, yvals, t], Smulti2.Diag.CNNVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_multi2, rfft(inputu_multi2), rfft(inputSu_multi2), nfft...)
        outv_multi2, inputv_multi2, inputSv_multi2 = paddingv(winv.*vmulti2[xvals, yvals, t], Smulti2.Diag.CNNVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_multi2, rfft(inputv_multi2), rfft(inputSv_multi2), nfft...)

        totalu_multi2_windowbefore += outu_multi2
        totalv_multi2_windowbefore += outv_multi2

        # batched 3 day ###########################

        umulti3_, vmulti3_, _ = ShallowWaters.add_halo(Float64.(winu.*umulti3[:,:,t]), Float64.(winv.*vmulti3[:,:,t]), Float64.(etamulti3[:,:,t]), zeros(128,128), Smulti3);
        ShallowWaters.CNN_momentum(umulti3_, vmulti3_, Smulti3)

        outu_multi3, inputu_multi3, inputSu_multi3 = paddingu(winu.*umulti3[xvals, yvals, t], Smulti3.Diag.CNNVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_multi3, rfft(inputu_multi3), rfft(inputSu_multi3), nfft...)
        outv_multi3, inputv_multi3, inputSv_multi3 = paddingv(winv.*vmulti3[xvals, yvals, t], Smulti3.Diag.CNNVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_multi3, rfft(inputv_multi3), rfft(inputSv_multi3), nfft...)

        totalu_multi3_windowbefore += outu_multi3
        totalv_multi3_windowbefore += outv_multi3

        # batched 10 day #####################

        umulti10_, vmulti10_, _ = ShallowWaters.add_halo(Float64.(winu.*umulti10[:,:,t]), Float64.(winv.*vmulti10[:,:,t]), Float64.(etamulti10[:,:,t]), zeros(128,128), Smulti10);
        ShallowWaters.CNN_momentum(umulti10_, vmulti10_, Smulti10)

        outu_multi10, inputu_multi10, inputSu_multi10 = paddingu(winu.*umulti10[xvals, yvals, t], Smulti10.Diag.CNNVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_multi10, rfft(inputu_multi10), rfft(inputSu_multi10), nfft...)
        outv_multi10, inputv_multi10, inputSv_multi10 = paddingv(winv.*vmulti10[xvals, yvals, t], Smulti10.Diag.CNNVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_multi10, rfft(inputv_multi10), rfft(inputSv_multi10), nfft...)

        totalu_multi10_windowbefore += outu_multi10
        totalv_multi10_windowbefore += outv_multi10

        # batched 20 day #####################

        umulti20_, vmulti20_, _ = ShallowWaters.add_halo(Float64.(winu.*umulti20[:,:,t]), Float64.(winv.*vmulti20[:,:,t]), Float64.(etamulti20[:,:,t]), zeros(128,128), Smulti20);
        ShallowWaters.CNN_momentum(umulti20_, vmulti20_, Smulti20)
        outu_multi20, inputu_multi20, inputSu_multi20 = paddingu(winu.*umulti20[xvals, yvals, t], Smulti20.Diag.CNNVars.S_u[xvals,yvals], nfft[1])
        fft2pow2radial!(outu_multi20, rfft(inputu_multi20), rfft(inputSu_multi20), nfft...)
        outv_multi20, inputv_multi20, inputSv_multi20 = paddingv(winv.*vmulti20[xvals, yvals, t], Smulti20.Diag.CNNVars.S_v[xvals,yvals], nfft[1])
        fft2pow2radial!(outv_multi20, rfft(inputv_multi20), rfft(inputSv_multi20), nfft...)

        totalu_multi20_windowbefore += outu_multi20
        totalv_multi20_windowbefore += outv_multi20

    end

end