function load_models()

    uhrcgall = cat(load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_imfilter_3years_postspinup_dailysaves_correctedsetup.jld2")[1],
        load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_imfilter_7years_startfrom3yearpostspinup_weeklysaves.jld2")[1][:,:,2:end]; dims=3
    );
    vhrcgall = cat(load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_imfilter_3years_postspinup_dailysaves_correctedsetup.jld2")[2],
        load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_imfilter_7years_startfrom3yearpostspinup_weeklysaves.jld2")[2][:,:,2:end]; dims=3
    );
    etahrcgall = cat(load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_imfilter_3years_postspinup_dailysaves_correctedsetup.jld2")[3],
        load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_imfilter_7years_startfrom3yearpostspinup_weeklysaves.jld2")[3][:,:,2:end]; dims=3
    );

    uhrall = cat(
        ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day1-766saves/u.nc", "u"),
        ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day766-end/u.nc", "u")[:,:,2:end],
        ncread("./dissipation_constant/spinup_files/1024_7years_startfrom3yearpostspinup_weeklysaves/u.nc", "u")[:,:,2:end]; dims=3
    );

    vhrall = cat(
        ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day1-766saves/v.nc", "v"),
        ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day766-end/v.nc", "v")[:,:,2:end],
        ncread("./dissipation_constant/spinup_files/1024_7years_startfrom3yearpostspinup_weeklysaves/v.nc", "v")[:,:,2:end]; dims=3
    );

    etahrall = cat(
        ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day1-766saves/eta.nc", "eta"),
        ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day766-end/eta.nc", "eta")[:,:,2:end],
        ncread("./dissipation_constant/spinup_files/1024_7years_startfrom3yearpostspinup_weeklysaves/eta.nc", "eta")[:,:,2:end]; dims=3
    );

    uofflinegelu = ncread("./dissipation_constant/results/128_offlineparam_postspinup_cginitcond_3days_gelu_8hoursaves/u.nc", "u");
    vofflinegelu = ncread("./dissipation_constant/results/128_offlineparam_postspinup_cginitcond_3days_gelu_8hoursaves/v.nc", "v");
    etaofflinegelu = ncread("./dissipation_constant/results/128_offlineparam_postspinup_cginitcond_3days_gelu_8hoursaves/eta.nc", "eta");

    uofflinerelu = ncread("./dissipation_constant/results/128_offlineparam_postspinup_cginitcond_3days_relu_8hoursaves/u.nc", "u");
    vofflinerelu = ncread("./dissipation_constant/results/128_offlineparam_postspinup_cginitcond_3days_relu_8hoursaves/v.nc", "v");
    etaofflinerelu = ncread("./dissipation_constant/results/128_offlineparam_postspinup_cginitcond_3days_relu_8hoursaves/eta.nc", "eta");

    # 3 year

    unoparam = ncread("./dissipation_constant/results/128_noparam_postspinup_cginitcond_3years_dailysaves/u.nc", "u");
    vnoparam = ncread("./dissipation_constant/results/128_noparam_postspinup_cginitcond_3years_dailysaves/v.nc", "v");
    etanoparam = ncread("./dissipation_constant/results/128_noparam_postspinup_cginitcond_3years_dailysaves/eta.nc", "eta");

    u1s = ncread("./dissipation_constant/results/result_online_gelu_stateweights_1dayoptimization_startfromoffline_3years_dailysaves/u.nc", "u");
    v1s = ncread("./dissipation_constant/results/result_online_gelu_stateweights_1dayoptimization_startfromoffline_3years_dailysaves/v.nc", "v");
    eta1s = ncread("./dissipation_constant/results/result_online_gelu_stateweights_1dayoptimization_startfromoffline_3years_dailysaves/eta.nc", "eta");

    u5s = ncread("./dissipation_constant/results/result_online_stateweights_5dayoptimization_startfrom1daystate_3years_dailysaves/u.nc", "u");
    v5s = ncread("./dissipation_constant/results/result_online_stateweights_5dayoptimization_startfrom1daystate_3years_dailysaves/v.nc", "v");
    eta5s = ncread("./dissipation_constant/results/result_online_stateweights_5dayoptimization_startfrom1daystate_3years_dailysaves/eta.nc", "eta");

    u10s = ncread("./dissipation_constant/results/result_online_stateweights_10dayoptimization_startfrom5day_3years_dailysaves/u.nc", "u");
    v10s = ncread("./dissipation_constant/results/result_online_stateweights_10dayoptimization_startfrom5day_3years_dailysaves/v.nc", "v");
    eta10s = ncread("./dissipation_constant/results/result_online_stateweights_10dayoptimization_startfrom5day_3years_dailysaves/eta.nc", "eta");

    u20s = ncread("./dissipation_constant/results/result_online_stateweights_20dayoptimization_startfrom10day_fixedcfl_3years_dailysaves/u.nc", "u");
    v20s = ncread("./dissipation_constant/results/result_online_stateweights_20dayoptimization_startfrom10day_fixedcfl_3years_dailysaves/v.nc", "v");
    eta20s = ncread("./dissipation_constant/results/result_online_stateweights_20dayoptimization_startfrom10day_fixedcfl_3years_dailysaves/eta.nc", "eta");

    u20scD = ncread("./dissipation_constant/results/result_online_state_pluscD_weights_20dayoptimization_startfrom20day_3years_dailysaves/u.nc", "u");
    v20scD = ncread("./dissipation_constant/results/result_online_state_pluscD_weights_20dayoptimization_startfrom20day_3years_dailysaves/v.nc", "v");
    eta20scD = ncread("./dissipation_constant/results/result_online_state_pluscD_weights_20dayoptimization_startfrom20day_3years_dailysaves/eta.nc", "eta");

    u30s = ncread("./dissipation_constant/results/result_online_stateweights_30dayoptimization_startfrom20day_fixedcfl_3years_dailysaves/u.nc", "u");
    v30s = ncread("./dissipation_constant/results/result_online_stateweights_30dayoptimization_startfrom20day_fixedcfl_3years_dailysaves/v.nc", "v");
    eta30s = ncread("./dissipation_constant/results/result_online_stateweights_30dayoptimization_startfrom20day_fixedcfl_3years_dailysaves/eta.nc", "eta");

    umulti1 = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_35initdays_startfrommulti3_3years_dailysaves/u.nc", "u");
    vmulti1 = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_35initdays_startfrommulti3_3years_dailysaves/v.nc", "v");
    etamulti1 = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_35initdays_startfrommulti3_3years_dailysaves/eta.nc", "eta");

    umulti1more = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_1:2:89initdays_startfrommulti3_fewerinitconds_3years_dailysaves/u.nc", "u");
    vmulti1more = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_1:2:89initdays_startfrommulti3_fewerinitconds_3years_dailysaves/v.nc", "v");
    etamulti1more = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_1:2:89initdays_startfrommulti3_fewerinitconds_3years_dailysaves/eta.nc", "eta");

    umulti2 = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/u.nc", "u");
    vmulti2 = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/v.nc", "v");
    etamulti2 = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/eta.nc", "eta");

    umulti3 = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_3-15-30-40-50-60-80-85initdays_startfrom20daystate_3years_dailysaves/u.nc", "u");
    vmulti3 = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_3-15-30-40-50-60-80-85initdays_startfrom20daystate_3years_dailysaves/v.nc", "v");
    etamulti3 = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_3-15-30-40-50-60-80-85initdays_startfrom20daystate_3years_dailysaves/eta.nc", "eta");

    umulti3more = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/u.nc", "u");
    vmulti3more = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/v.nc", "v");
    etamulti3more = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/eta.nc", "eta");

    umulti5 = ncread("./dissipation_constant/results/result_online_multistateweights_5dayoptimization_3-30-50-80initdays_fixedcfl_constdiss_3years_dailysaves/u.nc", "u");
    vmulti5 = ncread("./dissipation_constant/results/result_online_multistateweights_5dayoptimization_3-30-50-80initdays_fixedcfl_constdiss_3years_dailysaves/v.nc", "v");
    etamulti5 = ncread("./dissipation_constant/results/result_online_multistateweights_5dayoptimization_3-30-50-80initdays_fixedcfl_constdiss_3years_dailysaves/eta.nc", "eta");

    umulti10 = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_3years_dailysaves/u.nc", "u");
    vmulti10 = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_3years_dailysaves/v.nc", "v");
    etamulti10 = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_3years_dailysaves/eta.nc", "eta");

    umulti20 = ncread("./dissipation_constant/results/result_online_multistateweights_20dayoptimization_5-25-45-65initdays_startfrom20daystate_3years_dailysaves/u.nc", "u");
    vmulti20 = ncread("./dissipation_constant/results/result_online_multistateweights_20dayoptimization_5-25-45-65initdays_startfrom20daystate_3years_dailysaves/v.nc", "v");
    etamulti20 = ncread("./dissipation_constant/results/result_online_multistateweights_20dayoptimization_5-25-45-65initdays_startfrom20daystate_3years_dailysaves/eta.nc", "eta");

    uzb = ncread("./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves/u.nc", "u");
    vzb = ncread("./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves/v.nc", "v");
    etazb = ncread("./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves/eta.nc", "eta");

    # 10 year

    unoparam10 = ncread("./dissipation_constant/results/128_noparam_postspinup_10years_weeklysaves/u.nc", "u");
    vnoparam10 = ncread("./dissipation_constant/results/128_noparam_postspinup_10years_weeklysaves/v.nc", "v");
    etanoparam10 = ncread("./dissipation_constant/results/128_noparam_postspinup_10years_weeklysaves/eta.nc", "eta");

    u5s10 = ncread("./dissipation_constant/results/result_online_stateweights_5dayoptimization_startfrom1day_10years_weeklysaves/u.nc", "u");
    v5s10 = ncread("./dissipation_constant/results/result_online_stateweights_5dayoptimization_startfrom1day_10years_weeklysaves/v.nc", "v");
    eta5s10 = ncread("./dissipation_constant/results/result_online_stateweights_5dayoptimization_startfrom1day_10years_weeklysaves/eta.nc", "eta");

    u10s10 = ncread("./dissipation_constant/results/result_online_stateweights_10dayoptimization_startfrom5day_10years_weeklysaves/u.nc", "u");
    v10s10 = ncread("./dissipation_constant/results/result_online_stateweights_10dayoptimization_startfrom5day_10years_weeklysaves/v.nc", "v");
    eta10s10 = ncread("./dissipation_constant/results/result_online_stateweights_10dayoptimization_startfrom5day_10years_weeklysaves/eta.nc", "eta");

    u20s10 = ncread("./dissipation_constant/results/result_online_stateweights_20dayoptimization_startfrom10day_fixedcfl_10years_weeklysaves/u.nc", "u");
    v20s10 = ncread("./dissipation_constant/results/result_online_stateweights_20dayoptimization_startfrom10day_fixedcfl_10years_weeklysaves/v.nc", "v");
    eta20s10 = ncread("./dissipation_constant/results/result_online_stateweights_20dayoptimization_startfrom10day_fixedcfl_10years_weeklysaves/eta.nc", "eta");

    u30s10 = ncread("./dissipation_constant/results/result_online_stateweights_30dayoptimization_startfrom20day_fixedcfl_10years_weeklysaves/u.nc", "u");
    v30s10 = ncread("./dissipation_constant/results/result_online_stateweights_30dayoptimization_startfrom20day_fixedcfl_10years_weeklysaves/v.nc", "v");
    eta30s10 = ncread("./dissipation_constant/results/result_online_stateweights_30dayoptimization_startfrom20day_fixedcfl_10years_weeklysaves/eta.nc", "eta");

    umulti110 = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_35initdays_startfrommulti3_10years_weeklysaves/u.nc", "u");
    vmulti110 = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_35initdays_startfrommulti3_10years_weeklysaves/v.nc", "v");
    etamulti110 = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_35initdays_startfrommulti3_10years_weeklysaves/eta.nc", "eta");

    umulti1more10 = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_1:2:89initdays_startfrommulti3_fewerinitconds_10years_weeklysaves/u.nc", "u");
    vmulti1more10 = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_1:2:89initdays_startfrommulti3_fewerinitconds_10years_weeklysaves/v.nc", "v");
    etamulti1more10 = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_1:2:89initdays_startfrommulti3_fewerinitconds_10years_weeklysaves/eta.nc", "eta");

    umulti210 = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_10years_weeklysaves/u.nc", "u");
    vmulti210 = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_10years_weeklysaves/v.nc", "v");
    etamulti210 = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_10years_weeklysaves/eta.nc", "eta");

    umulti310 = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_3-15-30-40-50-60-80-85initdays_startfrom20daystate_10years_weeklysaves/u.nc", "u");
    vmulti310 = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_3-15-30-40-50-60-80-85initdays_startfrom20daystate_10years_weeklysaves/v.nc", "v");
    etamulti310 = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_3-15-30-40-50-60-80-85initdays_startfrom20daystate_10years_weeklysaves/eta.nc", "eta");

    umulti3more10 = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_10years_weeklysaves/u.nc", "u");
    vmulti3more10 = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_10years_weeklysaves/v.nc", "v");
    etamulti3more10 = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_10years_weeklysaves/eta.nc", "eta");

    umulti1010 = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_10years_weeklysaves/u.nc", "u");
    vmulti1010 = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_10years_weeklysaves/v.nc", "v");
    etamulti1010 = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_10years_weeklysaves/eta.nc", "eta");

    umulti2010 = ncread("./dissipation_constant/results/result_online_multistateweights_20dayoptimization_5-25-45-65initdays_startfrom20daystate_10years_weeklysaves/u.nc", "u");
    vmulti2010 = ncread("./dissipation_constant/results/result_online_multistateweights_20dayoptimization_5-25-45-65initdays_startfrom20daystate_10years_weeklysaves/v.nc", "v");
    etamulti2010 = ncread("./dissipation_constant/results/result_online_multistateweights_20dayoptimization_5-25-45-65initdays_startfrom20daystate_10years_weeklysaves/eta.nc", "eta");

    uzb10 = ncread("./dissipation_constant/spinup_files/ZB20_10yearspostspinup_weeklysaves/u.nc", "u");
    vzb10 = ncread("./dissipation_constant/spinup_files/ZB20_10yearspostspinup_weeklysaves/v.nc", "v");
    etazb10 = ncread("./dissipation_constant/spinup_files/ZB20_10yearspostspinup_weeklysaves/eta.nc", "eta");

    # another ten year (for just multi2 to check further stability)

    umulti210more = ncread("./dissipation_constant/results/result_online_multistate_2dayoptimization_further10years_weeklysaves/u.nc", "u");
    vmulti210more = ncread("./dissipation_constant/results/result_online_multistate_2dayoptimization_further10years_weeklysaves/v.nc", "v");
    etamulti210more = ncread("./dissipation_constant/results/result_online_multistate_2dayoptimization_further10years_weeklysaves/eta.nc", "eta");

    umulti310more = ncread("./dissipation_constant/results/result_online_multistate_3dayoptimization_further10years_weeklysaves/u.nc", "u");
    vmulti310more = ncread("./dissipation_constant/results/result_online_multistate_3dayoptimization_further10years_weeklysaves/v.nc", "v");
    etamulti310more = ncread("./dissipation_constant/results/result_online_multistate_3dayoptimization_further10years_weeklysaves/eta.nc", "eta");

    # the following didn't work as loss functions

    ukespec = ncread("./dissipation_constant/results/maybeneed/128_online_gelu_kespecweights_3dayoptimization_startfrom5daystate_3years_dailysaves/u.nc", "u");
    vkespec = ncread("./dissipation_constant/results/maybeneed/128_online_gelu_kespecweights_3dayoptimization_startfrom5daystate_3years_dailysaves/v.nc", "v");
    etakespec = ncread("./dissipation_constant/results/maybeneed/128_online_gelu_kespecweights_3dayoptimization_startfrom5daystate_3years_dailysaves/eta.nc", "eta");

    ukespecpd = ncread("./dissipation_constant/results/maybeneed/128_online_gelu_kespecpdweights_3dayoptimization_startfrom5daystate_3years_dailysaves/u.nc", "u");
    vkespecpd = ncread("./dissipation_constant/results/maybeneed/128_online_gelu_kespecpdweights_3dayoptimization_startfrom5daystate_3years_dailysaves/v.nc", "v");
    etakespecpd = ncread("./dissipation_constant/results/maybeneed/128_online_gelu_kespecpdweights_3dayoptimization_startfrom5daystate_3years_dailysaves/eta.nc", "eta");

    ufourier = ncread("./dissipation_constant/results/maybeneed/128_online_gelu_fourierweights_3dayoptimization_startfrom5daystate_3years_dailysaves/u.nc", "u");
    vfourier = ncread("./dissipation_constant/results/maybeneed/128_online_gelu_fourierweights_3dayoptimization_startfrom5daystate_3years_dailysaves/v.nc", "v");
    etafourier = ncread("./dissipation_constant/results/maybeneed/128_online_gelu_fourierweights_3dayoptimization_startfrom5daystate_3years_dailysaves/eta.nc", "eta");

    uhybrid = ncread("./dissipation_constant/results/maybeneed/128_online_gelu_hybridweights_3dayoptimization_startfrom5daystate_3years_dailysaves/u.nc", "u");
    vhybrid = ncread("./dissipation_constant/results/maybeneed/128_online_gelu_hybridweights_3dayoptimization_startfrom5daystate_3years_dailysaves/v.nc", "v");
    etahybrid = ncread("./dissipation_constant/results/maybeneed/128_online_gelu_hybridweights_3dayoptimization_startfrom5daystate_3years_dailysaves/eta.nc", "eta");

    # relu activation function

    u1daystaterelu = ncread("./dissipation_constant/results/result_online_reluactivation_stateweights_1dayoptimization_startfromoffline_madnlp_30days_8hoursaves/u.nc", "u");
    v1daystaterelu = ncread("./dissipation_constant/results/result_online_reluactivation_stateweights_1dayoptimization_startfromoffline_madnlp_30days_8hoursaves/v.nc", "v");
    eta1daystaterelu = ncread("./dissipation_constant/results/result_online_reluactivation_stateweights_1dayoptimization_startfromoffline_madnlp_30days_8hoursaves/eta.nc", "eta");

    u5daystaterelu = ncread("./dissipation_constant/results/128_online_reluactivation_stateweights_5dayoptimization_startfrom1daystate_madnlp_30days_8hoursaves/u.nc", "u");
    v5daystaterelu = ncread("./dissipation_constant/results/128_online_reluactivation_stateweights_5dayoptimization_startfrom1daystate_madnlp_30days_8hoursaves/v.nc", "v");
    eta5daystaterelu = ncread("./dissipation_constant/results/128_online_reluactivation_stateweights_5dayoptimization_startfrom1daystate_madnlp_30days_8hoursaves/eta.nc", "eta");

    ukespecpd1dayrelu = ncread("./dissipation_constant/results/128_online_reluactivation_kespecpdweights_1dayoptimization_startfrom1daystate_madnlp_30days_8hoursaves/u.nc", "u");
    vkespecpd1dayrelu = ncread("./dissipation_constant/results/128_online_reluactivation_kespecpdweights_1dayoptimization_startfrom1daystate_madnlp_30days_8hoursaves/v.nc", "v");
    etakespecpd1dayrelu = ncread("./dissipation_constant/results/128_online_reluactivation_kespecpdweights_1dayoptimization_startfrom1daystate_madnlp_30days_8hoursaves/eta.nc", "eta");

    ker = ImageFiltering.Kernel.gaussian((30e3/3750));
    # imfilter(hru[:,:,j], reflect(ker))

    # adding the above loads into a single file, making it easier to plot
    # entries 1 - 1096 correspond to daily saves over the first three years, 1097-1451 correspond to the next 7 years at weekly saves
    # these do not contain the first three years of weekly saves in the 10 year runs, since those are already accounted for
    # multi2 and multi3 have an additional 10 years of weekly saves 

    unoparamall = cat(unoparam, unoparam10[:,:,157:end]; dims=3);
    vnoparamall = cat(vnoparam, vnoparam10[:,:,157:end]; dims=3);
    etanoparamall = cat(etanoparam, etanoparam10[:,:,157:end]; dims=3);

    uzball = cat(uzb, uzb10[:,:,157:end]; dims=3);
    vzball = cat(vzb, vzb10[:,:,157:end]; dims=3);
    etazball = cat(etazb, etazb10[:,:,157:end]; dims=3);

    u10sall = cat(u10s, u10s10[:,:,157:end]; dims=3);
    v10sall = cat(v10s, v10s10[:,:,157:end]; dims=3);
    eta10sall = cat(eta10s, eta10s10[:,:,157:end]; dims=3);

    u20sall = cat(u20s, u20s10[:,:,157:end]; dims=3);
    v20sall = cat(v20s, v20s10[:,:,157:end]; dims=3);
    eta20sall = cat(eta20s, eta20s10[:,:,157:end]; dims=3);

    u30sall = cat(u30s, u30s10[:,:,157:end]; dims=3);
    v30sall = cat(v30s, v30s10[:,:,157:end]; dims=3);
    eta30sall = cat(eta30s, eta30s10[:,:,157:end]; dims=3);

    umulti2all = cat(umulti2, umulti210[:,:,157:end], umulti210more[:,:,2:end]; dims=3);
    vmulti2all = cat(vmulti2, vmulti210[:,:,157:end], vmulti210more[:,:,2:end]; dims=3);
    etamulti2all = cat(etamulti2, etamulti210[:,:,157:end], etamulti210more[:,:,2:end]; dims=3);

    umulti3all = cat(umulti3more, umulti3more10[:,:,157:end], umulti310more[:,:,2:end]; dims=3);
    vmulti3all = cat(vmulti3more, vmulti3more10[:,:,157:end], vmulti310more[:,:,2:end]; dims=3);
    etamulti3all = cat(etamulti3more, etamulti3more10[:,:,157:end], etamulti310more[:,:,2:end]; dims=3);

    umulti10all = cat(umulti10, umulti1010[:,:,157:end]; dims=3);
    vmulti10all = cat(vmulti10, vmulti1010[:,:,157:end]; dims=3);
    etamulti10all = cat(etamulti10, etamulti1010[:,:,157:end]; dims=3);

    umulti20all = cat(umulti20, umulti2010[:,:,157:end]; dims=3);
    vmulti20all = cat(vmulti20, vmulti2010[:,:,157:end]; dims=3);
    etamulti20all = cat(etamulti20, etamulti2010[:,:,157:end]; dims=3);

end

function load_Sfiles()

    true_S = load_object("./dissipation_constant/computing_trueS/trueS_fromtendencies_withrk4_SuSv_first3years_dailysaves_032526.jld2");
    Suhr = true_S[1];
    Svhr = true_S[2];

    approx_S = load_object("./dissipation_constant/computing_trueS/approxS_nonlinearadvec_nottendencies_hasextraDelta_SuSv_first3years_040726.jld2");
    Suapprox = approx_S[1];
    Svapprox = approx_S[2];

    # visc_hr = load_object("./dissipation_constant/alternate_S_files/viscosity_hr_3years_dailysaves_MuMv.jld2");
    visc_cg = load_object("./dissipation_constant/alternate_S_files/viscosity_cghr_3years_dailysaves_MuMv.jld2");

    # bd_hr = load_object("./dissipation_constant/alternate_S_files/bottomdrag_hr_3years_dailysaves_BuBv.jld2");
    bd_cg = load_object("./dissipation_constant/alternate_S_files/bottomdrag_cghr_3years_dailysaves_BuBv.jld2");

    # advec_hr = load_object("./nonlinear_advec_fromhrstates_advu_advv.jld2");
    advec_cg = load_object("./dissipation_constant/alternate_S_files/nonlinear_advec_fromcghrstates_advu_advv.jld2");

    tend_cg = load_object("./dissipation_constant/computing_trueS/coarsegrain_tendencies_doverline_rk1_dudvdeta.jld2");
    # tend_hr = load_object("./dissipation_constant/computing_trueS/highresolution_tendencies_rk1_dudvdeta.jld2");

    # coarse-grain the high-resolution viscosity terms
    # ker = ImageFiltering.Kernel.gaussian((30e3/3750));
    # tendufiltered = zeros(1023, 1024, 1096);
    # tendvfiltered = zeros(1024, 1023, 1096);

    # for t = 1:1096
    #     @views tendufiltered[:,:,t] .= imfilter(tend_hr[1][:,:,t], reflect(ker))
    #     @views tendvfiltered[:,:,t] .= imfilter(tend_hr[2][:,:,t], reflect(ker))
    # end

    tendu_hrdownsized = load_object("./dissipation_constant/computing_trueS/tendency_hrcg_overline(d)_rk1_dudv.jld2")[1];
    tendv_hrdownsized = load_object("./dissipation_constant/computing_trueS/tendency_hrcg_overline(d)_rk1_dudv.jld2")[2];

    Mu_hrdownsized = load_object("./dissipation_constant/alternate_S_files/viscosity_coarsegrainedhr_3years_dailysaves_overlineMuoverlineMv.jld2")[1];
    Mv_hrdownsized = load_object("./dissipation_constant/alternate_S_files/viscosity_coarsegrainedhr_3years_dailysaves_overlineMuoverlineMv.jld2")[2];

    Bu_hrdownsized = load_object("./dissipation_constant/alternate_S_files/bottomdrag_coarsegrained_hr_3years_dailysaves_overlineBuoverlineBv.jld2")[1];
    Bv_hrdownsized = load_object("./dissipation_constant/alternate_S_files/bottomdrag_coarsegrained_hr_3years_dailysaves_overlineBuoverlineBv.jld2")[2];

    Advecu_hrdownsized = load_object("./dissipation_constant/alternate_S_files/nonlinearadvec_coarsegrainedhr_3years_dailysaves_overlineAdvecuoverlineAdvecv.jld2")[1];
    Advecv_hrdownsized = load_object("./dissipation_constant/alternate_S_files/nonlinearadvec_coarsegrainedhr_3years_dailysaves_overlineAdvecuoverlineAdvecv.jld2")[2];

    P = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        cfl=.898,
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
        # nn_forcing_momentum=false,
        # nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=2,
        initial_cond="rest"
    );
    S = ShallowWaters.model_setup(P);

end