"""
Mostly figure generation, I just wanted to be able to run include("technical_paper.jl")
without all of this also running.
"""

function create_models()

    T = Float64
    Ndays = 1
    coarse_grained_hrstates = load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_3years_postspinup_dailysaves.jld2");
    uhrcg = coarse_grained_hrstates[1]
    vhrcg = coarse_grained_hrstates[2]
    etahrcg = coarse_grained_hrstates[3]

    # coarse_grained_hrstates = load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_90days_postspinup_8hoursaves.jld2");
    # uhrcg2 = coarse_grained_hrstates[1];

    Pnoparam = ShallowWaters.Parameter(T=T,
        output=true,
        # output_dt = 1,
        # output_dt=168,
        # output_dt=12600,
        L_ratio=1,
        g=9.81,
        H=500,
        # cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        # diffusion="Smagorinsky",        # this is the only new parameter to be adjusted in the new spinups
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
        Ndays=Ndays,
        # initial_cond="rest",
        initial_cond="ncfile",
        initpath="./dissipation_constant/spinup_files/10yearspinup_128_noslipbc_noforcing_float64prog"
        # initpath = "./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves/1024_postspinup_day1-766saves"
    );

    Snoparam = ShallowWaters.model_setup(Pnoparam);

    u0, v0, eta0, _ = ShallowWaters.add_halo(uhrcg[:,:,1],vhrcg[:,:,1],etahrcg[:,:,1],zeros(128,128),Snoparam);
    initial_cond = [u0, v0, eta0];

    Snoparam.Prog.u .= copy(initial_cond[1]);
    Snoparam.Prog.v .= copy(initial_cond[2]);
    Snoparam.Prog.η .= copy(initial_cond[3]);

    ShallowWaters.time_integration(Snoparam);

    Ndays = 10*365
    PZB = ShallowWaters.Parameter(T=T,
        output=true,
        output_dt=168,
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

    offlineweights = load_object("./tuned_weights/result_offline_150iterations_reluactivation_111925.jld2").solution
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
    Ndays = 10*365
    Ponline = ShallowWaters.Parameter(T=T,
        output=true,
        output_dt=168,
        # output_dt=24,
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
        diffusion="constant",        # this is the only new parameter to be adjusted in the new spinups
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
        Ndays=Ndays,
        initial_cond="ncfile",
        initpath="./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_10years_weeklysaves"
    );

    Sonline = ShallowWaters.model_setup(Ponline);

    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_3-15-30-40-50-60-80-85daystart_3dayoptimization_initialweights20daystate_30iterations.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_madnlp_states_5dayoptimization_startfrom1daystate_50iterations_geluactivation.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_state_10dayoptimization_startfrom5daystate_noeta_30iterations.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_state_30dayoptimzation_startfrom20day_constantdissipation_6iterations_8hourdata_200maxhistory_fixedcfl.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_state_20dayoptimzation_startfrom10day_constantdissipation_10iterations_8hourdata_200maxhistory_fixedcfl.jld2").solution

    # onlineweights = load_object("./result_multistate_5-20-35-50-65-75daystart_10dayoptimization_initialweights20daystate_fixedcfl_15iterations_constdissipation.jld2").solution
    # onlineweights = load_object("./result_multistate_5-25-45-65daystart_20dayoptimization_initialweights20daystate_fixedcfl_15iterations_constdissipation.jld2").solution
    # onlineweights = load_object("./result_multistate_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86daystart_3dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1-3-4-6-8-10-13-15-18-23-28-30-33-35-38-41-44-46-48-51-52-53-55-58-60-63-64-65-68-73-78-83-86-88-89daystart_1dayoptimizationinitialweightsmulti3daystate_20iterations.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1:2:89initdaystart_1dayoptimization_initialweightsmulti3daystate_fewerinitconds_15iterations.jld2").solution;
    onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-6-8-10-13-15-18-23-28-33-38-41-44-48-51-53-58-63-65-68-73-78-83-86-88daystart_2dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution;
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
    # Sonline.constants.cD = onlineweights[end]

    Sonline.Prog.u .= copy(initial_cond[1]);
    Sonline.Prog.v .= copy(initial_cond[2]);
    Sonline.Prog.η .= copy(initial_cond[3]);

    P = ShallowWaters.time_integration(Sonline);

    # name for run
    # result_online_state_pluscD_weights_constdiffusion_20dayoptimization_startfrom20day_3years_dailysaves

end

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

    uofflinegelu = ncread("./results/128_offlineparam_postspinup_cginitcond_3days_8hoursaves/u.nc", "u");
    vofflinegelu = ncread("./results/128_offlineparam_postspinup_cginitcond_3days_8hoursaves/v.nc", "v");
    etaofflinegelu = ncread("./results/128_offlineparam_postspinup_cginitcond_3days_8hoursaves/eta.nc", "eta");

    uofflinerelu = ncread("./results/128_offlineparam_postspinup_cginitcond_3days_relu_8hoursaves/u.nc", "u");
    vofflinerelu = ncread("./results/128_offlineparam_postspinup_cginitcond_3days_relu_8hoursaves/v.nc", "v");
    etaofflinerelu = ncread("./results/128_offlineparam_postspinup_cginitcond_3days_relu_8hoursaves/eta.nc", "eta");

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

    unoparam10 = ncread("./dissipation_constant/spinup_files/128_noparam_10year_postspinup_weeklysaves/u.nc", "u");
    vnoparam10 = ncread("./dissipation_constant/spinup_files/128_noparam_10year_postspinup_weeklysaves/v.nc", "v");
    etanoparam10 = ncread("./dissipation_constant/spinup_files/128_noparam_10year_postspinup_weeklysaves/eta.nc", "eta");

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

    u1daystaterelu = ncread("./results/128_online_reluactivation_stateweights_1dayoptimization_startfromoffline_madnlp_30days_8hoursaves/u.nc", "u");
    v1daystaterelu = ncread("./results/128_online_reluactivation_stateweights_1dayoptimization_startfromoffline_madnlp_30days_8hoursaves/v.nc", "v");
    eta1daystaterelu = ncread("./results/128_online_reluactivation_stateweights_1dayoptimization_startfromoffline_madnlp_30days_8hoursaves/eta.nc", "eta");

    u5daystaterelu = ncread("./results/128_online_reluactivation_stateweights_5dayoptimization_startfrom1daystate_madnlp_30days_8hoursaves/u.nc", "u");
    v5daystaterelu = ncread("./results/128_online_reluactivation_stateweights_5dayoptimization_startfrom1daystate_madnlp_30days_8hoursaves/v.nc", "v");
    eta5daystaterelu = ncread("./results/128_online_reluactivation_stateweights_5dayoptimization_startfrom1daystate_madnlp_30days_8hoursaves/eta.nc", "eta");

    ukespecpd1dayrelu = ncread("./results/128_online_reluactivation_kespecpdweights_1dayoptimization_startfrom1daystate_madnlp_30days_8hoursaves/u.nc", "u");
    vkespecpd1dayrelu = ncread("./results/128_online_reluactivation_kespecpdweights_1dayoptimization_startfrom1daystate_madnlp_30days_8hoursaves/v.nc", "v");
    etakespecpd1dayrelu = ncread("./results/128_online_reluactivation_kespecpdweights_1dayoptimization_startfrom1daystate_madnlp_30days_8hoursaves/eta.nc", "eta");

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

function offline_plots()

    # nc files
    # u_zb, uhr
    # jld2 files (my save states function)
    # states_noparam, states_nn (untrained), states_trainednn_pd, states_trainednn_pd65, states_trainednn_kespec, states_trainednn_states

    # to get coarse-grained states

    # t is timestep, and I saved every 8 hours up to 30 days
    # this means t can be anything between 1 (the initial condition) and 91 (the final step after 30 days)

    ###################################################################################

     # for showing offline instability
    t = 10
    fig = Figure(size=(950, 250), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etaofflinegelu[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\eta(3 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(etaofflinegelu[:,:,t])),maximum(abs.(etaofflinegelu[:,:,t])))
    );
    Colorbar(fig[1,2], hm1, label="m")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uofflinegelu[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u(3 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(uofflinegelu[:,:,t])),maximum(abs.(uofflinegelu[:,:,t])))
    );
    Colorbar(fig[1,4], hm2, label="m/s")

    ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    vofflinegelu[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"v(3 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(vofflinegelu[:,:,t])),maximum(abs.(vofflinegelu[:,:,t])))
    );
    Colorbar(fig[1,6], hm3, label="m/s")

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

end

function prognostic_plots()

    # Prognostic variables #############################################################

    # high versus low resolution eta
    t = 1
    fig = Figure(size=(775, 300), fontsize=15);
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uhr[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\mathbf{\eta}(3650 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(etahr[:,:,t])),maximum(abs.(etahr[:,:,t])))
    );
    Colorbar(fig[1,2], hm1, label="m")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uhrcg[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\overline{\mathbf{\eta}}(3650 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(etahrcg[:,:,t])),maximum(abs.(etahrcg[:,:,t])))
    );
    Colorbar(fig[1,4], hm2,label="m")

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
    t = 1096
    fig = Figure(size=(700, 550), fontsize=15);
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uhrcg[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\overline{u}(30 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[1,2], hm1, label="m/s")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    unoparam[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u(30 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[1,4], hm2, label="m/s")

    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uzb[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u_{\text{ZB20}}(30 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[2,2], hm3, label="m/s")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    umulti10[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u_{multi20}(30 \; \text{days}, x, y)"),
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

    # looking at u fields to see if additional state optimization helped
    # just u fields
    t = 7*52
    fig = Figure(size=(700, 550), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    u5s[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u_5(2 \; \text{years}, x, y)"),
    # colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    colorrange=(-maximum(abs.(uzb[:,:,t])),maximum(abs.(uzb[:,:,t])))
    );
    Colorbar(fig[1,2], hm1, label="m/s")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    u10s[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u_{10}(2 \; \text{years}, x, y)"),
    # colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    colorrange=(-maximum(abs.(uzb[:,:,t])),maximum(abs.(uzb[:,:,t])))
    );
    Colorbar(fig[1,4], hm2, label="m/s")

    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    u20s[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u_{20}(2 \; \text{years}, x, y)"),
    # colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    colorrange=(-maximum(abs.(uzb[:,:,t])),maximum(abs.(uzb[:,:,t])))    
    );
    Colorbar(fig[2,2], hm3, label="m/s")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    u30s[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u_{30}(2 \; \text{years}, x, y)"),
    # colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    colorrange=(-maximum(abs.(uzb[:,:,t])),maximum(abs.(uzb[:,:,t])))
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

    # just eta fields
    t = 90
    fig = Figure(size=(700, 525), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etahrcg[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\overline{\eta}(90 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(etahrcg[:,:,t])),maximum(abs.(etahrcg[:,:,t])))
    );
    Colorbar(fig[1,2], hm1, label="m")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etazb[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\eta_{\text{ZB20}}(90 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(etahrcg[:,:,t])),maximum(abs.(etahrcg[:,:,t])))
    );
    Colorbar(fig[1,4], hm2, label="m")

    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etamulti3_new[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\eta_{\text{multi}3}(90 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(etahrcg[:,:,t])),maximum(abs.(etahrcg[:,:,t])))
    );
    Colorbar(fig[2,2], hm3, label="m")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    eta30s[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\eta_{30}(90 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(etahrcg[:,:,t])),maximum(abs.(etahrcg[:,:,t])))
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

    # time-averaged eta fields
    fig = Figure(size=(1040, 520), fontsize=15);
    fig.layout.alignmode = Outside();
    Label(
        fig[0, 3],
        "10-year averaged sea-surface height",
        fontsize = 20,
        tellwidth = false
    )

    # change this if you want ten year or 3 year
    index = vcat(1:7:1096, 1097:1461)
    total = length(index)
    timeavg_hr = sum(etahrcgall[:,:,index],dims=3)[:,:,1] ./ total

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    timeavg_hr,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Coarse-grained 3.75 km"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[1,2], hm1, label="m")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(etanoparam10, dims=3)[:,:,1] ./ total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="No closure"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[1,4], hm1, label="m")

    ax2, hm2 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(etazb10, dims=3)[:,:,1] ./ total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="ZB20"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[1,6], hm1, label="m")

    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(etamulti210, dims=3)[:,:,1] ./ total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 2 day"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[2,2], hm3, label="m")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(etamulti3more10, dims=3)[:,:,1] ./ total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 3 day"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[2,4], hm4, label="m")

    ax4, hm4 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(etamulti1010, dims=3)[:,:,1] ./ total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 10 day"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[2,6], hm1, label="m")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # absolute difference in the time averaged ssh fields
    index = vcat(1:7:1096, 1097:1461)
    total = length(index)
    timeavg_hr = sum(etahrcgall[:,:,index],dims=3)[:,:,1] ./ total

    fig = Figure(size=(1040, 520), fontsize=15);
    fig.layout.alignmode = Outside();
    Label(
        fig[0, 3],
        L"Absolute difference between 10-year averaged $\overline{\eta}$ and parameterized models",
        fontsize = 20,
        tellwidth = false
    )

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    timeavg_hr,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Coarse-grained 3.75 km"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[1,2], hm1, label="m")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    abs.(sum(etanoparam10, dims=3)[:,:,1] ./ total .- timeavg_hr),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="No closure"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[1,4], hm1, label="m")

    ax2, hm2 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    abs.(sum(etazb10, dims=3)[:,:,1] ./ total .- timeavg_hr),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="ZB20"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[1,6], hm1, label="m")

    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    abs.(sum(etamulti210, dims=3)[:,:,1] ./ total .- timeavg_hr),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 2 day"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[2,2], hm3, label="m")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    abs.(sum(etamulti3more10, dims=3)[:,:,1] ./ total .- timeavg_hr),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 3 day"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[2,4], hm4, label="m")

    ax4, hm4 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    abs.(sum(etamulti1010, dims=3)[:,:,1] ./ total .- timeavg_hr),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 10 day"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[2,6], hm1, label="m")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # 1 + 5 state optimization versus ZB20

    t = 182
    fig = Figure(size=(700, 550), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    u5daystategelu[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u_{1+5}(90 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(uhrcg[:,:,31])),maximum(abs.(uhrcg[:,:,31])))
    );
    Colorbar(fig[1,2], hm1, label="m")

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
    t = [46, 91, 136, 273]
    fig = Figure(size=(700, 525), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ufourier[:,:,t[1]],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u_{1}(1 \text{ day}, x, y)"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t[2]])),maximum(abs.(uhrcg[:,:,t[2]])))
    );
    Colorbar(fig[1,2], hm1, label="m/s")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ufourier[:,:,t[2]],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u_{1}(3 \text{ days}, x, y)"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t[2]])),maximum(abs.(uhrcg[:,:,t[2]])))
    );
    Colorbar(fig[1,4], hm2, label="m/s")

    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ufourier[:,:,t[3]],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u_{1}(15 \text{ days}, x, y)"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t[2]])),maximum(abs.(uhrcg[:,:,t[2]])))
    );
    Colorbar(fig[2,2], hm3, label="m/s")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ufourier[:,:,t[4]],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u_{1}(30 \text{ days}, x, y)"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t[2]])),maximum(abs.(uhrcg[:,:,t[2]])))
    );
    Colorbar(fig[2,4], hm4, label="m/s")

    # ax1, hm1 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # u5daystategelu[:,:,t[1]],
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title=L"u_{1+5}(1 \text{ day}, x, y)"),
    # colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    # );
    # Colorbar(fig[2,2], hm1, label="m/s")

    # ax2, hm2 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # u5daystategelu[:,:,t[2]],
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title=L"u_{1+5}(3 \text{ days}, x, y)"),
    # colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    # );
    # Colorbar(fig[2,4], hm2, label="m/s")

    # ax3, hm3 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # u5daystategelu[:,:,t[3]],
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title=L"u_{1+5}(15 \text{ days}, x, y)"),
    # colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    # );
    # Colorbar(fig[2,6], hm3, label="m/s")

    # ax4, hm4 = heatmap(fig[2,7], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # u5daystategelu[:,:,t[4]],
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title=L"u_{1+5}(30 \text{ days}, x, y)"),
    # colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    # );
    # Colorbar(fig[2,8], hm4, label="m/s")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[2, 1] = GridLayout()
    gd = fig[2, 3] = GridLayout()
    # ge = fig[2, 1] = GridLayout()
    # gf = fig[2, 3] = GridLayout()
    # gg = fig[2, 5] = GridLayout()
    # gh = fig[2, 7] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)"], [ga, gb, gc, gd])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

end

function computing_vorticity()

    # ∂x!(dudx, u)
    # ∂y!(dudy, u)

    # ∂x!(dvdx, v)
    # ∂y!(dvdy, v)

    P = ShallowWaters.Parameter(T=T,
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
        zb_forcing_dissipation=false,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=2,
        initial_cond="rest"
    );
    S = ShallowWaters.model_setup(P);

    Phr = ShallowWaters.Parameter(T=T,
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
        N=1,
        α=2,
        nx=1024,
        Ndays=2,
        initial_cond="rest"
    );
    Shr = ShallowWaters.model_setup(Phr);


    # high-resolution vorticity
    ζhr = zeros(1025,1025,1461)
    dvdx = zeros(1027,1027)
    dudy = zeros(1027,1027)
    inv_scale = 1 / Shr.constants.scale
    for t in 1:1461
        uhr_, vhr_, _ = ShallowWaters.add_halo(uhrall[:,:,t],vhrall[:,:,t],etahrall[:,:,t],Shr)

        ShallowWaters.∂x!(dvdx, vhr_ .* inv_scale)
        ShallowWaters.∂y!(dudy, uhr_ .* inv_scale)

        @views @inbounds ζhr[:,:,t] .= (dvdx[2:end-1, 2:end-1] .- dudy[2:end-1, 2:end-1]) ./ 3750
    end

    # 30 km relative vorticity probability distribution over the first three years
    ζzb = zeros(129,129,1096);
    ζnoparam = zeros(129,129,1096);
    ζmulti3more = zeros(129,129,1096);
    ζmulti2 = zeros(129,129,1096);
    ζmulti1more = zeros(129,129,1096);
    ζmulti1 = zeros(129,129,1096);
    ζ20s = zeros(129,129,1096);
    ζ30s = zeros(129,129,1096);
    ζ10s = zeros(129,129,1096);
    ζmulti3 = zeros(129,129,1096);
    ζmulti10 = zeros(129,129,1096);
    ζmulti20 = zeros(129,129,1096);

    dvdx = zeros(131,131)
    dudy = zeros(131,131)
    inv_scale = 1 / S.constants.scale
    for t = 1:1096

        umulti2_, vmulti2_, _ = ShallowWaters.add_halo(umulti2[:,:,t], vmulti2[:,:,t], etamulti2[:,:,t], S);
        ShallowWaters.∂x!(dvdx, vmulti2_ .* inv_scale)
        ShallowWaters.∂y!(dudy, umulti2_ .* inv_scale)
        @views @inbounds ζmulti2[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        umulti3_, vmulti3_, _ = ShallowWaters.add_halo(umulti3[:,:,t], vmulti3[:,:,t], etamulti3[:,:,t], S);
        ShallowWaters.∂x!(dvdx, vmulti3_ .* inv_scale)
        ShallowWaters.∂y!(dudy, umulti3_ .* inv_scale)
        @views @inbounds ζmulti3[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        umulti3more_, vmulti3more_, _ = ShallowWaters.add_halo(umulti3more[:,:,t], vmulti3more[:,:,t], etamulti3more[:,:,t], S);
        ShallowWaters.∂x!(dvdx, vmulti3more_ .* inv_scale)
        ShallowWaters.∂y!(dudy, umulti3more_ .* inv_scale)
        @views @inbounds ζmulti3more[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        umulti10_, vmulti10_, _ = ShallowWaters.add_halo(umulti10[:,:,t], vmulti10[:,:,t], etamulti10[:,:,t], S);
        ShallowWaters.∂x!(dvdx, vmulti10_ .* inv_scale)
        ShallowWaters.∂y!(dudy, umulti10_ .* inv_scale)
        @views @inbounds ζmulti10[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        umulti20_, vmulti20_, _ = ShallowWaters.add_halo(umulti20[:,:,t], vmulti20[:,:,t], etamulti20[:,:,t], S);
        ShallowWaters.∂x!(dvdx, vmulti20_ .* inv_scale)
        ShallowWaters.∂y!(dudy, umulti20_ .* inv_scale)
        @views @inbounds ζmulti20[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        uzb_, vzb_, _ = ShallowWaters.add_halo(uzb[:,:,t],vzb[:,:,t],etazb[:,:,t],S);
        ShallowWaters.∂x!(dvdx, vzb_ .* inv_scale)
        ShallowWaters.∂y!(dudy, uzb_ .* inv_scale)
        @views @inbounds ζzb[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        unoparam_, vnoparam_, _ = ShallowWaters.add_halo(unoparam[:,:,t],vnoparam[:,:,t],etanoparam[:,:,t],S);
        ShallowWaters.∂x!(dvdx, vnoparam_ .* inv_scale)
        ShallowWaters.∂y!(dudy, unoparam_ .* inv_scale)
        @views @inbounds ζnoparam[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        u10s_, v10s_, _ = ShallowWaters.add_halo(u10s[:,:,t],v10s[:,:,t],eta10s[:,:,t], S);
        ShallowWaters.∂x!(dvdx, v10s_ .* inv_scale)
        ShallowWaters.∂y!(dudy, u10s_ .* inv_scale)
        @views @inbounds ζ10s[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        u20s_, v20s_, _ = ShallowWaters.add_halo(u20s[:,:,t],v20s[:,:,t],eta20s[:,:,t], S);
        ShallowWaters.∂x!(dvdx, v20s_ .* inv_scale)
        ShallowWaters.∂y!(dudy, u20s_ .* inv_scale)
        @views @inbounds ζ20s[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        u30s_, v30s_, _ = ShallowWaters.add_halo(u30s[:,:,t],v30s[:,:,t],eta30s[:,:,t], S);
        ShallowWaters.∂x!(dvdx, v30s_ .* inv_scale)
        ShallowWaters.∂y!(dudy, u30s_ .* inv_scale)
        @views @inbounds ζ30s[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

    end

    # Relative vorticity plots in the last seven years
    ζzb10 = zeros(129,129,522);
    ζnoparam10 = zeros(129,129,522);
    ζmulti3more10 = zeros(129,129,522);
    ζmulti210 = zeros(129,129,522);
    ζmulti1more10 = zeros(129,129,522);
    ζ20s10 = zeros(129,129,522);
    ζ30s10 = zeros(129,129,522);
    ζ10s10 = zeros(129,129,522);
    ζmulti310 = zeros(129,129,522);
    ζmulti1010 = zeros(129,129,522);
    ζmulti2010 = zeros(129,129,522);
    dvdx = zeros(131,131)
    dudy = zeros(131,131)
    inv_scale = 1 / S.constants.scale
    for t = 1:522

        umulti210_, vmulti210_, _ = ShallowWaters.add_halo(umulti210[:,:,t], vmulti210[:,:,t], etamulti210[:,:,t], S);
        ShallowWaters.∂x!(dvdx, vmulti210_ .* inv_scale)
        ShallowWaters.∂y!(dudy, umulti210_ .* inv_scale)
        @views @inbounds ζmulti210[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        umulti310_, vmulti310_, _ = ShallowWaters.add_halo(umulti310[:,:,t], vmulti310[:,:,t], etamulti310[:,:,t], S);
        ShallowWaters.∂x!(dvdx, vmulti310_ .* inv_scale)
        ShallowWaters.∂y!(dudy, umulti310_ .* inv_scale)
        @views @inbounds ζmulti310[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        umulti3more10_, vmulti3more10_, _ = ShallowWaters.add_halo(umulti3more10[:,:,t], vmulti3more10[:,:,t], etamulti3more10[:,:,t], S);
        ShallowWaters.∂x!(dvdx, vmulti3more10_ .* inv_scale)
        ShallowWaters.∂y!(dudy, umulti3more10_ .* inv_scale)
        @views @inbounds ζmulti3more10[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        umulti1010_, vmulti1010_, _ = ShallowWaters.add_halo(umulti1010[:,:,t], vmulti1010[:,:,t], etamulti1010[:,:,t], S);
        ShallowWaters.∂x!(dvdx, vmulti1010_ .* inv_scale)
        ShallowWaters.∂y!(dudy, umulti1010_ .* inv_scale)
        @views @inbounds ζmulti1010[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        umulti2010_, vmulti2010_, _ = ShallowWaters.add_halo(umulti2010[:,:,t], vmulti2010[:,:,t], etamulti2010[:,:,t], S);
        ShallowWaters.∂x!(dvdx, vmulti2010_ .* inv_scale)
        ShallowWaters.∂y!(dudy, umulti2010_ .* inv_scale)
        @views @inbounds ζmulti2010[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        uzb10_, vzb10_, _ = ShallowWaters.add_halo(uzb10[:,:,t],vzb10[:,:,t],etazb10[:,:,t],S);
        ShallowWaters.∂x!(dvdx, vzb10_ .* inv_scale)
        ShallowWaters.∂y!(dudy, uzb10_ .* inv_scale)
        @views @inbounds ζzb10[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        unoparam10_, vnoparam10_, _ = ShallowWaters.add_halo(unoparam10[:,:,t],vnoparam10[:,:,t],etanoparam10[:,:,t],S);
        ShallowWaters.∂x!(dvdx, vnoparam10_ .* inv_scale)
        ShallowWaters.∂y!(dudy, unoparam10_ .* inv_scale)
        @views @inbounds ζnoparam10[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        u10s10_, v10s10_, _ = ShallowWaters.add_halo(u10s10[:,:,t],v10s10[:,:,t],eta10s10[:,:,t], S);
        ShallowWaters.∂x!(dvdx, v10s10_ .* inv_scale)
        ShallowWaters.∂y!(dudy, u10s10_ .* inv_scale)
        @views @inbounds ζ10s10[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        u20s10_, v20s10_, _ = ShallowWaters.add_halo(u20s10[:,:,t],v20s10[:,:,t],eta20s10[:,:,t], S);
        ShallowWaters.∂x!(dvdx, v20s10_ .* inv_scale)
        ShallowWaters.∂y!(dudy, u20s10_ .* inv_scale)
        @views @inbounds ζ20s10[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

        u30s10_, v30s10_, _ = ShallowWaters.add_halo(u30s10[:,:,t],v30s10[:,:,t],eta30s10[:,:,t], S);
        ShallowWaters.∂x!(dvdx, v30s10_ .* inv_scale)
        ShallowWaters.∂y!(dudy, u30s10_ .* inv_scale)
        @views @inbounds ζ30s10[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000

    end

    ζhrcg = zeros(129,129,1461);
    dvdx = zeros(131,131)
    dudy = zeros(131,131)
    inv_scale = 1 / S.constants.scale
    for t = 1:1461
        uhrcg_, vhrcg_, _ = ShallowWaters.add_halo(uhrcgall[:,:,t],vhrcgall[:,:,t],etahrcgall[:,:,t],S);
        ShallowWaters.∂x!(dvdx, vhrcg_ .* inv_scale)
        ShallowWaters.∂y!(dudy, uhrcg_ .* inv_scale)
        @views @inbounds ζhrcg[:,:,t] .= (dvdx[2:end-1,2:end-1] .- dudy[2:end-1,2:end-1]) ./ 30000
    end

end

function vorticity_plots()

    # three year 
    hr = kde(vec(ζhr[:,:,1:1096]); npoints=16384);
    hrcg = kde(vec(ζhrcg[:,:,1:1096]));
    zb = kde(vec(ζzb));
    noparam = kde(vec(ζnoparam));
    multi2 = kde(vec( ζmulti2));
    multi3more = kde(vec(ζmulti3more));
    multi20 = kde(vec(ζmulti20));
    thirty = kde(vec(ζ30s));

    fig = Figure(size = (650, 420));
    ax = Axis(
        fig[1, 1],
        xlabel = "Relative Vorticity",
        ylabel = "Probability Density",
        title = "Three-year relative vorticity probability density"
        # yscale = log10
    )
    lines!(ax, hr.x, hr.density, label="3.75 km resolution", color=:black)
    lines!(ax, hrcg.x, hrcg.density, label="Coarse-grained 3.75 km resolution",color=:gray)
    lines!(ax, zb.x, zb.density, label="ZB20")
    lines!(ax, noparam.x, noparam.density, label="No closure")
    lines!(ax, multi2.x, multi2.density, label="Ensemble 2 day")
    lines!(ax, multi3more.x, multi3more.density, label="Ensemble 3 day")
    lines!(ax, multi20.x, multi20.density, label="Ensemble 20 day")#,linestyle=:dash)
    lines!(ax, thirty.x, thirty.density, label="30 day")#, linestyle=:dashdot)

    Legend(fig[1,2], ax)
    xlims!(ax, -0.000007, 0.000007)
    # ylims!(ax, 1, 10^4)

    # 10 year figure

    hr10 = kde(vec(cat(ζhr[:,:,1:7:1096], ζhr[:,:,1097:end]; dims=3)); npoints=16384);
    hrcg10 = kde(vec(cat(ζhrcg[:,:,1:7:1096], ζhrcg[:,:,1097:end]; dims=3)));
    zb10 = kde(vec(ζzb10));
    noparam10 = kde(vec(ζnoparam10));
    multi210 = kde(vec(ζmulti210));
    multi3more10 = kde(vec(ζmulti3more10));
    multi2010 = kde(vec(ζmulti2010));
    thirty10 = kde(vec(ζ30s10));

    fig = Figure(size = (650, 420));
    ax = Axis(
        fig[1, 1],
        xlabel = "Relative Vorticity",
        ylabel = "Probability Density",
        title = "Ten-year relative vorticity probability density"
        # yscale = log10
    )
    lines!(ax, hr10.x, hr10.density, label="3.75 km resolution", color=:black)
    lines!(ax, hrcg10.x, hrcg10.density, label="Coarse-grained 3.75 km resolution",color=:gray)
    lines!(ax, zb10.x, zb10.density, label="ZB20")
    lines!(ax, noparam10.x, noparam10.density, label="No closure")
    lines!(ax, multi210.x, multi210.density, label="Ensemble 2 day")
    lines!(ax, multi3more10.x, multi3more10.density, label="Ensemble 3 day")
    lines!(ax, multi2010.x, multi2010.density, label="Ensemble 20 day")#,linestyle=:dash)
    lines!(ax, thirty10.x, thirty10.density, label="30 day")#, linestyle=:dashdot)
    Legend(fig[1,2], ax)
    xlims!(ax, -0.000007, 0.000007)
    # ylims!(ax, 1, 10^4)

    # plot of the vorticity itself
    fig = Figure(size=(1040, 400), fontsize=15);
    t = 1096
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 1025),
    LinRange(0, 3840, 1025),
    ζhrch[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Fitlered, coarse-grained 3.75 km resolution vorticity"),# title=L"\overline{\zeta}(3 \; \text{years}, x, y), \; 3.75 \; \text{km}"),
    colorrange=(-maximum(abs.(ζhrcg[:,:,t])),maximum(abs.(ζhcgr[:,:,t])))
    );
    Colorbar(fig[1,2], hm1, label="1/s")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζnoparam[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="30 km resolution vorticity"),# title=L"\zeta(3 \; \text{years}, x, y), "),
    colorrange=(-maximum(abs.(ζhrcg[:,:,t])),maximum(abs.(ζhrcg[:,:,t])))
    );
    Colorbar(fig[1,4], hm1, label="1/s")

    ax2, hm2 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζzb,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{\text{ZB20}}(3 \; \text{years}, x, y)"),
    colorrange=(-maximum(abs.(ζhrcg)),maximum(abs.(ζhrcg)))
    );
    Colorbar(fig[1,6], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti2,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{\text{multi}2}(3 \; \text{years}, x, y)"),
    colorrange=(-maximum(abs.(ζhrcg)),maximum(abs.(ζhrcg)))
    );
    Colorbar(fig[2,2], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti3more,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{\text{multi}3}(3 \; \text{years}, x, y)"),
    colorrange=(-maximum(abs.(ζhrcg)),maximum(abs.(ζhrcg)))
    );
    Colorbar(fig[2,4], hm1, label="1/s")

    ax3, hm3 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti10,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{\text{multi}10}(3 \; \text{years}, x, y)"),
    colorrange=(-maximum(abs.(ζhrcg)),maximum(abs.(ζhrcg)))
    );
    Colorbar(fig[2,6], hm1, label="1/s")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # examining the vorticity of the results that diverge (during first three years)
    fig = Figure(size=(950, 475), fontsize=15);
    t = 1096
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζzb[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{ZB20}(2000 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(ζhrcg[:,:,t])),maximum(abs.(ζhrcg[:,:,t])))
    );
    Colorbar(fig[1,2], hm1, label="1/s")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζnoparam[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta(2000 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(ζhrcg[:,:,t])),maximum(abs.(ζhrcg[:,:,t])))
    );
    Colorbar(fig[1,4], hm1, label="1/s")

    ax2, hm2 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζ10s[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{10}(2000 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(ζhrcg[:,:,t])),maximum(abs.(ζhrcg[:,:,t])))
    );
    Colorbar(fig[1,6], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζ20s[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{20}(2000 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(ζhrcg[:,:,t])),maximum(abs.(ζhrcg[:,:,t])))
    );
    Colorbar(fig[2,2], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζ30s[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{30}(2000 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(ζhrcg[:,:,t])),maximum(abs.(ζhrcg[:,:,t])))
    );
    Colorbar(fig[2,4], hm1, label="1/s")

    ax3, hm3 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti20[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{\text{multi}20}(2000 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(ζhrcg[:,:,t])),maximum(abs.(ζhrcg[:,:,t])))
    );
    Colorbar(fig[2,6], hm1, label="1/s")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # vorticity plots in the last 7 years

    fig = Figure(size=(1040, 520), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζzb,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\overline{\zeta}(3 \; \text{years}, x, y)"),
    colorrange=(-maximum(abs.(ζzb)),maximum(abs.(ζzb)))
    );
    Colorbar(fig[1,2], hm1, label="1/s")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζnoparam,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta(3 \; \text{years}, x, y)"),
    colorrange=(-maximum(abs.(ζzb)),maximum(abs.(ζzb)))
    );
    Colorbar(fig[1,4], hm1, label="1/s")

    ax2, hm2 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζ10s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{\text{ZB20}}(3 \; \text{years}, x, y)"),
    colorrange=(-maximum(abs.(ζzb)),maximum(abs.(ζzb)))
    );
    Colorbar(fig[1,6], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti2,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{\text{multi}2}(3 \; \text{years}, x, y)"),
    colorrange=(-maximum(abs.(ζzb)),maximum(abs.(ζzb)))
    );
    Colorbar(fig[2,2], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti3more,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{\text{multi}3}(3 \; \text{years}, x, y)"),
    colorrange=(-maximum(abs.(ζzb)),maximum(abs.(ζzb)))
    );
    Colorbar(fig[2,4], hm1, label="1/s")

    ax3, hm3 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti10,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{\text{multi}10}(3 \; \text{years}, x, y)"),
    colorrange=(-maximum(abs.(ζzb)),maximum(abs.(ζzb)))
    );
    Colorbar(fig[2,6], hm1, label="1/s")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # examining the vorticity of the results that diverge over the last seven years
    fig = Figure(size=(950, 475), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζzb[:,:,522],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{ZB20}(2000 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(ζzb[:,:,522])),maximum(abs.(ζzb[:,:,522])))
    );
    Colorbar(fig[1,2], hm1, label="1/s")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζnoparam[:,:,522],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta(2000 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(ζzb[:,:,522])),maximum(abs.(ζzb[:,:,522])))
    );
    Colorbar(fig[1,4], hm1, label="1/s")

    # ax2, hm2 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # ζ10s,
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title=L"\zeta_{10}(2000 \; \text{days}, x, y)"),
    # colorrange=(-maximum(abs.(ζzb)),maximum(abs.(ζzb)))
    # );
    # Colorbar(fig[1,6], hm1, label="1/s")

    # ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # ζ20s,
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title=L"\zeta_{20}(2000 \; \text{days}, x, y)"),
    # colorrange=(-maximum(abs.(ζzb)),maximum(abs.(ζzb)))
    # );
    # Colorbar(fig[2,2], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζ30s[:,:,522],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{30}(2000 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(ζzb[:,:,522])),maximum(abs.(ζzb[:,:,522])))
    );
    Colorbar(fig[2,4], hm1, label="1/s")

    ax3, hm3 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti20[:,:,522],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\zeta_{\text{multi}20}(2000 \; \text{days}, x, y)"),
    colorrange=(-maximum(abs.(ζzb[:,:,522])),maximum(abs.(ζzb[:,:,522])))
    );
    Colorbar(fig[2,6], hm1, label="1/s")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # plot of the vorticity itself, for poster with different labels
    fig = Figure(size=(1098, 520), fontsize=15);

    Label(
        fig[0, 3],
        L"\zeta(3 \text{ years}, x, y)",
        fontsize = 20,
        tellwidth = false
    )

    t = 1096
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 1024),
    LinRange(0, 3840, 1024),
    ζhrcg[:,:, t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Coarse-grained 3.75 km"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    Colorbar(fig[1,2], hm1, label="1/s")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζnoparam[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="No closure"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    Colorbar(fig[1,4], hm1, label="1/s")

    ax2, hm2 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζzb[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="ZB20"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    Colorbar(fig[1,6], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti2[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 2 day"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    Colorbar(fig[2,2], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti3more[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 3 day"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    Colorbar(fig[2,4], hm1, label="1/s")

    ax3, hm3 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti10[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 10 day"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    Colorbar(fig[2,6], hm1, label="1/s")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

end

function energy_plots()

    # Energy ############################################################################

    # high-resolution versus coarse-grained high resolution energy
    t = 1098
    fig = Figure(size=(800, 350), fontsize=15);
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uhr[:,1:end-1,t].^2 .+ vhr[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="3.75 km resolution E(10 days, x, y)"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[1,2], hm1, label=L"(m/s)^2")

    ax1, hm1 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Coarse-grained E(10 days, x, y)"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[1,4], hm1, label=L"(m/s)^2")

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
    t = 1096
    fig = Figure(size=(800, 650), fontsize=15);
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Coarsened, filtered 3.75 km"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[1,2], hm1)

    ax1, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (unoparam[:,1:end-1,t].^2 .+ vnoparam[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="No closure"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[1,4], hm1)

    ax1, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uzb[:,1:end-1,t].^2 .+ vzb[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="ZB closure"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[2,2], hm1)

    ax1, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (u30s[:,1:end-1,t].^2 .+ v30s[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Online, 20 day state"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[2,4], hm1, label=L"(m/s)^2")

    # ax1, hm5 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # (uofflinegelu[:,1:end-1,10].^2 .+ vofflinegelu[1:end-1,:,10].^2),
    # colormap=:amp,
    # axis=(xlabel="km", ylabel="km", title="30 km resolution E, offline closure after 3 days"),
    # colorrange=(0,
    # maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    # );
    # Colorbar(fig[3,2], hm1)

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

    # same as the above but without the coarse-grained energy
    # cg, zb, nn, no param
    t = 200
    fig = Figure(size=(900, 800), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(unoparam[:,:,t]).^2 .+ sum(vnoparam[1:end-1,:,t]).^2,
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30km resolution E(10 days, x, y), no closure"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[1,2], hm1)

    ax1, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uzb[:,1:end-1,t].^2 .+ vzb[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30 km resolution E(10 days, x, y), ZB closure"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[1,4], hm2)

    ax1, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (u20s[:,1:end-1,t].^2 .+ v20s[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30 km resolution E(10 days, x, y), online closure"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[2,2], hm3)

    ax1, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (u30s[:,1:end-1,t].^2 .+ v30s[1:end-1,:,t].^2),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="30 km resolution E(3 days, x, y), offline closure"),
    colorrange=(0,
    maximum(abs.(uhrcg[:,1:end-1,t].^2 .+ vhrcg[1:end-1,:,t].^2)))
    );
    Colorbar(fig[2,4], hm5)

    # looking at the energy in the models past the three year mark
    t = 286
    fig = Figure(size=(900, 800), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    unoparam[:,1:end-1,t].^2 .+ vnoparam[1:end-1,:,t].^2,
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="E(2000 days, x, y), no closure"),
    colorrange=(0,
        maximum(uzb[:,1:end-1,t].^2 .+ vzb[1:end-1,:,t].^2))
    );
    Colorbar(fig[1,2], hm1)

    ax1, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uzb[:,1:end-1,t].^2 .+ vzb[1:end-1,:,t].^2,
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="E(2000 days, x, y), ZB20 closure"),
    colorrange=(0,
        maximum(uzb[:,1:end-1,t].^2 .+ vzb[1:end-1,:,t].^2))
    );
    Colorbar(fig[1,4], hm2)

    ax1, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    umulti3more[:,1:end-1,t].^2 .+ vmulti3more[1:end-1,:,t].^2,
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="E(2000 days, x, y), batched 3 day day closure"),
    colorrange=(0,
        maximum(uzb[:,1:end-1,t].^2 .+ vzb[1:end-1,:,t].^2))
    );
    Colorbar(fig[2,2], hm3)

    ax1, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    umulti10[:,1:end-1,t].^2 .+ vmulti10[1:end-1,:,t].^2,
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="E(2000 days, x, y), batched 20 day closure"),
    colorrange=(0,
        maximum(uzb[:,1:end-1,t].^2 .+ vzb[1:end-1,:,t].^2))
    );
    Colorbar(fig[2,4], hm5)

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

    # spatially averaged energy over integration (3 year integrations on all)

    N = 1096

    zb          = zeros(Float64, N)
    noparam    = zeros(Float64, N)
    fiveday    = zeros(Float64, N)
    tenday     = zeros(Float64, N)
    twentyday  = zeros(Float64, N)
    thirtyday  = zeros(Float64, N)
    multi2      = zeros(Float64, N)
    multi3      = zeros(Float64, N)
    multi3more = zeros(Float64, N)
    multi5      = zeros(Float64, N)
    multi10     = zeros(Float64, N)
    multi20     = zeros(Float64, N)
    for j = 1:1096

        zb[j] = sum(abs2, uzb[:,:,j]) + sum(abs2, vzb[:,:,j])
        noparam[j] = sum(abs2, unoparam[:,:,j]) + sum(abs2, vnoparam[:,:,j])

        fiveday[j] = sum(abs2, u5s[:,:,j])  + sum(abs2, v5s[:,:,j])
        tenday[j] = sum(abs2, u10s[:,:,j]) + sum(abs2, v10s[:,:,j])
        twentyday[j] = sum(abs2, u20s[:,:,j]) + sum(abs2, v20s[:,:,j])
        thirtyday[j] = sum(abs2, u30s[:,:,j]) + sum(abs2, v30s[:,:,j])
        multi2[j] = sum(abs2, umulti2[:,:,j]) + sum(abs2, vmulti2[:,:,j])
        multi3[j] = sum(abs2, umulti3[:,:,j]) + sum(abs2, vmulti3[:,:,j])
        multi3more[j] = sum(abs2, umulti3more[:,:,j]) + sum(abs2, vmulti3more[:,:,j])
        multi5[j] = sum(abs2, umulti5[:,:,j]) + sum(abs2, vmulti5[:,:,j])
        multi10[j] = sum(abs2, umulti10[:,:,j]) + sum(abs2, vmulti10[:,:,j])
        multi20[j] = sum(abs2, umulti20[:,:,j]) + sum(abs2, vmulti20[:,:,j])
        # appendix stuff
        # push!(kespec, sum(ukespec[:,1:end-1,j].^2 .+ vkespec[1:end-1,:,j].^2))
        # push!(hybrid, sum(uhybrid[:,1:end-1,j].^2 .+ vhybrid[1:end-1,:,j].^2))
        # push!(fourier, sum(ufourier[:,1:end-1,j].^2 .+ vfourier[1:end-1,:,j].^2))
        # push!(kespecpd, sum(ukespecpd[:,1:end-1,j].^2 .+ vkespecpd[1:end-1,:,j].^2))
        # push!(relu1day, sum(uonline1dayrelu[:,1:end-1,j].^2 .+ vonline1dayrelu[1:end-1,:,j].^2))
        # push!(relu5day, sum(uonline5dayrelu[:,1:end-1,j].^2 .+ vonline5dayrelu[1:end-1,:,j].^2))
        # push!(reluKEspec, sum(uonlinekespecpdrelu[:,1:end-1,j].^2 .+ vonlinekespecpdrelu[1:end-1,:,j].^2))

    end

    N = 522

    zb10 = zeros(Float64, N)
    noparam10 = zeros(Float64, N)
    fiveday10 = zeros(Float64, N)
    tenday10 = zeros(Float64, N)
    twentyday10 = zeros(Float64, N)
    thirtyday10 = zeros(Float64, N)
    multi210 = zeros(Float64, N)
    multi310 = zeros(Float64, N)
    multi3more10 = zeros(Float64, N)
    multi510 = zeros(Float64, N)
    multi1010 = zeros(Float64, N)
    multi2010 = zeros(Float64, N)
    for j = 1:522

        zb10[j] = sum(abs2, uzb10[:,:,j]) + sum(abs2, vzb10[:,:,j])
        noparam10[j] = sum(abs2, unoparam10[:,:,j]) + sum(abs2, vnoparam10[:,:,j])

        fiveday10[j] = sum(abs2, u5s10[:,:,j]) + sum(abs2, v5s10[:,:,j])
        tenday10[j] = sum(abs2, u10s10[:,:,j]) + sum(abs2, v10s10[:,:,j])
        twentyday10[j] = sum(abs2, u20s10[:,:,j]) + sum(abs2, v20s10[:,:,j])
        thirtyday10[j] = sum(abs2, u30s10[:,:,j]) + sum(abs2, v30s10[:,:,j])

        multi210[j] = sum(abs2, umulti210[:,:,j]) + sum(abs2, vmulti210[:,:,j])
        multi310[j] = sum(abs2, umulti310[:,:,j]) + sum(abs2, vmulti310[:,:,j])
        multi3more10[j] = sum(abs2, umulti3more10[:,:,j]) + sum(abs2, vmulti3more10[:,:,j])
        # multi510[j] = sum(abs2, umulti510[:,:,j]) + sum(abs2, vmulti510[:,:,j])
        multi1010[j] = sum(abs2, umulti1010[:,:,j]) + sum(abs2, vmulti1010[:,:,j])
        multi2010[j] = sum(abs2, umulti2010[:,:,j]) + sum(abs2, vmulti2010[:,:,j])
        # appendix stuff
        # push!(kespec, sum(ukespec10[:,1:end-1,j].^2 .+ vkespec10[1:end-1,:,j].^2))
        # push!(hybrid, sum(uhybrid10[:,1:end-1,j].^2 .+ vhybrid10[1:end-1,:,j].^2))
        # push!(fourier, sum(ufourier10[:,1:end-1,j].^2 .+ vfourier10[1:end-1,:,j].^2))
        # push!(kespecpd, sum(ukespecpd10[:,1:end-1,j].^2 .+ vkespecpd10[1:end-1,:,j].^2))
        # push!(relu1day, sum(uonline1dayrelu10[:,1:end-1,j].^2 .+ vonline1dayrelu10[1:end-1,:,j].^2))
        # push!(relu5day, sum(uonline5dayrelu10[:,1:end-1,j].^2 .+ vonline5dayrelu10[1:end-1,:,j].^2))
        # push!(reluKEspec, sum(uonlinekespecpdrelu10[:,1:end-1,j].^2 .+ vonlinekespecpdrelu10[1:end-1,:,j].^2))

    end

    hrcg = zeros(Float64, 1461)
    for j = 1:1461
            hrcg[j] = sum(abs2, uhrcgall[:,:,j]) + sum(abs2, vhrcgall[:,:,j])
    end

    multi210p = zeros(Float64, 522)
    multi310p = zeros(Float64, 522)
    for j = 1:522
            multi210p[j] = sum(abs2, umulti210more[:,:,j]) + sum(abs2, vmulti210more[:,:,j])
            multi310p[j] = sum(abs2, umulti310more[:,:,j]) + sum(abs2, vmulti310more[:,:,j])
    end

    multi2all = cat(multi210, multi210p[2:end]; dims=1);
    multi3all = cat(multi3more10, multi310p[2:end]; dims=1);

    # 3 year figure
    fig = Figure(size=(1000, 500), fontsize=15);
    ax = Axis(fig[1,1],
            xlabel="Day",
            ylabel="Energy",
            title="Spatially averaged energy over 3 years"
    )
    lines!(ax, LinRange(0, 3*365, 1096),  cghr ./ (128^2), label="Coarse-grained 3.75 km resolution", color=:black)
    lines!(ax, LinRange(0, 3*365, 1096), noparam./ (128^2), label="30 km resolution, no closure", color=:gray)
    lines!(ax, LinRange(0, 3*365, 1096), zb./ (128^2), label="ZB20")
    # lines!(fig[1,1], LinRange(0, 3*365, 1096), fiveday./ (128^2), label="Online closure, 5 day")
    # lines!(ax, LinRange(0, 3*365, 1096), tenday./ (128^2), label="Online closure, 10 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentyday ./ (128^2), label="Online closure, 20 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentydaycD ./ (128^2), label="Online closure, 20 day with BD coeff")
    lines!(ax, LinRange(0, 3*365, 1096), thirtyday ./ (128^2), label="Online closure, 30 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi1 ./ (128^2), label="Online closure, batched 1 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi1_more ./ (128^2), label="Online closure, batched 1 day")
    lines!(ax, LinRange(0, 3*365, 1096), multi2 ./ (128^2), label="Online closure, batched 2 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi3more ./ (128^2), label="Online closure, batched 3 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi5 ./ (128^2), label="Online closure, batched 5 day", color=:mediumorchid)
    # lines!(ax, LinRange(0, 3*365, 1096), multi10 ./ (128^2), label="Online closure, batched 10 day")#, color=:teal)
    lines!(ax, LinRange(0, 3*365, 1096), multi20 ./ (128^2), label="Online closure, batched 20 day", color=:red3)
    Legend(fig[1, 2], ax)

    # 3 year relative error figure
    fig = Figure(size=(1000, 500), fontsize=15);
    
    ax = Axis(fig[1,1],
            xlabel="Day",
            # ylabel="Energy",
            title="Relative error in spatially averaged energy"
    )
    lines!(ax, LinRange(0, 3*365, 1096), (abs.(noparam./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    lines!(ax, LinRange(0, 3*365, 1096), (abs.(zb./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{ZB} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    # lines!(fig[1,1], LinRange(0, 3*365, 1096), fiveday./ (128^2), label="Online closure, 5 day")
    # lines!(ax, LinRange(0, 3*365, 1096), tenday./ (128^2), label="Online closure, 10 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentyday ./ (128^2), label="Online closure, 20 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentydaycD ./ (128^2), label="Online closure, 20 day with BD coeff")
    lines!(ax, LinRange(0, 3*365, 1096), (abs.( thirtyday./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{30} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    # lines!(ax, LinRange(0, 3*365, 1096), multi1 ./ (128^2), label="Online closure, batched 1 day")
    # lines!(ax, LinRange(0, 3*365, 1096), (abs.( multi1./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi1} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    # lines!(ax, LinRange(0, 3*365, 1096), (abs.( multi1_more./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi1} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    lines!(ax, LinRange(0, 3*365, 1096), (abs.( multi2./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi2} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    # lines!(ax, LinRange(0, 3*365, 1096), multi2 ./ (128^2), label="Online closure, batched 2 day")
    lines!(ax, LinRange(0, 3*365, 1096), (abs.( multi3more./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi3} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    # lines!(ax, LinRange(0, 3*365, 1096), multi5 ./ (128^2), label="Online closure, batched 5 day", color=:mediumorchid)
    lines!(ax, LinRange(0, 3*365, 1096), (abs.( multi10./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi10} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")#, color=:teal)
    lines!(ax, LinRange(0, 3*365, 1096), multi20 ./ (128^2), label="Online closure, batched 20 day", color=:red3)
    Legend(fig[1, 2], ax)

    # combined three year figure plus the relative error figure

    fig = Figure(size=(1000, 500), fontsize=15);
    ax = Axis(fig[1,1],
            xlabel="Day",
            ylabel="Energy",
            title="Spatially averaged energy over 3 years"
    )
    lines!(ax, LinRange(0, 3*365, 1096), noparam./ (128^2), label="30 km resolution, no closure")
    lines!(ax, LinRange(0, 3*365, 1096), zb./ (128^2), label="ZB20")
    # lines!(fig[1,1], LinRange(0, 3*365, 1096), fiveday./ (128^2), label="Online closure, 5 day")
    # lines!(ax, LinRange(0, 3*365, 1096), tenday./ (128^2), label="Online closure, 10 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentyday ./ (128^2), label="Online closure, 20 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentydaycD ./ (128^2), label="Online closure, 20 day with BD coeff")
    lines!(ax, LinRange(0, 3*365, 1096), thirtyday ./ (128^2), label="Online closure, 30 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi1 ./ (128^2), label="Online closure, batched 1 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi1_more ./ (128^2), label="Online closure, batched 1 day")
    lines!(ax, LinRange(0, 3*365, 1096), multi2 ./ (128^2), label="Online closure, batched 2 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi3more ./ (128^2), label="Online closure, batched 3 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi5 ./ (128^2), label="Online closure, batched 5 day", color=:mediumorchid)
    lines!(ax, LinRange(0, 3*365, 1096), multi10 ./ (128^2), label="Online closure, batched 10 day")#, color=:teal)
    lines!(ax, LinRange(0, 3*365, 1096), multi20 ./ (128^2), label="Online closure, batched 20 day", color=:red3)
    lines!(ax, LinRange(0, 3*365, 1096),  cghr ./ (128^2), label="Coarse-grained HR", color=:darkorchid)

    Legend(fig[1, 2], ax)

    ax2 = Axis(fig[2,1],
            xlabel="Day",
            # ylabel="Energy",
            title="Relative error in spatially averaged energy"
    )
    lines!(ax2, LinRange(0, 3*365, 1096), (abs.(noparam./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    lines!(ax2, LinRange(0, 3*365, 1096), (abs.(zb./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{ZB} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    # lines!(fig[1,1], LinRange(0, 3*365, 1096), fiveday./ (128^2), label="Online closure, 5 day")
    # lines!(ax, LinRange(0, 3*365, 1096), tenday./ (128^2), label="Online closure, 10 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentyday ./ (128^2), label="Online closure, 20 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentydaycD ./ (128^2), label="Online closure, 20 day with BD coeff")
    lines!(ax2, LinRange(0, 3*365, 1096), (abs.( thirtyday./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{30} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    # lines!(ax, LinRange(0, 3*365, 1096), multi1 ./ (128^2), label="Online closure, batched 1 day")
    # lines!(ax, LinRange(0, 3*365, 1096), (abs.( multi1./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi1} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    # lines!(ax, LinRange(0, 3*365, 1096), (abs.( multi1_more./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi1} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    lines!(ax2, LinRange(0, 3*365, 1096), (abs.( multi2./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi2} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    # lines!(ax, LinRange(0, 3*365, 1096), multi2 ./ (128^2), label="Online closure, batched 2 day")
    lines!(ax2, LinRange(0, 3*365, 1096), (abs.( multi3more./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi3} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    # lines!(ax, LinRange(0, 3*365, 1096), multi5 ./ (128^2), label="Online closure, batched 5 day", color=:mediumorchid)
    lines!(ax2, LinRange(0, 3*365, 1096), (abs.( multi10./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi10} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")#, color=:teal)
    lines!(ax2, LinRange(0, 3*365, 1096), (abs.( multi20./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi20} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}", color=:red3)
    Legend(fig[2, 2], ax2)

    # 10 year figure
    cghr10_forplotting = cat(cghr10[1:7:1096], cghr10[1097:end]; dims=1)
    fig = Figure(size=(1000, 500), fontsize=15);
    ax = Axis(fig[1,1],
            xlabel="Day",
            ylabel="Energy",
            title="Spatially averaged energy over 10 years"
    )
    lines!(ax, LinRange(0, 10*365, 522),  cghr10_forplotting ./ (128^2), label="Filtered, coarse-grained 3.75 km")
    lines!(ax, LinRange(0, 10*365, 522),  zb10 ./ (128^2), label="ZB20")
    lines!(ax, LinRange(0, 10*365, 522), noparam10 ./ (128^2), label="30 km resolution, no closure")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), fiveday./ (128^2), label="Online closure, 5 day")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), tenday./ (128^2), label="Online closure, 10 day")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), twentyday./ (128^2), label="Online closure, 20 day")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), thirtyday./ (128^2), label="Online closure, 30 day")
    lines!(fig[1,1], LinRange(0, 10*365, 522), multi210 ./ (128^2), label="Online closure, ensemble 2 day")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), multi3 ./ (128^2), label="Online closure, batched 3 day, fewer initial conditions")
    lines!(ax, LinRange(0, 10*365, 522), multi3more10 ./ (128^2), label="Online closure, ensemble 3 day")
    lines!(ax, LinRange(0, 10*365, 522), multi1010 ./ (128^2), label="Online closure, ensemble 10 day")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), multi20 ./ (128^2), label="Online closure, batched 20 day", color=:red3)
    Legend(fig[1, 2], ax)

    # 10 year relative error figure
    fig = Figure(size=(1000, 500), fontsize=15);
    ax = Axis(fig[1,1],
            xlabel="Day",
            ylabel="Energy",
            title="Spatially averaged energy over 3 years"
    )
    lines!(ax, LinRange(0, 3*365, 1096), (abs.(noparam./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    lines!(ax, LinRange(0, 3*365, 1096), (abs.(zb./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{ZB} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    # lines!(fig[1,1], LinRange(0, 3*365, 1096), fiveday./ (128^2), label="Online closure, 5 day")
    # lines!(ax, LinRange(0, 3*365, 1096), tenday./ (128^2), label="Online closure, 10 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentyday ./ (128^2), label="Online closure, 20 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentydaycD ./ (128^2), label="Online closure, 20 day with BD coeff")
    # lines!(ax, LinRange(0, 3*365, 1096), (abs.( thirtyday./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{30} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    # lines!(ax, LinRange(0, 3*365, 1096), multi1 ./ (128^2), label="Online closure, batched 1 day")
    lines!(ax, LinRange(0, 3*365, 1096), (abs.( multi1./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi1} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    # lines!(ax, LinRange(0, 3*365, 1096), (abs.( multi1_more./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi1} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    lines!(ax, LinRange(0, 3*365, 1096), (abs.( multi2./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi2} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    # lines!(ax, LinRange(0, 3*365, 1096), multi2 ./ (128^2), label="Online closure, batched 2 day")
    lines!(ax, LinRange(0, 3*365, 1096), (abs.( multi3_new./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi3} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")
    # lines!(ax, LinRange(0, 3*365, 1096), multi5 ./ (128^2), label="Online closure, batched 5 day", color=:mediumorchid)
    lines!(ax, LinRange(0, 3*365, 1096), (abs.( multi10./ (128^2) .- cghr./ (128^2) ) ./ (cghr ./ (128^2)) ), label=L"|\mathcal{E}_{multi10} - \overline{\mathcal{E}}| / \overline{\mathcal{E}}")#, color=:teal)
    # lines!(ax, LinRange(0, 3*365, 1096), multi20 ./ (128^2), label="Online closure, batched 20 day", color=:red3)
    Legend(fig[1, 2], ax)

    fig = Figure(size=(1200, 500), fontsize=15);
    lines(fig[1,1], LinRange(0, 3*365, 1096), cghr ./ (128^2), label="KE spectrum", 
        axis=(
            xlabel="Day",
            ylabel="Energy",
            title="Spatially averaged energy over 3 years"
        )
    )
    # lines!(fig[1,1], LinRange(0, 3*365, 1096), noparam./ (128^2), label="30 km resolution, no closure")
    lines!(fig[1,1], LinRange(0, 3*365, 1096), zb./ (128^2), label="ZB20")
    lines!(fig[1,1], LinRange(0, 3*365, 1096), hybrid./ (128^2), label="Hybrid")
    lines!(fig[1,1], LinRange(0, 3*365, 1096), kespec./ (128^2), label="KE spectrum")
    lines!(fig[1,1], LinRange(0, 3*365, 1096), kespecpd./ (128^2), label="KE spectrum percent-difference")
    lines!(fig[1,1], LinRange(0, 3*365, 1096), fourier./ (128^2), label="Fourier")
    axislegend(position = (0,1))

    # 3 year, 10 year, 20 year

    fig = Figure(size=(1000, 750), fontsize=15);
    ax = Axis(fig[1,1],
            xlabel="Day",
            ylabel="Energy",
            title="Spatially averaged energy over 3 years"
    )
    lines!(ax, LinRange(0, 3*365, 1096),  hrcg[1:1096] ./ (128^2), label="Filtered, coarse-grained 3.75 km", color=:black)
    lines!(ax, LinRange(0, 3*365, 1096), noparam./ (128^2), label="30 km resolution, no closure",color=:gray)
    lines!(ax, LinRange(0, 3*365, 1096), zb./ (128^2), label="ZB20", color=:red)
    # lines!(fig[1,1], LinRange(0, 3*365, 1096), fiveday./ (128^2), label="Online closure, 5 day")
    # lines!(ax, LinRange(0, 3*365, 1096), tenday./ (128^2), label="Online closure, 10 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentyday ./ (128^2), label="Online closure, 20 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentydaycD ./ (128^2), label="Online closure, 20 day with BD coeff")
    lines!(ax, LinRange(0, 3*365, 1096), thirtyday ./ (128^2), label="Online closure, 30 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi1 ./ (128^2), label="Online closure, batched 1 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi1_more ./ (128^2), label="Online closure, batched 1 day")
    lines!(ax, LinRange(0, 3*365, 1096), multi2 ./ (128^2), label="Online closure, ensemble 2 day", color=:teal)
    lines!(ax, LinRange(0, 3*365, 1096), multi3more ./ (128^2), label="Online closure, ensemble 3 day",color=:darkorchid)
    # lines!(ax, LinRange(0, 3*365, 1096), multi5 ./ (128^2), label="Online closure, batched 5 day", color=:mediumorchid)
    lines!(ax, LinRange(0, 3*365, 1096), multi10 ./ (128^2), label="Online closure, ensemble 10 day")#, color=:teal)
    lines!(ax, LinRange(0, 3*365, 1096), multi20 ./ (128^2), label="Online closure, ensemble 20 day",color=:green)
    Legend(fig[1, 2], ax)

    hrcg_10 = cat(hrcg[1:7:1096], hrcg[1097:end];dims=1)
    ax3 = Axis(fig[2,1],
        xlabel="Day",
        # ylabel="Energy",
        title="Spatially averaged energy over 10 years"
    )
    lines!(ax3, LinRange(0, 10*365, 522),  hrcg_10 ./ (128^2), label="Filtered, coarse-grained 3.75 km",color=:black)
    lines!(ax3, LinRange(0, 10*365, 522), noparam10 ./ (128^2), label="30 km resolution, no closure",color=:gray)
    lines!(ax3, LinRange(0, 10*365, 522),  zb10 ./ (128^2), label="ZB20",color=:red)
    # lines!(fig[1,1], LinRange(0, 10*365, 522), fiveday./ (128^2), label="Online closure, 5 day")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), tenday./ (128^2), label="Online closure, 10 day")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), twentyday./ (128^2), label="Online closure, 20 day")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), thirtyday./ (128^2), label="Online closure, 30 day")
    lines!(ax3, LinRange(0, 10*365, 522), multi210 ./ (128^2), label="Online closure, ensemble 2 day",color=:teal)
    # lines!(fig[1,1], LinRange(0, 10*365, 522), multi3 ./ (128^2), label="Online closure, batched 3 day, fewer initial conditions")
    lines!(ax3, LinRange(0, 10*365, 522), multi3more10 ./ (128^2), label="Online closure, ensemble 3 day",color=:darkorchid)
    lines!(ax3, LinRange(0, 10*365, 522), multi1010 ./ (128^2), label="Online closure, ensemble 10 day", color=:green)
    # lines!(fig[1,1], LinRange(0, 10*365, 522), multi20 ./ (128^2), label="Online closure, batched 20 day", color=:red3)
    Legend(fig[2, 2], ax3)

    ax4 = Axis(fig[3,1],
        xlabel="Day",
        # ylabel="Energy",
        title="Spatially averaged energy over 20 years"
    )
    lines!(ax4, LinRange(0, 20*365, 1043), multi2all ./ (128^2), label="Online closure, ensemble 2 day",color=:teal)
    lines!(ax4, LinRange(0, 20*365, 1043), multi3all ./ (128^2), label="Online closure, ensemble 3 day",color=:darkorchid)
    # lines!(ax3, LinRange(0, 10*365, 522),  cghr10_forplotting ./ (128^2), label="Filtered, coarse-grained 3.75 km",color=:red)

    Legend(fig[3, 2], ax4)

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


    # Diverging results

    fig = Figure(size=(1000, 750), fontsize=15);
    ax = Axis(fig[1,1],
            xlabel="Day",
            ylabel="Energy",
            title="Spatially averaged energy over 3 years"
    )
    lines!(ax, LinRange(0, 3*365, 1096),  hrcg[1:1096] ./ (128^2), label="Filtered, coarse-grained 3.75 km", color=:black)
    lines!(ax, LinRange(0, 3*365, 1096), noparam./ (128^2), label="30 km resolution, no closure",color=:gray)
    lines!(ax, LinRange(0, 3*365, 1096), zb./ (128^2), label="ZB20", color=:red)
    # lines!(ax, LinRange(0, 3*365, 1096), fiveday./ (128^2), label="Online closure, 5 day")
    lines!(ax, LinRange(0, 3*365, 1096), tenday./ (128^2), label="Online closure, 10 day")
    lines!(ax,LinRange(0, 3*365, 1096), twentyday ./ (128^2), label="Online closure, 20 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentydaycD ./ (128^2), label="Online closure, 20 day with BD coeff")
    lines!(ax, LinRange(0, 3*365, 1096), thirtyday ./ (128^2), label="Online closure, 30 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi1 ./ (128^2), label="Online closure, batched 1 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi1_more ./ (128^2), label="Online closure, batched 1 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi2 ./ (128^2), label="Online closure, ensemble 2 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi3more ./ (128^2), label="Online closure, ensemble 3 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi5 ./ (128^2), label="Online closure, batched 5 day", color=:mediumorchid)
    # lines!(ax, LinRange(0, 3*365, 1096), multi10 ./ (128^2), label="Online closure, ensemble 10 day")#, color=:teal)
    lines!(ax, LinRange(0, 3*365, 1096), multi20 ./ (128^2), label="Online closure, ensemble 20 day")
    Legend(fig[1, 2], ax)

    ax3 = Axis(fig[2,1],
        xlabel="Day",
        # ylabel="Energy",
        title="Spatially averaged energy over ten years"
    )
    lines!(ax3, LinRange(0, 10*365, 522),  hrcg_10 ./ (128^2), label="Filtered, coarse-grained 3.75 km",color=:black)
    lines!(ax3, LinRange(0, 10*365, 522), noparam10 ./ (128^2), label="30 km resolution, no closure",color=:gray)
    lines!(ax3, LinRange(0, 10*365, 522),  zb10 ./ (128^2), label="ZB20",color=:red)
    # lines!(fig[1,1], LinRange(0, 10*365, 522), fiveday./ (128^2), label="Online closure, 5 day")
    lines!(ax3, LinRange(0, 10*365, 522), tenday10./ (128^2), label="Online closure, 10 day")
    lines!(ax3, LinRange(0, 10*365, 522), twentyday10./ (128^2), label="Online closure, 20 day")
    lines!(ax3, LinRange(0, 10*365, 522), thirtyday10./ (128^2), label="Online closure, 30 day")
    lines!(ax3, LinRange(0, 10*365, 522), multi2010 ./ (128^2), label="Online closure, ensemble 20 day")

    Legend(fig[2, 2], ax3)

    ga = fig[1, 1] = GridLayout()
    gb = fig[2, 1] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)"], [ga, gb])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end


end

function spectrum_plots()

    # KE spectrum #############################################################

    # to get coarse-grained states
    ker = ImageFiltering.Kernel.gaussian((30e3/3750))
    # imfilter(hru[:,:,j], reflect(ker))

    # for first three years
    totalstates = 1096 # saved every 8 hours (this is when the hr cg and low resolution match up)

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
        xticks=[700, 100, 30, 10, 2],
        title="3-year averaged kinetic energy spectrum"
    )
    ax2 = Axis(fig[2,1],
        xscale=log10,
        yscale=log10,
        xlabel="Wavelength (km)",
        ylabel=L"\text{KE(k) } (m^3/s^2)", xreversed=true,
        xticks=[700, 100, 30, 10, 2],
        title="3-year averaged kinetic energy spectrum"
    )
    lines!(ax, lr_wl[2:end], (up_cghr_avg[2:end] + vp_cghr_avg[2:end])/totalstates, label="Filtered, coarse-grained 3.75km resolution", color=:black)
    lines!(ax, lr_wl[2:end], (up_zb_avg[2:end] + vp_zb_avg[2:end])/totalstates, label="ZB20", color=:red)
    lines!(ax, lr_wl[2:end], (up_noparam_avg[2:end] + vp_noparam_avg[2:end])/totalstates, label="30 km resolution, no closure", color=:gray)

    # lines!(ax, lr_wl[2:end], (up_multi1_avg[2:end] + vp_multi1_avg[2:end])/1096, label="Online NN closure, batched 1 day")#, linestyle=:dash)
    lines!(ax, lr_wl[2:end], (up_multi2_avg[2:end] + vp_multi2_avg[2:end])/totalstates, label="Ensemble 2 day")#, linestyle=:dash)
    lines!(ax, lr_wl[2:end], (up_multi3more_avg[2:end] + vp_multi3more_avg[2:end])/totalstates, label="Ensemble 3 day")#, linestyle=:dashdot)
    lines!(ax, lr_wl[2:end], (up_multi10_avg[2:end] + vp_multi10_avg[2:end])/totalstates, label="Ensemble 10 day")#,linestyle=:dot)
    # lines!(ax, lr_wl[2:end], (up_multi20_avg[2:end] + vp_multi20_avg[2:end])/totalstates, label="Ensemble 20 day")#, color=:red)

    lines!(ax2, lr_wl[2:end], (up_cghr_avg[2:end] + vp_cghr_avg[2:end])/totalstates, label="Filtered, coarse-grained 3.75km resolution", color=:black)
    lines!(ax2, lr_wl[2:end], (up_zb_avg[2:end] + vp_zb_avg[2:end])/totalstates, label="ZB20", color=:red)
    lines!(ax2, lr_wl[2:end], (up_noparam_avg[2:end] + vp_noparam_avg[2:end])/totalstates, label="30 km resolution, no closure", color=:gray)
    lines!(ax2, lr_wl[2:end], (up_5day_avg[2:end] + vp_5day_avg[2:end])/1096, label="5 day")
    lines!(ax2, lr_wl[2:end], (up_10day_avg[2:end] + vp_10day_avg[2:end])/1096, label="10 day")
    lines!(ax2, lr_wl[2:end], (up_20day_avg[2:end] + vp_20day_avg[2:end])/1096, label="20 day")
    lines!(ax2, lr_wl[2:end], (up_30day_avg[2:end] + vp_30day_avg[2:end])/1096, label="30 day")#, linestyle=:dashdot)

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

# using the definition 
#       T(k_x, k_y) = Re( F(u)^* F(S_x) + F(v)^* F(S_y) )
# where ^* is the complex conjugate. It's not clear to me if S should come from the same
# timestep or the prior, because the prior is what when into computing u and v
function computing_ketransfer()
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

    S10 = ShallowWaters.model_setup(Ponline);
    S20 = ShallowWaters.model_setup(Ponline);
    S5 = ShallowWaters.model_setup(Ponline);
    Shrcg = ShallowWaters.model_setup(Ponline);
    SZB = ShallowWaters.model_setup(PZB);
    S30 = ShallowWaters.model_setup(Ponline);

    Smulti3 = ShallowWaters.model_setup(Ponline);
    Smulti5 = ShallowWaters.model_setup(Ponline);

    Smulti10 = ShallowWaters.model_setup(Ponline);
    Smulti20 = ShallowWaters.model_setup(Ponline);

    Smulti1 = ShallowWaters.model_setup(Ponline);
    Smulti1more = ShallowWaters.model_setup(Ponline);
    Smulti2 = ShallowWaters.model_setup(Ponline);

    # original setup
    # uhr = ncread("./dissipation_constant/spinup_files/1024_postspinup_noslip_5years_061824/u.nc", "u");
    # vhr = ncread("./dissipation_constant/spinup_files/1024_postspinup_noslip_5years_061824/v.nc", "v");

    # coarse_grained_hrstates = load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_3years_postspinup_dailysaves.jld2");
    # uhrcg = coarse_grained_hrstates[1];
    # vhrcg = coarse_grained_hrstates[2];

    # Suhr = load_object("./dissipation_constant/results/true_S_SuSv_fromtimederivatives_new.jld2")[1];
    # Svhr = load_object("./dissipation_constant/results/true_S_SuSv_fromtimederivatives_new.jld2")[2];

    # corrected setup
    uhr1 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day1-766saves/u.nc", "u");
    vhr1 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day1-766saves/v.nc", "v");

    uhr2 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day766-end/u.nc", "u");
    vhr2 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day766-end/v.nc", "v");

    coarse_grained_hrstates = load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_imfilter_3years_postspinup_dailysaves_correctedsetup.jld2");
    uhrcg = coarse_grained_hrstates[1];
    vhrcg = coarse_grained_hrstates[2];

    Suhr2 = load_object("./dissipation_constant/results/true_S_SuSv_correctedsetup_fromtimederivatives.jld2")[1]
    Svhr2 = load_object("./dissipation_constant/results/true_S_SuSv_correctedsetup_fromtimederivatives.jld2")[2]

    u10day = ncread("./dissipation_constant/results/128_online_stateweights_10dayoptimization_startfrom5daystate_noeta_oneyear/u.nc", "u");
    v10day = ncread("./dissipation_constant/results/128_online_stateweights_10dayoptimization_startfrom5daystate_noeta_oneyear/v.nc", "v");
    eta10day = ncread("./dissipation_constant/results/128_online_stateweights_10dayoptimization_startfrom5daystate_noeta_oneyear/eta.nc", "eta");

    u5day = ncread("./dissipation_constant/results/128_online_gelu_stateweights_5dayoptimization_startfrom1daystate_3years_dailysaves/u.nc", "u");
    v5day = ncread("./dissipation_constant/results/128_online_gelu_stateweights_5dayoptimization_startfrom1daystate_3years_dailysaves/v.nc", "v");
    eta5day = ncread("./dissipation_constant/results/128_online_gelu_stateweights_5dayoptimization_startfrom1daystate_3years_dailysaves/eta.nc", "eta");

    u20day = ncread("./dissipation_constant/results/result_online_stateweights_20dayoptimization_startfrom10day_fixedcfl_3years_dailysaves/u.nc", "u");
    v20day = ncread("./dissipation_constant/results/result_online_stateweights_20dayoptimization_startfrom10day_fixedcfl_3years_dailysaves/v.nc", "v");
    eta20day = ncread("./dissipation_constant/results/result_online_stateweights_20dayoptimization_startfrom10day_fixedcfl_3years_dailysaves/eta.nc", "eta");

    u30day = ncread("./dissipation_constant/results/result_online_stateweights_30dayoptimization_startfrom20day_fixedcfl_3years_dailysaves/u.nc", "u");
    v30day = ncread("./dissipation_constant/results/result_online_stateweights_30dayoptimization_startfrom20day_fixedcfl_3years_dailysaves/v.nc", "v");
    eta30day = ncread("./dissipation_constant/results/result_online_stateweights_30dayoptimization_startfrom20day_fixedcfl_3years_dailysaves/eta.nc", "eta");

    uzb_ = ncread("./dissipation_constant/results/128_ZBparam_postspinup_cginitcond_3years_dailysaves/u.nc", "u");
    vzb_ = ncread("./dissipation_constant/results/128_ZBparam_postspinup_cginitcond_3years_dailysaves/v.nc", "v");
    etazb_ = ncread("./dissipation_constant/results/128_ZBparam_postspinup_cginitcond_3years_dailysaves/eta.nc", "eta");

    umulti3 = ncread("./dissipation_constant/results/128_online_multistateweights_3dayoptimization_3-25-30-40-50-60-80initdays_startfrom20daystate_3years_dailysaves/u.nc", "u");
    vmulti3 = ncread("./dissipation_constant/results/128_online_multistateweights_3dayoptimization_3-25-30-40-50-60-80initdays_startfrom20daystate_3years_dailysaves/v.nc", "v");
    etamulti3 = ncread("./dissipation_constant/results/128_online_multistateweights_3dayoptimization_3-25-30-40-50-60-80initdays_startfrom20daystate_3years_dailysaves/eta.nc", "eta");

    umulti5 = ncread("./dissipation_constant/results/128_multistateweights_5dayoptimization_3-30-50-80initdays_fixedcfl_constdiss_3years_dailysaves/u.nc", "u")
    vmulti5 = ncread("./dissipation_constant/results/128_multistateweights_5dayoptimization_3-30-50-80initdays_fixedcfl_constdiss_3years_dailysaves/v.nc", "v")
    etamulti5 = ncread("./dissipation_constant/results/128_multistateweights_5dayoptimization_3-30-50-80initdays_fixedcfl_constdiss_3years_dailysaves/eta.nc", "eta")

    umulti10 = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_3years_dailysaves/u.nc", "u");
    vmulti10 = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_3years_dailysaves/v.nc", "v");
    etamulti10 = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_3years_dailysaves/eta.nc", "eta");

    umulti20 = ncread("./dissipation_constant/results/result_online_multistateweights_20dayoptimization_5-25-45-65initdays_startfrom20daystate_3years_dailysaves/u.nc", "u");
    vmulti20 = ncread("./dissipation_constant/results/result_online_multistateweights_20dayoptimization_5-25-45-65initdays_startfrom20daystate_3years_dailysaves/v.nc", "v");
    etamulti20 = ncread("./dissipation_constant/results/result_online_multistateweights_20dayoptimization_5-25-45-65initdays_startfrom20daystate_3years_dailysaves/eta.nc", "eta");

    umulti1 = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_35initdays_startfrommulti3_3years_dailysaves/u.nc", "u");
    vmulti1 = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_35initdays_startfrommulti3_3years_dailysaves/v.nc", "v");
    etamulti1 = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_35initdays_startfrommulti3_3years_dailysaves/eta.nc", "eta");

    umulti1more = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_1:2:89initdays_startfrommulti3_fewerinitconds_3years_dailysaves/u.nc", "u");
    vmulti1more = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_1:2:89initdays_startfrommulti3_fewerinitconds_3years_dailysaves/v.nc", "v");
    etamulti1more = ncread("./dissipation_constant/results/result_online_multistateweights_1dayoptimization_1:2:89initdays_startfrommulti3_fewerinitconds_3years_dailysaves/eta.nc", "eta");

    umulti2 = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/u.nc", "u");
    vmulti2 = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/v.nc", "v");
    etamulti2 = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/eta.nc", "eta");

    lr_freq = 1/30 .* freq(periodogram(umulti1[:,:,10]; radialavg=true, radialsum=false));
    nfft = nextfastfft(size(uhrcg[:,:,1]))

    totalu_hrcg = zeros(65)
    totalv_hrcg = zeros(65)

    totalu_hrcg2 = zeros(65)
    totalv_hrcg2 = zeros(65)

    totalu_5 = zeros(65)
    totalv_5 = zeros(65)

    totalu_10 = zeros(65)
    totalv_10 = zeros(65)

    totalu_20 = zeros(65)
    totalv_20 = zeros(65)
    
    totalu_30 = zeros(65)
    totalv_30 = zeros(65)

    totalu_multi1 = zeros(65)
    totalv_multi1 = zeros(65)

    totalu_multi1more = zeros(65)
    totalv_multi1more = zeros(65)

    totalu_multi2 = zeros(65)
    totalv_multi2 = zeros(65)

    totalu_multi3 = zeros(65)
    totalv_multi3 = zeros(65)

    totalu_multi5 = zeros(65)
    totalv_multi5 = zeros(65)

    totalu_multi10 = zeros(65)
    totalv_multi10 = zeros(65)

    totalu_multi20 = zeros(65)
    totalv_multi20 = zeros(65)

    totalu_ZB = zeros(65)
    totalv_ZB = zeros(65)

    totalu_3 = zeros(65)
    totalv_3 = zeros(65)

    totalstates = 1096
    for t = 1:totalstates

        # coarse-grained high-resolution computation, corrected setup
        if t in 1:766
            Suhr, Svhr = compute_hrS(uhr1[:,:,t], vhr1[:,:, t])
        else
            Suhr, Svhr = compute_hrS(uhr2[:,:,t-766+1], vhr2[:,:, t-766+1])
        end

        outu_hrcg, inputu_hrcg, inputSu_hrcg = paddingu(uhrcg[:, :, t], Suhr, nfft[1])
        fft2pow2radial!(outu_hrcg, rfft(inputu_hrcg), rfft(inputSu_hrcg), nfft...)
        outv_hrcg, inputv_hrcg, inputSv_hrcg = paddingv(vhrcg[:, :, t], Svhr, nfft[1])
        fft2pow2radial!(outv_hrcg, rfft(inputv_hrcg), rfft(inputSv_hrcg), nfft...)

        totalu_hrcg += outu_hrcg
        totalv_hrcg += outv_hrcg

        # outu_hrcg2, inputu_hrcg2, inputSu_hrcg2 = paddingu(uhrcg2[:, :, t], Suhr2[:,:,t], nfft[1])
        # fft2pow2radial!(outu_hrcg2, rfft(inputu_hrcg2), rfft(inputSu_hrcg2), nfft...)
        # outv_hrcg2, inputv_hrcg2, inputSv_hrcg2 = paddingv(vhrcg2[:, :, t], Svhr2[:,:,t], nfft[1])
        # fft2pow2radial!(outv_hrcg2, rfft(inputv_hrcg2), rfft(inputSv_hrcg2), nfft...)

        # totalu_hrcg2 += outu_hrcg2
        # totalv_hrcg2 += outv_hrcg2

        # ZB20 ############################

        # ShallowWaters.ZB_momentum(uzb, vzb, SZB, SZB.Diag);
        # uzb, vzb, etazb = ShallowWaters.add_halo(Float64.(uzb_[:,:,t]), Float64.(vzb_[:,:,t]), Float64.(etazb_[:,:,t]), zeros(128,128), SZB);

        # outu_ZB, inputu_ZB, inputSu_ZB = paddingu(uzb_[:, :, t], SZB.Diag.ZBVars.S_u, nfft[1])
        # fft2pow2radial!(outu_ZB, rfft(inputu_ZB), rfft(inputSu_ZB), nfft...)
        # outv_ZB, inputv_ZB, inputSv_ZB = paddingv(vzb_[:, :, t], SZB.Diag.ZBVars.S_v, nfft[1])
        # fft2pow2radial!(outv_ZB, rfft(inputv_ZB), rfft(inputSv_ZB), nfft...)

        # totalu_ZB += outu_ZB
        # totalv_ZB += outv_ZB

        # 5 day optimization ################

        # u5, v5, eta5 = ShallowWaters.add_halo(Float64.(u5day[:,:,t]), Float64.(v5day[:,:,t]), Float64.(eta5day[:,:,t]), zeros(128,128), S5);
        # ShallowWaters.CNN_momentum(u5, v5, S5);
        # outu_5, inputu_5, inputSu_5 = paddingu(u5day[:, :, t], S5.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_5, rfft(inputu_5), rfft(inputSu_5), nfft...)
        # outv_5, inputv_5, inputSv_5 = paddingv(v5day[:, :, t], S5.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_5, rfft(inputv_5), rfft(inputSv_5), nfft...)

        # totalu_5 += outu_5
        # totalv_5 += outv_5

        # 20 day optimization ##############

        # u20, v20, _ = ShallowWaters.add_halo(Float64.(u20day[:,:,t]), Float64.(v20day[:,:,t]), Float64.(eta20day[:,:,t]), zeros(128,128), S20);
        # ShallowWaters.CNN_momentum(u20, v20, S20)
        # outu_20, inputu_20, inputSu_20 = paddingu(u20day[:, :, t], S20.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_20, rfft(inputu_20), rfft(inputSu_20), nfft...)
        # outv_20, inputv_20, inputSv_20 = paddingv(v20day[:, :, t], S20.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_20, rfft(inputv_20), rfft(inputSv_20), nfft...)

        # totalu_20 += outu_20
        # totalv_20 += outv_20

        # 30 day optimization ###############

        # u30, v30, _ = ShallowWaters.add_halo(Float64.(u30day[:,:,t]), Float64.(v30day[:,:,t]), Float64.(eta30day[:,:,t]), zeros(128,128), S30);
        # ShallowWaters.CNN_momentum(u30, v30, S30)
        # outu_30, inputu_30, inputSu_30 = paddingu(u30day[:, :, t], S30.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_30, rfft(inputu_30), rfft(inputSu_30), nfft...)
        # outv_30, inputv_30, inputSv_30 = paddingv(v30day[:, :, t], S30.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_30, rfft(inputv_30), rfft(inputSv_30), nfft...)

        # totalu_30 += outu_30
        # totalv_30 += outv_30

        # batched 1 day wtih fewer initial conditions (35) ################

        # umulti1_, vmulti1_, _ = ShallowWaters.add_halo(Float64.(umulti1[:,:,t]), Float64.(vmulti1[:,:,t]), Float64.(etamulti1[:,:,t]), zeros(128,128), Smulti1);  
        # ShallowWaters.CNN_momentum(umulti1_, vmulti1_, Smulti1)
        # outu_multi1, inputu_multi1, inputSu_multi1 = paddingu(umulti1[:, :, t], Smulti1.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_multi1, rfft(inputu_multi1), rfft(inputSu_multi1), nfft...)
        # outv_multi1, inputv_multi1, inputSv_multi1 = paddingv(vmulti1[:, :, t], Smulti1.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_multi1, rfft(inputv_multi1), rfft(inputSv_multi1), nfft...)

        # totalu_multi1 += outu_multi1
        # totalv_multi1 += outv_multi1

        # batched 1 day with a ton of initial conditions #############################

        # umulti1more_, vmulti1more_, _ = ShallowWaters.add_halo(Float64.(umulti1more[:,:,t]), Float64.(vmulti1more[:,:,t]), Float64.(etamulti1more[:,:,t]), zeros(128,128), Smulti1more);
        # ShallowWaters.CNN_momentum(umulti1more_, vmulti1more_, Smulti1more)
        # outu_multi1more, inputu_multi1more, inputSu_multi1more = paddingu(umulti1more[:, :, t], Smulti1more.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_multi1more, rfft(inputu_multi1more), rfft(inputSu_multi1more), nfft...)
        # outv_multi1more, inputv_multi1more, inputSv_multi1more = paddingv(vmulti1more[:, :, t], Smulti1more.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_multi1more, rfft(inputv_multi1more), rfft(inputSv_multi1more), nfft...)

        # totalu_multi1more += outu_multi1more
        # totalv_multi1more += outv_multi1more

        # batched 2 day ##################################

        # umulti2_, vmulti2_, _ = ShallowWaters.add_halo(Float64.(umulti2[:,:,t]), Float64.(vmulti2[:,:,t]), Float64.(etamulti2[:,:,t]), zeros(128,128), Smulti2);
        # ShallowWaters.CNN_momentum(umulti2_, vmulti2_, Smulti2)
        # outu_multi2, inputu_multi2, inputSu_multi2 = paddingu(umulti2[:, :, t], Smulti2.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_multi2, rfft(inputu_multi2), rfft(inputSu_multi2), nfft...)
        # outv_multi2, inputv_multi2, inputSv_multi2 = paddingv(vmulti2[:, :, t], Smulti2.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_multi2, rfft(inputv_multi2), rfft(inputSv_multi2), nfft...)

        # totalu_multi2 += outu_multi2
        # totalv_multi2 += outv_multi2

        # batched 3 day ###########################

        # umulti3_, vmulti3_, _ = ShallowWaters.add_halo(Float64.(umulti3[:,:,t]), Float64.(vmulti3[:,:,t]), Float64.(etamulti3[:,:,t]), zeros(128,128), Smulti3);
        # ShallowWaters.CNN_momentum(umulti3_, vmulti3_, Smulti3)
        # outu_multi3, inputu_multi3, inputSu_multi3 = paddingu(umulti3[:, :, t], Smulti3.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_multi3, rfft(inputu_multi3), rfft(inputSu_multi3), nfft...)
        # outv_multi3, inputv_multi3, inputSv_multi3 = paddingv(vmulti3[:, :, t], Smulti3.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_multi3, rfft(inputv_multi3), rfft(inputSv_multi3), nfft...)

        # totalu_multi3 += outu_multi3
        # totalv_multi3 += outv_multi3

        # batched 5 day ######################

        # umulti5_, vmulti5_, _ = ShallowWaters.add_halo(Float64.(umulti5[:,:,t]), Float64.(vmulti5[:,:,t]), Float64.(etamulti5[:,:,t]), zeros(128,128), Smulti5);
        # ShallowWaters.CNN_momentum(umulti5_, vmulti5_, Smulti5)
        # outu_multi5, inputu_multi5, inputSu_multi5 = paddingu(umulti5[:, :, t], Smulti5.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_multi5, rfft(inputu_multi5), rfft(inputSu_multi5), nfft...)
        # outv_multi5, inputv_multi5, inputSv_multi5 = paddingv(vmulti5[:, :, t], Smulti5.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_multi5, rfft(inputv_multi5), rfft(inputSv_multi5), nfft...)

        # totalu_multi5 += outu_multi5
        # totalv_multi5 += outv_multi5

        # batched 10 day #####################

        # umulti10_, vmulti10_, _ = ShallowWaters.add_halo(Float64.(umulti10[:,:,t]), Float64.(vmulti10[:,:,t]), Float64.(etamulti10[:,:,t]), zeros(128,128), Smulti10);
        # ShallowWaters.CNN_momentum(umulti10_, vmulti10_, Smulti10)
        # outu_multi10, inputu_multi10, inputSu_multi10 = paddingu(umulti10[:, :, t], Smulti10.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_multi10, rfft(inputu_multi10), rfft(inputSu_multi10), nfft...)
        # outv_multi10, inputv_multi10, inputSv_multi10 = paddingv(vmulti10[:, :, t], Smulti10.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_multi10, rfft(inputv_multi10), rfft(inputSv_multi10), nfft...)

        # totalu_multi10 += outu_multi10
        # totalv_multi10 += outv_multi10

        # batched 20 day #####################

        # umulti20_, vmulti20_, _ = ShallowWaters.add_halo(Float64.(umulti20[:,:,t]), Float64.(vmulti20[:,:,t]), Float64.(etamulti20[:,:,t]), zeros(128,128), Smulti20);
        # ShallowWaters.CNN_momentum(umulti20_, vmulti20_, Smulti20)
        # outu_multi20, inputu_multi20, inputSu_multi20 = paddingu(umulti20[:, :, t], Smulti20.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_multi20, rfft(inputu_multi20), rfft(inputSu_multi20), nfft...)
        # outv_multi20, inputv_multi20, inputSv_multi20 = paddingv(vmulti20[:, :, t], Smulti20.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_multi20, rfft(inputv_multi20), rfft(inputSv_multi20), nfft...)

        # totalu_multi20 += outu_multi20
        # totalv_multi20 += outv_multi20

    end
end

function ketransfer_plots()

    totalu_hrcg = load_object("./dissipation_constant/spinup_files/hrcg_ketransfer_fromapprox_correctedsetup_uv_3years.jld2")[1];
    totalv_hrcg = load_object("./dissipation_constant/spinup_files/hrcg_ketransfer_fromapprox_correctedsetup_uv_3years.jld2")[2];

    totalu_ZB = load_object("./dissipation_constant/results/zb20_ketransfer_uv_3years.jld2")[1];
    totalv_ZB = load_object("./dissipation_constant/results/zb20_ketransfer_uv_3years.jld2")[2];

    totalu_30 = load_object("./dissipation_constant/results/30day_stateloss_ketransfer_uv_3years.jld2")[1];
    totalv_30 = load_object("./dissipation_constant/results/30day_stateloss_ketransfer_uv_3years.jld2")[2];

    totalu_20 = load_object("./dissipation_constant/results/20day_stateloss_ketransfer_uv_3years.jld2")[1];
    totalv_20 = load_object("./dissipation_constant/results/20day_stateloss_ketransfer_uv_3years.jld2")[2];
    
    totalu_5 = load_object("./dissipation_constant/results/5day_ketransfer_uv_3years.jld2")[1];
    totalv_5 = load_object("./dissipation_constant/results/5day_ketransfer_uv_3years.jld2")[2];

    totalu_multi3_old = load_object("./dissipation_constant/results/multi3day_stateloss_3-25-30-40-50-60-80-initdays_ketransfer_uv_3years.jld2")[1];
    totalv_multi3_old = load_object("./dissipation_constant/results/multi3day_stateloss_3-25-30-40-50-60-80-initdays_ketransfer_uv_3years.jld2")[2];

    totalu_multi3_new = load_object("./dissipation_constant/results/multi3day_stateloss_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_ketransfer_uv_3years.jld2")[1]
    totalv_multi3_new = load_object("./dissipation_constant/results/multi3day_stateloss_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_ketransfer_uv_3years.jld2")[2]

    totalu_multi5 = load_object("./dissipation_constant/results/multi5day_stateloss_3-30-50-80-initdays_ketransfer_uv_3years.jld2")[1];
    totalv_multi5 = load_object("./dissipation_constant/results/multi5day_stateloss_3-30-50-80-initdays_ketransfer_uv_3years.jld2")[2];

    totalu_multi10 = load_object("./dissipation_constant/results/multi10day_stateloss_5-20-35-50-65-75initdays_ketransfer_uv_3years.jld2")[1];
    totalv_multi10 = load_object("./dissipation_constant/results/multi10day_stateloss_5-20-35-50-65-75initdays_ketransfer_uv_3years.jld2")[2];

    totalu_multi20 = load_object("./dissipation_constant/results/multi20day_stateloss_5-25-45-65initdays_ketransfer_uv_3years.jld2")[1];
    totalv_multi20 = load_object("./dissipation_constant/results/multi20day_stateloss_5-25-45-65initdays_ketransfer_uv_3years.jld2")[2];

    lr_freq = 1/30 .* freq(periodogram(u20s[:,:,10]; radialavg=true, radialsum=false));
    nfft = nextfastfft(size(uhrcg[:,:,1]))

    fig = Figure(size=(800, 300), fontsize=15);
    ax = Axis(fig[1,1],
        xscale = log10,
        xlabel="Wavenumber (1/km)",
        ylabel="KE(k)",
        title="Kinetic energy transfer"
    )
    lines!(ax, lr_freq.*(totalu_hrcg + totalv_hrcg) / 1096, label="Subgrid forcing")
    lines!(ax, lr_freq.*(totalu_ZB + totalv_ZB) / 1098, label="ZB20")
    lines!(ax, lr_freq.*(totalu_multi1 + totalv_multi1) / 1096, label="Batched 1 day optimization")
    # lines!(ax, lr_freq.*(totalu_multi1more + totalv_multi1more) / 1096, label="Batched 1 day optimization, more init conditions")
    lines!(ax, lr_freq.*(totalu_multi2 + totalv_multi2) / 1096, label="Batched 2 day optimization")

    lines!(ax, lr_freq.*(totalu_30 + totalv_30) / 1096, label="30 day optimization")
    lines!(ax, lr_freq.*(totalu_20 + totalv_20) / 1096, label="20 day optimization")
    # lines!(ax, lr_freq.*(totalu_multi3_old + totalv_multi3_old) / 1096, label="Multi 3 day state optimization old")
    lines!(ax, lr_freq.*(totalu_multi3_new + totalv_multi3_new) / 1096, label="Multi 3 day state optimization")
    lines!(ax, lr_freq.*(totalu_multi5 + totalv_multi5) / 1096, label="Multi 5 day state optimization")
    lines!(ax, lr_freq.*(totalu_multi10 + totalv_multi10) / 1096, label="Multi 10 day state optimization")
    lines!(ax, lr_freq.*(totalu_multi20 + totalv_multi20) / 1096, label="Multi 20 day state optimization")
    Legend(fig[1,2], ax)

end

function parameterization_S_plots()

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
        zb_forcing_dissipation=false,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=2,
        initial_cond="rest"
    );
    S = ShallowWaters.model_setup(P);

    Scghr = deepcopy(S);
    S5 = deepcopy(S);
    S10 = deepcopy(S);
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
    S20 = deepcopy(S);
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
    S30 = deepcopy(S);
    current = 1
    for m in (S30.Diag.CNNVars.model_Su, S30.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_state_30dayoptimzation_startfrom30day6iterations_constantdissipation_15iterations_21totaliterations_8hourdata_200maxhistory_fixedcfl.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    Smulti2 = deepcopy(S);
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

    Smulti3 = deepcopy(S);
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

    Smulti10 = deepcopy(S);
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

    Smulti20 = deepcopy(S);
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

    Soffline = deepcopy(S);
    current = 1
    for m in (Soffline.Diag.CNNVars.model_Su, Soffline.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/result_offline_150iterations_geluactivation_111925.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    Szb = deepcopy(S);
    Scghr = deepcopy(S);

    # plot within the first three years
    t = 1096

    # coarse-grained high-resolution computation, corrected setup
    if t in 1:766
        Suhr, Svhr = compute_hrS(uhr1[:,:,t], vhr1[:,:, t])
    else
        Suhr, Svhr = compute_hrS(uhr2[:,:,t-766+1], vhr2[:,:, t-766+1])
    end

    uzb_, vzb_, etazb_ = ShallowWaters.add_halo(Float64.(uzb[:,:,t]), Float64.(vzb[:,:,t]), Float64.(etazb[:,:,t]), zeros(128,128), S);
    u10s_, v10s_, eta10s_ = ShallowWaters.add_halo(Float64.(u10s[:,:,t]), Float64.(v10s[:,:,t]), Float64.(eta10s[:,:,t]), zeros(128,128), S);
    u20s_, v20s_, eta20s_ = ShallowWaters.add_halo(Float64.(u20s[:,:,t]), Float64.(v20s[:,:,t]), Float64.(eta20s[:,:,t]), zeros(128,128), S);
    u30s_, v30s_, eta30s_ = ShallowWaters.add_halo(Float64.(u30s[:,:,t]), Float64.(v30s[:,:,t]), Float64.(eta30s[:,:,t]), zeros(128,128), S);
    # umulti1_, vmulti1_, etamulti1_ = ShallowWaters.add_halo(Float64.(umulti1[:,:,t]), Float64.(vmulti1[:,:,t]), Float64.(etamulti1[:,:,t]), zeros(128,128), S);
    # umulti1more_, vmulti1more_, etamulti1more_ = ShallowWaters.add_halo(Float64.(umulti1more[:,:,t]), Float64.(vmulti1more[:,:,t]), Float64.(etamulti1more[:,:,t]), zeros(128,128), S);
    umulti2_, vmulti2_, etamulti2_ = ShallowWaters.add_halo(Float64.(umulti2[:,:,t]), Float64.(vmulti2[:,:,t]), Float64.(etamulti2[:,:,t]), zeros(128,128), S);
    # umulti3_, vmulti3_, etamulti3_ = ShallowWaters.add_halo(Float64.(umulti3[:,:,t]), Float64.(vmulti3[:,:,t]), Float64.(etamulti3[:,:,t]), zeros(128,128), S);
    umulti3more_, vmulti3more_, etamulti3more_ = ShallowWaters.add_halo(Float64.(umulti3more[:,:,t]), Float64.(vmulti3more[:,:,t]), Float64.(etamulti3more[:,:,t]), zeros(128,128), S);
    umulti10_, vmulti10_, etamulti10_ = ShallowWaters.add_halo(Float64.(umulti10[:,:,t]), Float64.(vmulti10[:,:,t]), Float64.(etamulti10[:,:,t]), zeros(128,128), S);
    umulti20_, vmulti20_, etamulti20_ = ShallowWaters.add_halo(Float64.(umulti20[:,:,t]), Float64.(vmulti20[:,:,t]), Float64.(etamulti20[:,:,t]), zeros(128,128), S);

    ShallowWaters.ZB_momentum(Float64.(uzb_), Float64.(vzb_), Szb, Szb.Diag);
    ShallowWaters.CNN_momentum(u20s_, v20s_, S20);
    ShallowWaters.CNN_momentum(u30s_, v30s_, S30);
    ShallowWaters.CNN_momentum(u10s_, v10s_, S10);

    ShallowWaters.CNN_momentum(umulti2_, vmulti2_, Smulti2);
    ShallowWaters.CNN_momentum(umulti3more_, vmulti3more_, Smulti3);
    ShallowWaters.CNN_momentum(umulti10_, vmulti10_, Smulti10);
    ShallowWaters.CNN_momentum(umulti20_, vmulti20_, Smulti20);

    fig = Figure(size=(1040, 520), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Suhr,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\overline{S_u}(3 \; \text{years}, x, y)"),
    colorrange=(-maximum(abs.(Suhr)),maximum(abs.(Suhr)))
    );
    Colorbar(fig[1,2], hm1, label="1/s")

    ax1, hm1 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Szb.Diag.ZBVars.S_u,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"S_u \text{, ZB20}"),
    colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u)),maximum(abs.(Szb.Diag.ZBVars.S_u)))
    );
    Colorbar(fig[1,4], hm1, label="1/s")

    ax2, hm2 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti2.Diag.CNNVars.S_u,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"S_u \text{, ensemble 2 day}"),
    colorrange=(-maximum(abs.(Suhr)),maximum(abs.(Suhr)))
    );
    Colorbar(fig[1,6], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti3.Diag.CNNVars.S_u,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"S_u \text{, ensemble 3 day}"),
    colorrange=(-maximum(abs.(Suhr)),maximum(abs.(Suhr)))
    );
    Colorbar(fig[2,2], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    S20.Diag.CNNVars.S_u,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"S_u \text{, 10 day}"),
    colorrange=(-maximum(abs.(Suhr)),maximum(abs.(Suhr)))
    );
    Colorbar(fig[2,4], hm1, label="1/s")

    ax3, hm3 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    S30.Diag.CNNVars.S_u,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"S_u \text{, 30 day}"),
    colorrange=(-maximum(abs.(Suhr)),maximum(abs.(Suhr)))
    );
    Colorbar(fig[2,6], hm1, label="1/s")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # comparing based on cghr snapshot results

    t = 1096
    uhrcg_, vhrcg_, etahrcg_ = ShallowWaters.add_halo(Float64.(uhrcg1[:,:,t]), Float64.(vhrcg1[:,:,t]), Float64.(etahrcg1[:,:,t]), zeros(128,128), S);

    ShallowWaters.ZB_momentum(uhrcg_, vhrcg_, Szb, Szb.Diag);
    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, S20);
    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, S30);
    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, S10);

    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, Smulti2);
    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, Smulti3);
    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, Smulti10);
    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, Smulti20);
    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, Soffline);

    # S_u
    fig = Figure(size=(1040, 520), fontsize=15);

    Label(
        fig[0, 3],
        L"S_u(3 \text{ years}, x, y)",
        fontsize = 20,
        tellwidth = false
    )

    ax2, hm2 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Szb.Diag.ZBVars.S_u,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="ZB20"),
    colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u)),maximum(abs.(Szb.Diag.ZBVars.S_u)))
    );
    Colorbar(fig[1,2], hm1, label=L"m/s^2")

    ax1, hm1 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Soffline.Diag.CNNVars.S_u,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Offline-learned NN"),
    colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u)),maximum(abs.(Szb.Diag.ZBVars.S_u)))
    );
    Colorbar(fig[1,4], hm1, label=L"m/s^2")

    ax2, hm2 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti2.Diag.CNNVars.S_u,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 2 day"),
    colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u)),maximum(abs.(Szb.Diag.ZBVars.S_u)))
    );
    Colorbar(fig[1,6], hm1, label=L"m/s^2")

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti3.Diag.CNNVars.S_u,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 3 day"),
    colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u)),maximum(abs.(Szb.Diag.ZBVars.S_u)))
    );
    Colorbar(fig[2,2], hm1, label=L"m/s^2")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti10.Diag.CNNVars.S_u,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 10 day"),
    colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u)),maximum(abs.(Szb.Diag.ZBVars.S_u)))
    );
    Colorbar(fig[2,4], hm1, label=L"m/s^2")

    ax3, hm3 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    S30.Diag.CNNVars.S_u,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="30 day"),
    colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u)),maximum(abs.(Szb.Diag.ZBVars.S_u)))
    );
    Colorbar(fig[2,6], hm1, label=L"m/s^2")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # S_v
   fig = Figure(size=(1040, 520), fontsize=15);

    Label(
        fig[0, 3],
        L"S_v(3 \text{ years}, x, y)",
        fontsize = 20,
        tellwidth = false
    )

    ax2, hm2 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Szb.Diag.ZBVars.S_v,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="ZB20"),
    colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_v)),maximum(abs.(Szb.Diag.ZBVars.S_v)))
    );
    Colorbar(fig[1,2], hm1, label=L"m/s^2")

    ax1, hm1 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Soffline.Diag.CNNVars.S_v,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Offline-learned NN"),
    colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_v)),maximum(abs.(Szb.Diag.ZBVars.S_v)))
    );
    Colorbar(fig[1,4], hm1, label=L"m/s^2")

    ax2, hm2 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti2.Diag.CNNVars.S_v,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 2 day"),
    colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_v)),maximum(abs.(Szb.Diag.ZBVars.S_v)))
    );
    Colorbar(fig[1,6], hm1, label=L"m/s^2")

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti3.Diag.CNNVars.S_v,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 3 day"),
    colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_v)),maximum(abs.(Szb.Diag.ZBVars.S_v)))
    );
    Colorbar(fig[2,2], hm1, label=L"m/s^2")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti10.Diag.CNNVars.S_v,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 10 day"),
    colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_v)),maximum(abs.(Szb.Diag.ZBVars.S_v)))
    );
    Colorbar(fig[2,4], hm1, label=L"m/s^2")

    ax3, hm3 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    S30.Diag.CNNVars.S_v,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="30 day"),
    colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_v)),maximum(abs.(Szb.Diag.ZBVars.S_v)))
    );
    Colorbar(fig[2,6], hm1, label=L"m/s^2")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end


end

function appendix_plots()

     ## Appendix figures

    # comparing the offline results relu versus gelu
    t = [4, 7, 10]
    fig = Figure(size=(1040, 520), fontsize=15);

    ax2, hm2 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    vofflinegelu[:,:,t[1]],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"v_{\text{offline}}(1 \text{ day}, x, y), \; \text{GELU}"),
    colorrange=(-maximum(abs.(vhrcg[:,:,t[1]])),maximum(abs.(vhrcg[:,:,t[1]])))
    );
    Colorbar(fig[1,2], hm2, label="m/s")

    ax3, hm3 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    vofflinegelu[:,:,t[2]],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"v_{\text{offline}}(2 \text{ days}, x, y), \text{ GELU}"),
    colorrange=(-maximum(abs.(vhrcg[:,:,t[2]])),maximum(abs.(vhrcg[:,:,t[2]])))
    );
    Colorbar(fig[1,4], hm3, label="m/s")

    ax4, hm4 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    vofflinegelu[:,:,t[3]],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"v_{\text{offline}}(3 \text{ days}, x, y), \; \text{GELU}"),
    colorrange=(-maximum(abs.(vhrcg[:,:,t[3]])),maximum(abs.(vhrcg[:,:,t[3]])))
    );
    Colorbar(fig[1,6], hm4, label="m/s")

    ax2, hm2 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    vofflinerelu[:,:,t[1]],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"v_{\text{offline}}(1 \text{ day}, x, y), \; \text{ReLU}"),
    colorrange=(-maximum(abs.(vhrcg[:,:,t[1]])),maximum(abs.(vhrcg[:,:,t[1]])))
    );
    Colorbar(fig[2,2], hm2, label="m/s")

    ax3, hm3 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    vofflinerelu[:,:,t[2]],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"v_{\text{offline}}(2 \text{ days}, x, y), \; \text{ReLU}"),
    colorrange=(-maximum(abs.(vhrcg[:,:,t[2]])),maximum(abs.(vhrcg[:,:,t[2]])))
    );
    Colorbar(fig[2,4], hm3, label="m/s")

    ax4, hm4 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    vofflinerelu[:,:,t[3]],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"v_{\text{offline}}(3 \text{ days}, x, y), \; \text{ReLU}"),
    colorrange=(-maximum(abs.(vhrcg[:,:,t[3]])),maximum(abs.(vhrcg[:,:,t[3]])))
    )
    Colorbar(fig[2,6], hm4, label="m/s")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end


    # comparing the online results relu versus gelu
    t = 91
    fig = Figure(size=(700, 525), fontsize=15);

    ax2, hm2 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    u1daystategelu[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u_1(30 \text{ days}, x, y), \; \text{GELU}"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[1,2], hm2, label="m/s")

    ax3, hm3 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    u5daystategelu[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u_{1 + 5}(30 \text{ days}, x, y), \text{ GELU}"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[1,4], hm3, label="m/s")

    # ax4, hm4 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # ukespecpd1day1daystartgelu[:,:,t],
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title=L"u_{\text{1 + KE pd}}(30 \text{ days}, x, y), \; \text{GELU}"),
    # colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    # );
    # Colorbar(fig[1,6], hm4, label="m/s")

    ax2, hm2 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    u1daystaterelu[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u_1(30 \text{ days}, x, y), \; \text{ReLU}"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[2,2], hm2, label="m/s")

    ax3, hm3 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    u5daystaterelu[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u_{1 + 5}(30 \text{ days}, x, y), \; \text{ReLU}"),
    colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    );
    Colorbar(fig[2,4], hm3, label="m/s")

    # ax4, hm4 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # ukespecpd1dayrelu[:,:,t],
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title=L"u_{\text{1 + KE pd}}(30 \text{ days}, x, y), \; \text{ReLU}"),
    # colorrange=(-maximum(abs.(uhrcg[:,:,t])),maximum(abs.(uhrcg[:,:,t])))
    # );
    # Colorbar(fig[2,6], hm4, label="m/s")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    # gc = fig[1, 5] = GridLayout()
    gc = fig[2, 1] = GridLayout()
    gd = fig[2, 3] = GridLayout()
    # gf = fig[2, 5] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", ], [ga, gb, gc, gd])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # looking at the loss functions that didn't work

    # some prognostic fields
    fig = Figure(size=(700, 525), fontsize=15);

    t = 90
    ax2, hm2 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etakespec[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\eta_{\text{KE spectrum}}(90 \text{ days}, x, y)"),
    colorrange=(-maximum(abs.(etahrcg[:,:,t])),maximum(abs.(etahrcg[:,:,t])))
    );
    Colorbar(fig[1,2], hm2, label="m/s")

    ax3, hm3 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etakespecpd[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\eta_{\text{KE spectrum pd}}(90 \text{ days}, x, y)"),
    colorrange=(-maximum(abs.(etahrcg[:,:,t])),maximum(abs.(etahrcg[:,:,t])))
    );
    Colorbar(fig[1,4], hm3, label="m/s")

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etafourier[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\eta_{\text{Fourier}}(90 \text{ days}, x, y)"),
    colorrange=(-maximum(abs.(etahrcg[:,:,t])),maximum(abs.(etahrcg[:,:,t])))
    );
    Colorbar(fig[2,2], hm4, label="m/s")

    ax2, hm2 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etahybrid[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\eta_{\text{hybrid}}(90 \text{ days}, x, y)"),
    colorrange=(-maximum(abs.(etahrcg[:,:,t])),maximum(abs.(etahrcg[:,:,t])))
    );
    Colorbar(fig[2,4], hm2, label="m/s")


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

end