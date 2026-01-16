"""
Mostly figure generation, I just wanted to be able to run include("technical_paper.jl")
without all of this also running.
"""

function create_models()

    T = Float64
    Ndays = 5*365
    coarse_grained_hrstates = load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_3years_postspinup_dailysaves.jld2");
    uhrcg = coarse_grained_hrstates[1];
    vhrcg = coarse_grained_hrstates[2];
    etahrcg = coarse_grained_hrstates[3];

    Pnoparam = ShallowWaters.Parameter(T=T,
        output=true,
        # output_dt=168,
        output_dt=12600,
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
        nn_forcing_dissipation=false,
        N=1,
        α=2,
        nx=128,
        Ndays=Ndays
        # initial_cond="rest"
        # initpath="./dissipation_constant/spinup_files/10yearspinup_128_noslipbc_noforcing_float64prog"
        # initial_cond="ncfile",
        # initpath = "./dissipation_constant/spinup_files/1024_postspinup_noslip_5years_061824/"
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
        Ndays=Ndays
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
    onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1-3-4-6-8-10-13-15-18-23-28-30-33-35-38-41-44-46-48-51-52-53-55-58-60-63-64-65-68-73-78-83-86-88-89daystart_1dayoptimizationinitialweightsmulti3daystate_20iterations.jld2").solution
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

    ShallowWaters.time_integration(Sonline)

    # name for run
    # result_online_state_pluscD_weights_constdiffusion_20dayoptimization_startfrom20day_3years_dailysaves

end

function load_models()


    coarse_grained_hrstates = load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_3years_postspinup_dailysaves.jld2");
    uhrcg = coarse_grained_hrstates[1];
    vhrcg = coarse_grained_hrstates[2];
    etahrcg = coarse_grained_hrstates[3];

    uhr = ncread("./dissipation_const/spinup_files/1024_postspinup_30days_8hoursaves/u.nc", "u");
    vhr = ncread("./dissipation_const/spinup_files/1024_postspinup_30days_8hoursaves/v.nc", "v");
    etahr = ncread("./dissipation_const/spinup_files/1024_postspinup_30days_8hoursaves/eta.nc", "eta");

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

    u1daystategelu = ncread("./dissipation_constant/results/result_online_gelu_stateweights_1dayoptimization_startfromoffline_3years_dailysaves/u.nc", "u");
    v1daystategelu = ncread("./dissipation_constant/results/result_online_gelu_stateweights_1dayoptimization_startfromoffline_3years_dailysaves/v.nc", "v");
    eta1daystategelu = ncread("./dissipation_constant/results/result_online_gelu_stateweights_1dayoptimization_startfromoffline_3years_dailysaves/eta.nc", "eta");

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

    umulti3_old = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_3-15-30-40-50-60-80-85initdays_startfrom20daystate_3years_dailysaves/u.nc", "u");
    vmulti3_old = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_3-15-30-40-50-60-80-85initdays_startfrom20daystate_3years_dailysaves/v.nc", "v");
    etamulti3_old = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_3-15-30-40-50-60-80-85initdays_startfrom20daystate_3years_dailysaves/eta.nc", "eta");

    umulti3_new = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/u.nc", "u");
    vmulti3_new = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/v.nc", "v");
    etamulti3_new = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/eta.nc", "eta");

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

    u5s = ncread("./dissipation_constant/results/result_online_stateweights_5dayoptimization_startfrom1day_10years_weeklysaves/u.nc", "u");
    v5s = ncread("./dissipation_constant/results/result_online_stateweights_5dayoptimization_startfrom1day_10years_weeklysaves/v.nc", "v");
    eta5s = ncread("./dissipation_constant/results/result_online_stateweights_5dayoptimization_startfrom1day_10years_weeklysaves/eta.nc", "eta");

    u10s = ncread("./dissipation_constant/results/result_online_stateweights_10dayoptimization_startfrom5day_10years_weeklysaves/u.nc", "u");
    v10s = ncread("./dissipation_constant/results/result_online_stateweights_10dayoptimization_startfrom5day_10years_weeklysaves/v.nc", "v");
    eta10s = ncread("./dissipation_constant/results/result_online_stateweights_10dayoptimization_startfrom5day_10years_weeklysaves/eta.nc", "eta");

    u20s = ncread("./dissipation_constant/results/result_online_stateweights_20dayoptimization_startfrom10day_fixedcfl_10years_weeklysaves/u.nc", "u");
    v20s = ncread("./dissipation_constant/results/result_online_stateweights_20dayoptimization_startfrom10day_fixedcfl_10years_weeklysaves/v.nc", "v");
    eta20s = ncread("./dissipation_constant/results/result_online_stateweights_20dayoptimization_startfrom10day_fixedcfl_10years_weeklysaves/eta.nc", "eta");

    u30s = ncread("./dissipation_constant/results/result_online_stateweights_30dayoptimization_startfrom20day_fixedcfl_10years_weeklysaves/u.nc", "u");
    v30s = ncread("./dissipation_constant/results/result_online_stateweights_30dayoptimization_startfrom20day_fixedcfl_10years_weeklysaves/v.nc", "v");
    eta30s = ncread("./dissipation_constant/results/result_online_stateweights_30dayoptimization_startfrom20day_fixedcfl_10years_weeklysaves/eta.nc", "eta");

    umulti3_old = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_3-15-30-40-50-60-80-85initdays_startfrom20daystate_10years_weeklysaves/u.nc", "u");
    vmulti3_old = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_3-15-30-40-50-60-80-85initdays_startfrom20daystate_10years_weeklysaves/v.nc", "v");
    etamulti3_old = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_3-15-30-40-50-60-80-85initdays_startfrom20daystate_10years_weeklysaves/eta.nc", "eta");

    umulti3_new = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_10years_weeklysaves/u.nc", "u");
    vmulti3_new = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_10years_weeklysaves/v.nc", "v");
    etamulti3_new = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_10years_weeklysaves/eta.nc", "eta");

    umulti10 = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_10years_weeklysaves/u.nc", "u");
    vmulti10 = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_10years_weeklysaves/v.nc", "v");
    etamulti10 = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_10years_weeklysaves/eta.nc", "eta");

    umulti20 = ncread("./dissipation_constant/results/result_online_multistateweights_20dayoptimization_5-25-45-65initdays_startfrom20daystate_10years_weeklysaves/u.nc", "u");
    vmulti20 = ncread("./dissipation_constant/results/result_online_multistateweights_20dayoptimization_5-25-45-65initdays_startfrom20daystate_10years_weeklysaves/v.nc", "v");
    etamulti20 = ncread("./dissipation_constant/results/result_online_multistateweights_20dayoptimization_5-25-45-65initdays_startfrom20daystate_10years_weeklysaves/eta.nc", "eta");

    uzb = ncread("./dissipation_constant/spinup_files/ZB20_10yearspostspinup_weeklysaves/u.nc", "u");
    vzb = ncread("./dissipation_constant/spinup_files/ZB20_10yearspostspinup_weeklysaves/v.nc", "v");
    etazb = ncread("./dissipation_constant/spinup_files/ZB20_10yearspostspinup_weeklysaves/eta.nc", "eta");

    unoparam = ncread("./dissipation_constant/spinup_files/128_noparam_10year_postspinup_weeklysaves/u.nc", "u");
    vnoparam = ncread("./dissipation_constant/spinup_files/128_noparam_10year_postspinup_weeklysaves/v.nc", "v");
    etanoparam = ncread("./dissipation_constant/spinup_files/128_noparam_10year_postspinup_weeklysaves/eta.nc", "eta");

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
    t = 52*5
    fig = Figure(size=(700, 525), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(etahrcg[:,:,t for t in 1:522])/522,
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
    t = 366
    fig = Figure(size=(900, 800), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (unoparam[:,1:end-1,t].^2 .+ vnoparam[1:end-1,:,t].^2),
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

    # spatially averaged energy over integration (3 year integrations on all)

    oneday = []
    fiveday = []
    kespec = []
    hybrid = []
    fourier = []
    kespecpd = []
    twentyday = []
    thirtyday = []
    # relu1day = []
    # relu5day = []
    # reluKEspec = []
    noparam = []
    tenday = []
    tendayeta = []
    tendaynoeta = []
    multi3_old = []
    multi3_new = []
    multi5 = []
    multi10 = []
    multi20 = []
    zb = []
    cghr = []
    twentydaycD = []
    for j = 1:1096
        # push!(oneday, sum(u1daystategelu[:,1:end-1,j].^2 .+ v1daystategelu[1:end-1,:,j].^2))
        push!(fiveday, sum(u5s[:,1:end-1,j].^2 .+ v5s[1:end-1,:,j].^2))
        push!(tenday, sum(u10s[:,1:end-1,j].^2 .+ v10s[1:end-1,:,j].^2))
        push!(twentyday, sum(u20s[:,1:end-1,j].^2 .+ v20s[1:end-1,:,j].^2))
        push!(twentydaycD, sum(u20scD[:,1:end-1,j].^2 .+ v20scD[1:end-1,:,j].^2))
        push!(thirtyday, sum(u30s[:,1:end-1,j].^2 .+ v30s[1:end-1,:,j].^2))
        push!(kespec, sum(ukespec[:,1:end-1,j].^2 .+ vkespec[1:end-1,:,j].^2))
        push!(hybrid, sum(uhybrid[:,1:end-1,j].^2 .+ vhybrid[1:end-1,:,j].^2))
        push!(fourier, sum(ufourier[:,1:end-1,j].^2 .+ vfourier[1:end-1,:,j].^2))
        push!(kespecpd, sum(ukespecpd[:,1:end-1,j].^2 .+ vkespecpd[1:end-1,:,j].^2))
        # push!(relu1day, sum(uonline1dayrelu[:,1:end-1,j].^2 .+ vonline1dayrelu[1:end-1,:,j].^2))
        # push!(relu5day, sum(uonline5dayrelu[:,1:end-1,j].^2 .+ vonline5dayrelu[1:end-1,:,j].^2))
        # push!(reluKEspec, sum(uonlinekespecpdrelu[:,1:end-1,j].^2 .+ vonlinekespecpdrelu[1:end-1,:,j].^2))
        push!(zb, sum(uzb[:,1:end-1,j].^2 .+ vzb[1:end-1,:,j].^2))
        push!(cghr, sum(uhrcg[:,1:end-1,j].^2 .+ vhrcg[1:end-1,:,j].^2))
        push!(noparam, sum(unoparam[:,1:end-1,j].^2 .+ vnoparam[1:end-1,:,j].^2))
        # push!(multi3_old, sum(umulti3_old[:,1:end-1,j].^2 .+ vmulti3_old[1:end-1,:,j].^2))
        push!(multi3_new, sum(umulti3_new[:,1:end-1,j].^2 .+ vmulti3_new[1:end-1,:,j].^2))
        push!(multi5, sum(umulti5[:,1:end-1,j].^2 .+ vmulti5[1:end-1,:,j].^2))
        push!(multi10, sum(umulti10[:,1:end-1,j].^2 .+ vmulti10[1:end-1,:,j].^2))
        push!(multi20, sum(umulti20[:,1:end-1,j].^2 .+ vmulti20[1:end-1,:,j].^2))
    end

    # 3 year figure
    fig = Figure(size=(1000, 500), fontsize=15);
    ax = Axis(fig[1,1],
            xlabel="Day",
            ylabel="Energy",
            title="Spatially averaged energy over 3 years"
    )
    lines!(ax, LinRange(0, 3*365, 1096),  cghr ./ (128^2), label="Coarse-grained HR")
    lines!(ax, LinRange(0, 3*365, 1096), noparam./ (128^2), label="30 km resolution, no closure")
    lines!(ax, LinRange(0, 3*365, 1096), zb./ (128^2), label="ZB20")
    # lines!(fig[1,1], LinRange(0, 3*365, 1096), fiveday./ (128^2), label="Online closure, 5 day")
    lines!(ax, LinRange(0, 3*365, 1096), tenday./ (128^2), label="Online closure, 10 day")
    lines!(ax,LinRange(0, 3*365, 1096), twentyday ./ (128^2), label="Online closure, 20 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentydaycD ./ (128^2), label="Online closure, 20 day with BD coeff")
    lines!(ax, LinRange(0, 3*365, 1096), thirtyday ./ (128^2), label="Online closure, 30 day")
    lines!(ax, LinRange(0, 3*365, 1096), multi3_new ./ (128^2), label="Online closure, multi 3 day")
    lines!(ax, LinRange(0, 3*365, 1096), multi5 ./ (128^2), label="Online closure, multi 5 day", color=:mediumorchid)
    lines!(ax, LinRange(0, 3*365, 1096), multi10 ./ (128^2), label="Online closure, multi 10 day", color=:teal)
    lines!(ax, LinRange(0, 3*365, 1096), multi20 ./ (128^2), label="Online closure, multi 20 day", color=:red3)
    Legend(fig[1, 2], ax)

    # 10 year figure
    fig = Figure(size=(1000, 500), fontsize=15);
    lines(fig[1,1], LinRange(0, 10*365, 522),  zb ./ (128^2), label="ZB20", 
        axis=(
            xlabel="Day",
            ylabel="Energy",
            title="Spatially averaged energy over 10 years"
        )
    )
    lines!(fig[1,1], LinRange(0, 10*365, 522), noparam ./ (128^2), label="30 km resolution, no closure")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), fiveday./ (128^2), label="Online closure, 5 day")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), tenday./ (128^2), label="Online closure, 10 day")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), twentyday./ (128^2), label="Online closure, 20 day")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), thirtyday./ (128^2), label="Online closure, 30 day")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), multi3_old ./ (128^2), label="Online closure, multi 3 day, fewer initial conditions")
    lines!(fig[1,1], LinRange(0, 10*365, 522), multi3_new ./ (128^2), label="Online closure, multi 3 day", color=:mediumorchid)
    lines!(fig[1,1], LinRange(0, 10*365, 522), multi10 ./ (128^2), label="Online closure, multi 10 day", color=:teal)
    # lines!(fig[1,1], LinRange(0, 10*365, 522), multi20 ./ (128^2), label="Online closure, multi 20 day", color=:red3)
    axislegend(position = (0,1))

    fig = Figure(size=(1000, 500), fontsize=15);
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

end

function spectrum_plots()

    # KE spectrum #############################################################

    # to get coarse-grained states
    ker = ImageFiltering.Kernel.gaussian((30e3/3750))
    # imfilter(hru[:,:,j], reflect(ker))

    uhr = ncread("./spinup_files/1024_postspinup_30days_8hoursaves/u.nc", "u");
    vhr = ncread("./spinup_files/1024_postspinup_30days_8hoursaves/v.nc", "v");
    etahr = ncread("./spinup_files/1024_postspinup_30days_8hoursaves/eta.nc", "eta");

    coarse_grained_hrstates = load_object("./spinup_files/1024_filtered_downsized_uveta_3years_postspinup_dailysaves.jld2");
    uhrcg = coarse_grained_hrstates[1];
    vhrcg = coarse_grained_hrstates[2];
    etahrcg = coarse_grained_hrstates[3];

    # coarse_grained_hrstates = load_object("./offline_files/1024_filtered_downsized_uveta_30days_postspinup_8hoursaves_112125.jld2");
    # uhrcg = coarse_grained_hrstates[1];
    # vhrcg = coarse_grained_hrstates[2];
    # etahrcg = coarse_grained_hrstates[3];

    filtered_hrstates = load_object("./offline_files/1024_filtered_uveta_imfilter_30days_postspinup_8hoursaves_112125.jld2");
    uhrfilter= filtered_hrstates[1];
    vhrfilter = filtered_hrstates[2];
    etahrfilter = filtered_hrstates[3];

    totalstates = 1096 # saved every 8 hours (this is when the hr cg and low resolution match up)
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

    up_multi3_old = zeros(65,totalstates)
    vp_multi3_old = zeros(65,totalstates)

    up_multi3_new = zeros(65,totalstates)
    vp_multi3_new = zeros(65,totalstates)

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

    for t = 1:totalstates

        # up_hr[:,t] = power(periodogram(uhr[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2
        # vp_hr[:,t] = power(periodogram(vhr[:,:,t]; radialavg=true, radialsum=false)) ./ 1024^2

        up_zb[:,t] = power(periodogram(uzb[:, :, t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_zb[:,t] = power(periodogram(vzb[:, :, t]; radialavg=true, radialsum=false)) ./ 128^2

        up_hrcg[:,t] = power(periodogram(uhrcg[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_hrcg[:,t] = power(periodogram(vhrcg[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_20day[:,t] = power(periodogram(u20s[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_20day[:,t] = power(periodogram(v20s[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_30day[:,t] = power(periodogram(u30s[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_30day[:,t] = power(periodogram(v30s[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_noparam[:,t] = power(periodogram(unoparam[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_noparam[:,t] = power(periodogram(vnoparam[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_gelu5day[:,t] = power(periodogram(u5s[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_gelu5day[:,t] = power(periodogram(v5s[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_multi10[:,t] = power(periodogram(umulti10[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_multi10[:,t] = power(periodogram(vmulti10[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_multi20[:,t] = power(periodogram(umulti20[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_multi20[:,t] = power(periodogram(vmulti20[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_multi3_old[:,t] = power(periodogram(umulti3_old[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_multi3_old[:,t] = power(periodogram(vmulti3_old[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

        up_multi3_new[:,t] = power(periodogram(umulti3_new[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2
        vp_multi3_new[:,t] = power(periodogram(vmulti3_new[:,:,t]; radialavg=true, radialsum=false)) ./ 128^2

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
    lines(fig[1,1], lr_wl[2:65], up_hrcg[2:65,t] + vp_hrcg[2:65,t], label="Filtered, coarse-grained 3.75km resolution", axis=(
            xscale=log10,
            yscale=log10,
            xlabel="Wavelength (km)",
            ylabel="KE(k)",
            xreversed=true,
            xticks=[700, 100, 30, 10, 2],
            title="KE spectrum")
    )
    # lines!(fig[1,1], lr_wl[2:end], up_gelu1day[2:end,t] + vp_gelu1day[2:end,t], label="Online NN closure, 1 day")
    lines!(fig[1,1], lr_wl[2:end], up_gelu5day[2:end,t] + vp_gelu5day[2:end,t], label="Online NN closure, 5 day")
    lines!(fig[1,1], lr_wl[2:end], up_20day[2:end,t] + vp_20day[2:end,t], label="Online NN closure, 20 day")
    # lines!(fig[1,1], lr_wl[2:end], up_geluKEspecpd[2:end,t] + up_geluKEspecpd[2:end,t], label="Online NN closure, KE spectrum PD")
    lines!(fig[1,1], lr_wl[2:end], up_noparam[2:end,t] + vp_noparam[2:end,t], label="30 km resolution, no closure")
    lines!(fig[1,1], lr_wl[2:end], up_zb[2:end,t] + vp_zb[2:end,t], label="ZB closure")
    lines!(fig[1,1], lr_wl[2:end], up_multi10[2:end,t] + vp_multi10[2:end,t], label="Online NN closure, multi 10 day")
    lines!(fig[1,1], lr_wl[2:end], up_multi20[2:end,t] + vp_multi20[2:end,t], label="Online NN closure, multi 20 day")
    lines!(fig[1,1], lr_wl[2:end], up_multi3_old[2:end,t] + vp_multi3_old[2:end,t], label="Online NN closure, multi 3 day old")
    lines!(fig[1,1], lr_wl[2:end], up_multi3_new[2:end,t] + vp_multi3_new[2:end,t], label="Online NN closure, multi 3 day new")

    axislegend(position = (0,0))

    #### comparing time-averaged ke spectra

    up_noparam_avg = zeros(65)
    vp_noparam_avg = zeros(65)

    up_zb_avg = zeros(65)
    vp_zb_avg = zeros(65)

    up_10day_avg = zeros(65)
    vp_10day_avg = zeros(65)
    
    up_10dayeta_avg = zeros(65)
    vp_10dayeta_avg = zeros(65)

    up_5dayeta_avg = zeros(65)
    vp_5dayeta_avg = zeros(65)

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

    up_gelu5day_avg = zeros(65)
    vp_gelu5day_avg = zeros(65)

    up_multi10_avg = zeros(65)
    vp_multi10_avg = zeros(65)

    up_multi20_avg = zeros(65)
    vp_multi20_avg = zeros(65)

    up_multi3_old_avg = zeros(65)
    vp_multi3_old_avg = zeros(65)

    up_multi3_new_avg = zeros(65)
    vp_multi3_new_avg = zeros(65)

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

        up_gelu5day_avg += up_gelu5day[:,t]
        vp_gelu5day_avg += vp_gelu5day[:,t]

        up_10day_avg += up_10day[:, t]
        vp_10day_avg += vp_10day[:, t]

        up_multi20_avg += up_multi20[:, t]
        vp_multi20_avg += vp_multi20[:, t]

        up_multi10_avg += up_multi10[:, t]
        vp_multi10_avg += vp_multi10[:, t]

        up_multi3_old_avg += up_multi3_old[:, t]
        vp_multi3_old_avg += vp_multi3_old[:, t]

        up_multi3_new_avg += up_multi3_new[:, t]
        vp_multi3_new_avg += vp_multi3_new[:, t]

        # up_5dayeta_avg += up_5dayeta[:, t]
        # vp_5dayeta_avg += vp_5dayeta[:, t]

        # up_10dayeta_avg += up_10dayeta[:, t]
        # vp_10dayeta_avg += vp_10dayeta[:, t]
        # up_geluKEspec_avg += up_geluKEspec[:,t]
        # vp_geluKEspec_avg += vp_geluKEspec[:,t]

        # up_geluKEspecpd_avg += up_geluKEspecpd[:,t]
        # vp_geluKEspecpd_avg += vp_geluKEspecpd[:,t]

        # up_geluhybrid_avg += up_geluhybrid[:,t]
        # vp_geluhybrid_avg += vp_geluhybrid[:,t]

        # up_gelufourier_avg += up_gelufourier[:,t]
        # vp_gelufourier_avg += vp_gelufourier[:,t]

    end

    fig = Figure(size=(1000, 500), fontsize=15);
    lines(fig[1,1], lr_wl[2:65], (up_cghr_avg[2:65] + vp_cghr_avg[2:65])/1098, label="Filtered, coarse-grained 3.75km resolution", axis=(
        xscale=log10,
        yscale=log10,
        xlabel="Wavelength (km)",
        ylabel="KE(k)", xreversed=true,
        xticks=[700, 100, 30, 10, 2],
        title="Time-averaged KE spectrum, 3 year spinup")
    )
    lines!(fig[1,1], lr_wl[2:end], (up_zb_avg[2:end] + vp_zb_avg[2:end])/1098, label="ZB20")
    lines!(fig[1,1], lr_wl[2:end], (up_noparam_avg[2:end] + vp_noparam_avg[2:end])/1098, label="30 km resolution, no closure")
    # lines!(fig[1,1], lr_wl[2:end], (up_relu1day_avg[2:end] + vp_gelu1day_avg[2:end])/totalstates, label="Online NN closure, 1 day")
    # lines!(fig[1,1], lr_wl[2:end], (up_gelu5day_avg[2:end] + vp_gelu5day_avg[2:end])/1096, label="Online NN closure, 5 day")
    # lines!(fig[1,1], lr_wl[2:end], (up_10day_avg[2:end] + vp_10day_avg[2:end])/1096, label="Online NN closure, 10 day")
    # lines!(fig[1,1], lr_wl[2:end], (up_20day_avg[2:end] + vp_20day_avg[2:end])/1096, label="Online NN closure, 20 day")
    lines!(fig[1,1], lr_wl[2:end], (up_30day_avg[2:end] + vp_30day_avg[2:end])/1096, label="Online NN closure, 30 day", linestyle=:dashdot)
    lines!(fig[1,1], lr_wl[2:end], (up_multi3_new_avg[2:end] + vp_multi3_new_avg[2:end])/1096, label="Online NN closure, multi 3 day", linestyle=:dash)
    lines!(fig[1,1], lr_wl[2:end], (up_multi10_avg[2:end] + vp_multi10_avg[2:end])/1096, label="Online NN closure, multi 10 day",linestyle=:dot)
    # lines!(fig[1,1], lr_wl[2:end], (up_multi20_avg[2:end] + vp_multi20_avg[2:end])/1096, label="Online NN closure, multi 20 day", color=:red)

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

    uhr = ncread("./dissipation_constant/spinup_files/1024_postspinup_noslip_5years_061824/u.nc", "u");
    vhr = ncread("./dissipation_constant/spinup_files/1024_postspinup_noslip_5years_061824/v.nc", "v");

    coarse_grained_hrstates = load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_3years_postspinup_dailysaves.jld2");
    uhrcg = coarse_grained_hrstates[1];
    vhrcg = coarse_grained_hrstates[2];
    etahrcg = coarse_grained_hrstates[3];

    Suhr = load_object("./dissipation_constant/results/true_S_SuSv_fromtimederivatives_new.jld2")[1]
    Svhr = load_object("./dissipation_constant/results/true_S_SuSv_fromtimederivatives_new.jld2")[2]

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

    lr_freq = 1/30 .* freq(periodogram(u20day[:,:,10]; radialavg=true, radialsum=false));
    nfft = nextfastfft(size(uhrcg[:,:,1]))

    totalu_hrcg = zeros(65)
    totalv_hrcg = zeros(65)

    totalu_5 = zeros(65)
    totalv_5 = zeros(65)

    totalu_10 = zeros(65)
    totalv_10 = zeros(65)

    totalu_20 = zeros(65)
    totalv_20 = zeros(65)
    
    totalu_30 = zeros(65)
    totalv_30 = zeros(65)

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

    totalstates = 1097
    for t = 1:totalstates

        # Suhr, Svhr = compute_hrS(uhr[:,:,t], vhr[:,:, t])

        # u5, v5, eta5 = ShallowWaters.add_halo(Float64.(u5day[:,:,t]), Float64.(v5day[:,:,t]), Float64.(eta5day[:,:,t]), zeros(128,128), S5);
        u, v, eta = ShallowWaters.add_halo(uhrcg[:,:,t], vhrcg[:,:,t], etahrcg[:,:,t], zeros(128,128), Shrcg);
        # uzb, vzb, etazb = ShallowWaters.add_halo(Float64.(uzb_[:,:,t]), Float64.(vzb_[:,:,t]), Float64.(etazb_[:,:,t]), zeros(128,128), SZB);
        # u20, v20, _ = ShallowWaters.add_halo(Float64.(u20day[:,:,t]), Float64.(v20day[:,:,t]), Float64.(eta20day[:,:,t]), zeros(128,128), S20);
        # u30, v30, _ = ShallowWaters.add_halo(Float64.(u30day[:,:,t]), Float64.(v30day[:,:,t]), Float64.(eta30day[:,:,t]), zeros(128,128), S30);
        # umulti3_, vmulti3_, _ = ShallowWaters.add_halo(Float64.(umulti3[:,:,t]), Float64.(vmulti3[:,:,t]), Float64.(etamulti3[:,:,t]), zeros(128,128), Smulti3);
        # umulti5_, vmulti5_, _ = ShallowWaters.add_halo(Float64.(umulti5[:,:,t]), Float64.(vmulti5[:,:,t]), Float64.(etamulti5[:,:,t]), zeros(128,128), Smulti5);

        # umulti10_, vmulti10_, _ = ShallowWaters.add_halo(Float64.(umulti10[:,:,t]), Float64.(vmulti10[:,:,t]), Float64.(etamulti10[:,:,t]), zeros(128,128), Smulti10);
        # umulti20_, vmulti20_, _ = ShallowWaters.add_halo(Float64.(umulti20[:,:,t]), Float64.(vmulti20[:,:,t]), Float64.(etamulti20[:,:,t]), zeros(128,128), Smulti20);

        # ShallowWaters.CNN_momentum(u20, v20, S20)
        # ShallowWaters.CNN_momentum(u30, v30, S30)
        # ShallowWaters.CNN_momentum(umulti3_, vmulti3_, Smulti3)
        # ShallowWaters.CNN_momentum(umulti5_, vmulti5_, Smulti5)
        # ShallowWaters.CNN_momentum(umulti10_, vmulti10_, Smulti10)
        # ShallowWaters.CNN_momentum(umulti20_, vmulti20_, Smulti20)

        # ShallowWaters.CNN_momentum(u, v, Shrcg);
        # ShallowWaters.CNN_momentum(u5, v5, S5);
        # ShallowWaters.CNN_momentum(u20, v20, S20);
        # ShallowWaters.ZB_momentum(uzb, vzb, SZB, SZB.Diag);

        outu_hrcg, inputu_hrcg, inputSu_hrcg = paddingu(uhrcg[:, :, t+1], Suhr[:,:,t], nfft[1])
        fft2pow2radial!(outu_hrcg, rfft(inputu_hrcg), rfft(inputSu_hrcg), nfft...)
        outv_hrcg, inputv_hrcg, inputSv_hrcg = paddingv(vhrcg[:, :, t+1], Svhr[:,:,t], nfft[1])
        fft2pow2radial!(outv_hrcg, rfft(inputv_hrcg), rfft(inputSv_hrcg), nfft...)

        totalu_hrcg += outu_hrcg
        totalv_hrcg += outv_hrcg

        # outu_5, inputu_5, inputSu_5 = paddingu(u5day[:, :, t], S5.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_5, rfft(inputu_5), rfft(inputSu_5), nfft...)
        # outv_5, inputv_5, inputSv_5 = paddingv(v5day[:, :, t], S5.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_5, rfft(inputv_5), rfft(inputSv_5), nfft...)

        # totalu_5 += outu_5
        # totalv_5 += outv_5

        # outu_20, inputu_20, inputSu_20 = paddingu(u20day[:, :, t], S20.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_20, rfft(inputu_20), rfft(inputSu_20), nfft...)
        # outv_20, inputv_20, inputSv_20 = paddingv(v20day[:, :, t], S20.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_20, rfft(inputv_20), rfft(inputSv_20), nfft...)

        # totalu_20 += outu_20
        # totalv_20 += outv_20

        # outu_30, inputu_30, inputSu_30 = paddingu(u30day[:, :, t], S30.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_30, rfft(inputu_30), rfft(inputSu_30), nfft...)
        # outv_30, inputv_30, inputSv_30 = paddingv(v30day[:, :, t], S30.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_30, rfft(inputv_30), rfft(inputSv_30), nfft...)

        # totalu_30 += outu_30
        # totalv_30 += outv_30

        # outu_multi3, inputu_multi3, inputSu_multi3 = paddingu(umulti3[:, :, t], Smulti3.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_multi3, rfft(inputu_multi3), rfft(inputSu_multi3), nfft...)
        # outv_multi3, inputv_multi3, inputSv_multi3 = paddingv(vmulti3[:, :, t], Smulti3.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_multi3, rfft(inputv_multi3), rfft(inputSv_multi3), nfft...)

        # totalu_multi3 += outu_multi3
        # totalv_multi3 += outv_multi3

        # outu_multi5, inputu_multi5, inputSu_multi5 = paddingu(umulti5[:, :, t], Smulti5.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_multi5, rfft(inputu_multi5), rfft(inputSu_multi5), nfft...)
        # outv_multi5, inputv_multi5, inputSv_multi5 = paddingv(vmulti5[:, :, t], Smulti5.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_multi5, rfft(inputv_multi5), rfft(inputSv_multi5), nfft...)

        # totalu_multi5 += outu_multi5
        # totalv_multi5 += outv_multi5

        # outu_multi10, inputu_multi10, inputSu_multi10 = paddingu(umulti10[:, :, t], Smulti10.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_multi10, rfft(inputu_multi10), rfft(inputSu_multi10), nfft...)
        # outv_multi10, inputv_multi10, inputSv_multi10 = paddingv(vmulti10[:, :, t], Smulti10.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_multi10, rfft(inputv_multi10), rfft(inputSv_multi10), nfft...)

        # totalu_multi10 += outu_multi10
        # totalv_multi10 += outv_multi10

        # outu_multi20, inputu_multi20, inputSu_multi20 = paddingu(umulti20[:, :, t], Smulti20.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_multi20, rfft(inputu_multi20), rfft(inputSu_multi20), nfft...)
        # outv_multi20, inputv_multi20, inputSv_multi20 = paddingv(vmulti20[:, :, t], Smulti20.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_multi20, rfft(inputv_multi20), rfft(inputSv_multi20), nfft...)

        # totalu_multi20 += outu_multi20
        # totalv_multi20 += outv_multi20

        # outu_3, inputu_3, inputSu_3 = paddingu(umulti3[:, :, t], S3.Diag.CNNVars.S_u, nfft[1])
        # fft2pow2radial!(outu_3, rfft(inputu_3), rfft(inputSu_3), nfft...)
        # outv_3, inputv_3, inputSv_3 = paddingv(vmulti3[:, :, t], S3.Diag.CNNVars.S_v, nfft[1])
        # fft2pow2radial!(outv_3, rfft(inputv_3), rfft(inputSv_3), nfft...)

        # totalu_3 += outu_3
        # totalv_3 += outv_3

        # outu_ZB, inputu_ZB, inputSu_ZB = paddingu(uzb_[:, :, t], SZB.Diag.ZBVars.S_u, nfft[1])
        # fft2pow2radial!(outu_ZB, rfft(inputu_ZB), rfft(inputSu_ZB), nfft...)
        # outv_ZB, inputv_ZB, inputSv_ZB = paddingv(vzb_[:, :, t], SZB.Diag.ZBVars.S_v, nfft[1])
        # fft2pow2radial!(outv_ZB, rfft(inputv_ZB), rfft(inputSv_ZB), nfft...)

        # totalu_ZB += outu_ZB
        # totalv_ZB += outv_ZB

    end
end

function ketransfer_plots()

    lr_freq = 1/30 .* freq(periodogram(u20s[:,:,10]; radialavg=true, radialsum=false));
    nfft = nextfastfft(size(uhrcg[:,:,1]))

    # totalu_hrcg = load_object("./dissipation_constant/results/hrcg_ketransfer_uv_3years.jld2")[1];
    # totalv_hrcg = load_object("./dissipation_constant/results/hrcg_ketransfer_uv_3years.jld2")[2];

    totalu_hrcg = load_object("./dissipation_constant/results/hrcg_ketransfer_fromtimederivative_uv_3years.jld2")[1];
    totalv_hrcg = load_object("./dissipation_constant/results/hrcg_ketransfer_fromtimederivative_uv_3years.jld2")[2];

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

    fig = Figure(size=(800, 300), fontsize=15);
    ax = Axis(fig[1,1],
        xscale = log10,
        xlabel="Wavenumber (1/km)",
        ylabel="KE(k)",
        title="Kinetic energy transfer"
    )
    lines!(ax, lr_freq.*(totalu_hrcg + totalv_hrcg) / 1096, label="Subgrid forcing")
    lines!(ax, lr_freq.*(totalu_ZB + totalv_ZB) / 1098, label="ZB20")
    lines!(ax, lr_freq.*(totalu_30 + totalv_30) / 1096, label="30 day optimization")
    lines!(ax, lr_freq.*(totalu_20 + totalv_20) / 1096, label="20 day optimization")
    # lines!(ax, lr_freq.*(totalu_multi3_old + totalv_multi3_old) / 1096, label="Multi 3 day state optimization old")
    lines!(ax, lr_freq.*(totalu_multi3_new + totalv_multi3_new) / 1096, label="Multi 3 day state optimization")
    lines!(ax, lr_freq.*(totalu_multi5 + totalv_multi5) / 1096, label="Multi 5 day state optimization")
    lines!(ax, lr_freq.*(totalu_multi10 + totalv_multi10) / 1096, label="Multi 10 day state optimization")
    lines!(ax, lr_freq.*(totalu_multi20 + totalv_multi20) / 1096, label="Multi 20 day state optimization")
    Legend(fig[1,2], ax)

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