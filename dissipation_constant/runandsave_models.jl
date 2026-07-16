function create_models()

    T = Float64
    Ndays = 3*365

    Shr = ShallowWaters.model_setup(T=T,
        output=true,
        output_vars=["u", "v", "η", "ζ", "du", "dv"],
        # output_dt = 1,
        # output_dt=168,
        # output_dt=12600,
        L_ratio=1,
        g=9.81,
        H=500,
        ϕ = 50.,
        wind_forcing_x="double_gyre",
        Fx0=0.12,
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
        nx=1024,
        Ndays=3*365-367,
        initial_cond="ncfile",
        initpath="./dissipation_constant/generalizability/1024_postspinup_newlatitude/1024_postspinup_newlatitude_days1-367"
    );

    ShallowWaters.time_integration(Shr)

    coarse_grained_hrstates = load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_imfilter_3years_postspinup_dailysaves_correctedsetup.jld2");
    uhrcg = coarse_grained_hrstates[1]
    vhrcg = coarse_grained_hrstates[2]
    etahrcg = coarse_grained_hrstates[3]

    # coarse_grained_hrstates = load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_90days_postspinup_8hoursaves.jld2");
    # uhrcg2 = coarse_grained_hrstates[1];

    Pnoparam = ShallowWaters.Parameter(T=T,
        output=true,
        output_vars=["u", "v", "η", "ζ"],
        # output_dt = 1,
        output_dt=168,
        # output_dt=12600,
        L_ratio=1,
        g=9.81,
        H=500,
        cfl=.898,
        ϕ = 45.,
        wind_forcing_x="double_gyre",
        Fx0=0.12,
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
        Ndays=Ndays,
        initial_cond="rest"
        # initpath="./dissipation_constant/spinup_files/1024_spinup_newlatandwindamp"
        # initpath="./dissipation_constant/spinup_files/1024_spinup_newwindamp/1024_spinup_newwindamp_days374-1095"
        # initpath="./dissipation_constant/spinup_files/1024_spinup_newlatitude/1024_spinup_newlatitude_days381-1095",
        # initpath="./dissipation_constant/spinup_files/10yearspinup_128_noslipbc_noforcing_float64prog"
        # initpath = "./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves/1024_postspinup_day1-766saves"
    );

    Snoparam = ShallowWaters.model_setup(Pnoparam);

    u0, v0, eta0, _ = ShallowWaters.add_halo(uhrcg[:,:,1],vhrcg[:,:,1],etahrcg[:,:,1],zeros(128,128),Snoparam);
    # u0, v0, eta0, _ = ShallowWaters.add_halo(uhrcg,vhrcg,etahrcg,zeros(128,128),Snoparam);

    initial_cond = [u0, v0, eta0];

    Snoparam.Prog.u .= copy(initial_cond[1]);
    Snoparam.Prog.v .= copy(initial_cond[2]);
    Snoparam.Prog.η .= copy(initial_cond[3]);

    ShallowWaters.time_integration(Snoparam);

    Ndays = 3*365
    PZB = ShallowWaters.Parameter(T=T,
        output=true,
        output_vars=["u", "v", "η", "ζ"],
        output_dt=24,
        L_ratio=1,
        g=9.81,
        H=500,
        ϕ = 45.,
        wind_forcing_x="double_gyre",
        Fx0=1.2,
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
        output_vars=["u", "v", "η", "ζ"],
        L_ratio=1,
        g=9.81,
        H=500,
        cfl=.898,
        wind_forcing_x="double_gyre",
        Fx0=0.12,
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

    offlineweights = load_object("./offline_relu_newercheck.jld2");
    # offlineweights = load_object("./dissipation_constant/tuned_weights/result_offline_150iterations_reluactivation_111925.jld2").solution
    # offlineweights = load_object("./dissipation_constant/tuned_weights/result_offline_150iterations_geluactivation_111925.jld2").solution
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
        output_vars=["u", "v", "η", "ζ"],
        # output_dt = 1,
        output_dt=168,
        # output_dt=12600,
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
        Ndays=Ndays,
        initial_cond="ncfile",
        initpath="./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_10years_weeklysaves"
    );

    Sonline = ShallowWaters.model_setup(Ponline);

    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_3-15-30-40-50-60-80-85daystart_3dayoptimization_initialweights20daystate_30iterations.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_madnlp_states_5dayoptimization_startfrom1daystate_50iterations_geluactivation.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_state_10dayoptimization_startfrom5daystate_noeta_30iterations.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_state_30dayoptimzation_startfrom20day_constantdissipation_6iterations_8hourdata_200maxhistory_fixedcfl.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_state_30dayoptimzation_startfrom30day6iterations_constantdissipation_15iterations_21totaliterations_8hourdata_200maxhistory_fixedcfl.jld2").solution;
    # onlineweights = load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_state_20dayoptimzation_startfrom10day_constantdissipation_10iterations_8hourdata_200maxhistory_fixedcfl.jld2").solution

    onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_5-20-35-50-65-75daystart_10dayoptimization_initialweights20daystate_fixedcfl_15iterations_constdissipation.jld2").solution
    # onlineweights = load_object("./result_multistate_5-25-45-65daystart_20dayoptimization_initialweights20daystate_fixedcfl_15iterations_constdissipation.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1-3-4-6-8-10-13-15-18-23-28-30-33-35-38-41-44-46-48-51-52-53-55-58-60-63-64-65-68-73-78-83-86-88-89daystart_1dayoptimizationinitialweightsmulti3daystate_20iterations.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1:2:89initdaystart_1dayoptimization_initialweightsmulti3daystate_fewerinitconds_15iterations.jld2").solution;

    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86daystart_3dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-6-8-10-13-15-18-23-28-33-38-41-44-48-51-53-58-63-65-68-73-78-83-86-88daystart_2dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution;
    # onlineweights = load_object("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/average_allresults.jld2")
    # onlineweights = load_object("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/average_allresults.jld2")
    
    # onlineweights = load_object("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/result_online_manyinitconds_initday40_10dayoptimzation_startfrom20dayoptimization_constantdissipation_25iterations_8hourdata_200maxhistory_fixedcfl.jld2").solution;
    # onlineweights = load_object("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/result_online_manyinitconds_initday33_3dayoptimzation_startfrom20dayoptimization_constantdissipation_40iterations_8hourdata_200maxhistory_fixedcfl.jld2").solution;
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

    # Sonline.Prog.u .= copy(initial_cond[1]);
    # Sonline.Prog.v .= copy(initial_cond[2]);
    # Sonline.Prog.η .= copy(initial_cond[3]);

    P = ShallowWaters.time_integration(Sonline);

    # name for run
    # result_multistate_10day_further10years_weeklysaves

end

function online_runs(j, init_day)

    T = Float64
    Ndays = 10 * 365

    coarse_grained_hrstates = load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_3years_postspinup_dailysaves.jld2");
    uhrcg = coarse_grained_hrstates[1]
    vhrcg = coarse_grained_hrstates[2]
    etahrcg = coarse_grained_hrstates[3]

    Pnoparam = ShallowWaters.Parameter(T=T,
        output=true,
        output_vars=["u", "v", "η", "ζ"],
        # output_dt = 1,
        output_dt=168,
        # output_dt=12600,
        L_ratio=1,
        g=9.81,
        H=500,
        cfl=.898,
        ϕ = 45.,
        wind_forcing_x="double_gyre",
        Fx0=0.12,
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
        Ndays=Ndays,
        initial_cond="rest"
    );

    Snoparam = ShallowWaters.model_setup(Pnoparam);

    u0, v0, eta0, _ = ShallowWaters.add_halo(uhrcg[:,:,init_day],vhrcg[:,:,init_day],etahrcg[:,:,init_day],zeros(128,128),Snoparam);

    initial_cond = [u0, v0, eta0];

    Ponline = ShallowWaters.Parameter(T=T,
        output=true,
        output_vars=["u", "v", "η", "ζ"],
        output_dt=168,
        L_ratio=1,
        g=9.81,
        H=500,
        cfl=.898,
        ϕ = 45.,
        wind_forcing_x="double_gyre",
        Fx0=0.12,
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

    Sonline = ShallowWaters.model_setup(Ponline);

    # initial days used for the 3 day optimizations
    # init_ts = [1,3,8,13,18,23,28,33,38,41,44,48,53,58,63,68,73,78,83,86]

    # inital days used for the 10 day optimizations
    init_ts = [2, 20, 30, 40, 50, 60, 70, 80]
    t = init_ts[j]
    onlineweights = load_object(make_filename2(t)).solution;

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
end

function make_filename2(t)
    # return "./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/" *
    #     "result_online_manyinitconds_initday$(t)_3dayoptimzation_startfrom20dayoptimization_constantdissipation_40iterations_8hourdata_200maxhistory_fixedcfl.jld2"
    return "./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/" *
          "result_online_manyinitconds_initday$(t)_10dayoptimzation_startfrom20dayoptimization_constantdissipation_25iterations_8hourdata_200maxhistory_fixedcfl.jld2"
end