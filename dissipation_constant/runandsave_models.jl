function create_models()

    T = Float64
    Ndays = 5*365

    Shr = ShallowWaters.model_setup(T=T,
        output=true,
        output_vars=["u", "v", "η", "ζ", "du", "dv"],
        # output_dt = 1,
        # output_dt=168,
        # output_dt=12600,
        L_ratio=1,
        g=9.81,
        H=500,
        ϕ = -60,
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
        Ndays=4*365,
        initial_cond="ncfile",
        initpath="./dissipation_constant/generalizability/1024_4yearspinup_latmin60"
    );

    ShallowWaters.time_integration(Shr)

    uhrcg = uhrcgall;
    vhrcg = vhrcgall;
    etahrcg = etahrcgall;

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

    Ndays = 10*365
    PZB = ShallowWaters.Parameter(T=T,
        output=true,
        output_vars=["u", "v", "η", "ζ"],
        # output_dt=168,
        L_ratio=1,
        g=9.81,
        H=500,
        ϕ = 45.,
        wind_forcing_x="double_gyre",
        Fx0=.12,
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
        Ndays=Ndays,
        initial_cond="rest"
        # initpath="./dissipation_constant/spinup_files/ZB20_10yearspostspinup_weeklysaves"
    );

    SZB = ShallowWaters.model_setup(PZB);
    SZB.Prog.u .= copy(initial_cond[1]);
    SZB.Prog.v .= copy(initial_cond[2]);
    SZB.Prog.η .= copy(initial_cond[3]);

    ShallowWaters.time_integration(SZB);

    # the offline problem is unstable even for short integrations, Ndays should be capped around 3 days
    Ndays = 3
    Poffline = ShallowWaters.Parameter(T=T,
        output=true,
        output_dt=8,
        output_vars=["u", "v", "η", "ζ"],
        L_ratio=1,
        g=9.81,
        H=500,
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

    offlineweights = load_object("./dissipation_constant/tuned_weights/result_offline_150iterations_geluactivation_111925.jld2").solution
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
    Ndays = 5
    Ponline = ShallowWaters.Parameter(T=T,
        output=true,
        output_vars=["u", "v", "η", "ζ"],
        # output_dt = 1,
        # output_dt=168,
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
        initial_cond="rest"
    );

    Sonline = ShallowWaters.model_setup(Ponline);

    # onlineweights = load_object("./dissipation_constant/tuned_weights/single_initial_condition/result_online_states_5dayoptimization_startfrom1daystate_50iterations_geluactivation.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/single_initial_condition/result_online_states_10dayoptimization_startfrom5daystate_noeta_30iterations.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/single_initial_condition/result_online_state_20dayoptimzation_startfrom10day_constantdissipation_10iterations_8hourdata_200maxhistory_fixedcfl.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/single_initial_condition/result_online_state_30dayoptimzation_startfrom20day_constantdissipation_6iterations_8hourdata_200maxhistory_fixedcfl.jld2").solution

    # onlineweights = load_object("./dissipation_constant/tuned_weights/ensemble_initial_condition/result_multistate_1-4-6-8-10-13-15-18-23-28-33-38-41-44-48-51-53-58-63-65-68-73-78-83-86-88daystart_2dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution;
    # onlineweights = load_object("./dissipation_constant/tuned_weights/ensemble_initial_condition/result_multistate_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86daystart_3dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/ensemble_initial_condition/result_multistate_5-20-35-50-65-75daystart_10dayoptimization_initialweights20daystate_fixedcfl_15iterations_constdissipation.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/ensemble_initial_condition/result_multistate_5-25-45-65daystart_20dayoptimization_initialweights20daystate_fixedcfl_15iterations_constdissipation.jld2").solution

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

    Sonline.Prog.u .= u0;
    Sonline.Prog.v .= v0;
    Sonline.Prog.η .= eta0;

    P = ShallowWaters.time_integration(Sonline);

    # name for run
    # result_multistate_10day_further10years_weeklysaves

end