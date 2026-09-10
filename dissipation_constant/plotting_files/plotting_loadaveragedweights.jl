function averaged_weights_models()

    # trying the averaging thing
    u3avg3 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/128_averagedweights_3dayoptimizations_3years_dailysaves/u.nc", "u");
    v3avg3 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/128_averagedweights_3dayoptimizations_3years_dailysaves/v.nc", "v");
    eta3avg3 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/128_averagedweights_3dayoptimizations_3years_dailysaves/eta.nc", "eta");

    u3avg7 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/128_averagedweights_3dayoptimizations_startfrom3years_7years_weeklysaves/u.nc", "u");
    v3avg7 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/128_averagedweights_3dayoptimizations_startfrom3years_7years_weeklysaves/v.nc", "v");
    eta3avg7 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/128_averagedweights_3dayoptimizations_startfrom3years_7years_weeklysaves/eta.nc", "eta");

    u3avgall = cat(u3avg3, u3avg7[:,:,2:end]; dims=3);
    v3avgall = cat(v3avg3, v3avg7[:,:,2:end]; dims=3);
    eta3avgall = cat(eta3avg3, eta3avg7[:,:,2:end]; dims=3);

    u10avg3 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/128_averagedweights_10dayoptimizations_3years_dailysaves/u.nc", "u");
    v10avg3 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/128_averagedweights_10dayoptimizations_3years_dailysaves/v.nc", "v");
    eta10avg3 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/128_averagedweights_10dayoptimizations_3years_dailysaves/eta.nc", "eta");

    u10avg7 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/128_averagedweights_10dayoptimizations_startfrom3year_7years_weeklysaves/u.nc", "u");
    v10avg7 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/128_averagedweights_10dayoptimizations_startfrom3year_7years_weeklysaves/v.nc", "v");
    eta10avg7 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/128_averagedweights_10dayoptimizations_startfrom3year_7years_weeklysaves/eta.nc", "eta");

    u10avgall = cat(u10avg3, u10avg7[:,:,2:end]; dims=3);
    v10avgall = cat(v10avg3, v10avg7[:,:,2:end]; dims=3);
    eta10avgall = cat(eta10avg3, eta10avg7[:,:,2:end]; dims=3);

    # inidividual weight integrations
    # 3 day
    u3single1 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday1_15years_weeklysaves/u.nc", "u");
    v3single1 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday1_15years_weeklysaves/v.nc", "v");
    eta3single1 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday1_15years_weeklysaves/eta.nc", "eta");

    u3single3 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday3_15years_weeklysaves/u.nc", "u");
    v3single3 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday3_15years_weeklysaves/v.nc", "v");
    eta3single3 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday3_15years_weeklysaves/eta.nc", "eta");

    u3single8 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday8_15years_weeklysaves/u.nc", "u");
    v3single8 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday8_15years_weeklysaves/v.nc", "v");
    eta3single8 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday8_15years_weeklysaves/eta.nc", "eta");

    u3single13 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday13_15years_weeklysaves/u.nc", "u");
    v3single13 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday13_15years_weeklysaves/v.nc", "v");
    eta3single13 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday13_15years_weeklysaves/eta.nc", "eta");

    u3single23 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday23_15years_weeklysaves/u.nc", "u");
    v3single23 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday23_15years_weeklysaves/v.nc", "v");
    eta3single23 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday23_15years_weeklysaves/eta.nc", "eta");

    u3single28 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday28_15years_weeklysaves/u.nc", "u");
    v3single28 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday28_15years_weeklysaves/v.nc", "v");
    eta3single28 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday28_15years_weeklysaves/eta.nc", "eta");

    u3single33 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday33_15years_weeklysaves/u.nc", "u");
    v3single33 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday33_15years_weeklysaves/v.nc", "v");
    eta3single33 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday33_15years_weeklysaves/eta.nc", "eta");

    u3single38 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday38_15years_weeklysaves/u.nc", "u");
    v3single38 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday38_15years_weeklysaves/v.nc", "v");
    eta3single38 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday38_15years_weeklysaves/eta.nc", "eta");

    u3single41 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday41_15years_weeklysaves/u.nc", "u");
    v3single41 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday41_15years_weeklysaves/v.nc", "v");
    eta3single41 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday41_15years_weeklysaves/eta.nc", "eta");

    u3single44 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday44_15years_weeklysaves/u.nc", "u");
    v3single44 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday44_15years_weeklysaves/v.nc", "v");
    eta3single44 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday44_15years_weeklysaves/eta.nc", "eta");

    u3single48 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday48_15years_weeklysaves/u.nc", "u");
    v3single48 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday48_15years_weeklysaves/v.nc", "v");
    eta3single48 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday48_15years_weeklysaves/eta.nc", "eta");

    u3single53 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday53_15years_weeklysaves/u.nc", "u");
    v3single53 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday53_15years_weeklysaves/v.nc", "v");
    eta3single53 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday53_15years_weeklysaves/eta.nc", "eta");

    u3single58 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday58_15years_weeklysaves/u.nc", "u");
    v3single58 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday58_15years_weeklysaves/v.nc", "v");
    eta3single58 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday58_15years_weeklysaves/eta.nc", "eta");

    u3single63 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday63_15years_weeklysaves/u.nc", "u");
    v3single63 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday63_15years_weeklysaves/v.nc", "v");
    eta3single63 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday63_15years_weeklysaves/eta.nc", "eta");

    u3single68 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday68_15years_weeklysaves/u.nc", "u");
    v3single68 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday68_15years_weeklysaves/v.nc", "v");
    eta3single68 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday68_15years_weeklysaves/eta.nc", "eta");

    u3single73 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday73_15years_weeklysaves/u.nc", "u");
    v3single73 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday73_15years_weeklysaves/v.nc", "v");
    eta3single73 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday73_15years_weeklysaves/eta.nc", "eta");

    u3single78 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday78_15years_weeklysaves/u.nc", "u");
    v3single78 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday78_15years_weeklysaves/v.nc", "v");
    eta3single78 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday78_15years_weeklysaves/eta.nc", "eta");

    u3single83 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday83_15years_weeklysaves/u.nc", "u");
    v3single83 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday83_15years_weeklysaves/v.nc", "v");
    eta3single83 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday83_15years_weeklysaves/eta.nc", "eta");

    u3single86 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday86_15years_weeklysaves/u.nc", "u");
    v3single86 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday86_15years_weeklysaves/v.nc", "v");
    eta3single86 = ncread("./dissipation_constant/manystates_singleinitcond_3dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/3day_single_initday86_15years_weeklysaves/eta.nc", "eta");

    # 10 day
    # day zero start
    u10single2 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday2_15years_weeklysaves/u.nc", "u");
    v10single2 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday2_15years_weeklysaves/v.nc", "v");
    eta10single2 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday2_15years_weeklysaves/eta.nc", "eta");

    u10single20 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday20_15years_weeklysaves/u.nc", "u");
    v10single20 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday20_15years_weeklysaves/v.nc", "v");
    eta10single20 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday20_15years_weeklysaves/eta.nc", "eta");

    u10single30 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday30_15years_weeklysaves/u.nc", "u");
    v10single30 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday30_15years_weeklysaves/v.nc", "v");
    eta10single30 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday30_15years_weeklysaves/eta.nc", "eta");

    u10single40 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday40_10years_weeklysaves/u.nc", "u");
    v10single40 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday40_10years_weeklysaves/v.nc", "v");
    eta10single40 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday40_10years_weeklysaves/eta.nc", "eta");
    
    u10single50 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday50_15years_weeklysaves/u.nc", "u");
    v10single50 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday50_15years_weeklysaves/v.nc", "v");
    eta10single50 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday50_15years_weeklysaves/eta.nc", "eta");

    u10single60 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday60_15years_weeklysaves/u.nc", "u");
    v10single60 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday60_15years_weeklysaves/v.nc", "v");
    eta10single60 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday60_15years_weeklysaves/eta.nc", "eta");
    
    u10single70 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday70_15years_weeklysaves/u.nc", "u");
    v10single70 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday70_15years_weeklysaves/v.nc", "v");
    eta10single70 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday70_15years_weeklysaves/eta.nc", "eta");

    u10single80 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday80_15years_weeklysaves/u.nc", "u");
    v10single80 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday80_15years_weeklysaves/v.nc", "v");
    eta10single80 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day0initcond/10day_single_initday80_15years_weeklysaves/eta.nc", "eta");

    # day 50 start
    u10single30start50 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day50initcond/10day_single_initday30_15years_weeklysaves/u.nc", "u");
    v10single30start50 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day50initcond/10day_single_initday30_15years_weeklysaves/v.nc", "v");
    eta10single30start50 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day50initcond/10day_single_initday30_15years_weeklysaves/eta.nc", "eta");

    u10single40start50 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day50initcond/10day_single_initday40_15years_weeklysaves/u.nc", "u");
    v10single40start50 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day50initcond/10day_single_initday40_15years_weeklysaves/v.nc", "v");
    eta10single40start50 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day50initcond/10day_single_initday40_15years_weeklysaves/eta.nc", "eta");

    u10single50start50 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day50initcond/10day_single_initday50_15years_weeklysaves/u.nc", "u");
    v10single50start50 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day50initcond/10day_single_initday50_15years_weeklysaves/v.nc", "v");
    eta10single50start50 = ncread("./dissipation_constant/manystates_singleinitcond_10dayoptimizations_allstartfrom20daysingle/individual_weight_integrations/day50initcond/10day_single_initday50_15years_weeklysaves/eta.nc", "eta");

end
