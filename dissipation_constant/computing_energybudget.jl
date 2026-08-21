"""
Like I did for the momentum budget I now want to compute the energy budget and see what it looks like.
Specifically I want to compute the energy transfer for each of the terms to see how they all contribute to
subgrid-scale forcing

I specifically want to look at each terms contribution to the subgrid-scale forcings.
To do this I need to compute each of the momentum budget components with (1) a high-resolution state and (2) the corresponding
coarse-grained high-resolution state. I then filter the result from (1) and look at the difference to see the contribution
to the subgrid-scale forcings. I'm going to start by doing this with just a single snapshot and may progress to computing an
averaged value.
"""
function compute_averagedtransfer()

    Plr = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        zb_forcing_momentum=false,
        zb_forcing_dissipation=false,
        zb_filtered=true,
        # nn_forcing_momentum=false,
        # nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=2,
        initial_cond="rest",
        # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
    );
    Slr = ShallowWaters.model_setup(Plr);

    Phr = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        zb_forcing_momentum=false,
        zb_forcing_dissipation=false,
        zb_filtered=true,
        # nn_forcing_momentum=false,
        # nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=1024,
        Ndays=2,
        initial_cond="rest",
        # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
    );
    Shr = ShallowWaters.model_setup(Phr);

    Plr2 = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=1e-12,
        H=500,
        cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        ω=1e-16,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        zb_forcing_momentum=false,
        zb_forcing_dissipation=false,
        zb_filtered=true,
        # nn_forcing_momentum=false,
        # nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=2,
        initial_cond="rest",
        # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
    );
    Slr2 = ShallowWaters.model_setup(Padveclr);

    Phr2 = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=1e-12,
        H=500,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        ω=1e-16,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        zb_forcing_momentum=false,
        zb_forcing_dissipation=false,
        zb_filtered=true,
        # nn_forcing_momentum=false,
        # nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=1024,
        Ndays=2,
        initial_cond="rest",
        # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
    );
    Shr2 = ShallowWaters.model_setup(Padvechr);

    tlr = 225 * Slr.grid.dtint
    thr = 1800 * Shr.grid.dtint

    tendency_transfer = zeros(65);
    advection_transfer = zeros(65);
    pressure_transfer = zeros(65);
    viscosity_transfer = zeros(65);
    coriolis_transfer = zeros(65);
    bottomdrag_transfer = zeros(65);

    for n = 1:1096

        Seulerlr = deepcopy(Slr);
        Seulerhr = deepcopy(Shr);
        Sadveclr = deepcopy(Slr2);
        Sadvechr = deepcopy(Shr2)

        ulr_, vlr_, etalr_ = ShallowWaters.add_halo(uhrcgall[:,:,n], vhrcgall[:,:,n], etahrcgall[:,:,n], Slr)
        uhr_, vhr_, etahr_ = ShallowWaters.add_halo(uhrall[:,:,n], vhrall[:,:,n], etahrall[:,:,n], Shr)

        Seulerlr.Prog.u = ulr_
        Seulerlr.Prog.v = vlr_
        Seulerlr.Prog.η = etalr_

        Seulerhr.Prog.u = uhr_
        Seulerhr.Prog.v = vhr_
        Seulerhr.Prog.η = etahr_

        Sadveclr.Prog.u = ulr_
        Sadveclr.Prog.v = vlr_
        Sadveclr.Prog.η = etalr_

        Sadvechr.Prog.u = uhr_
        Sadvechr.Prog.v = vhr_
        Sadvechr.Prog.η = etahr_

        dulr, dvlr, adv_ulr, adv_vlr, fulr, fvlr, detagdxlr, detagdylr, bottomdragulr, bottomdragvlr, visculr, viscvlr, Fxlr = compute_mom_budget(Seulerlr, Sadveclr, tlr, n);
        duhr, dvhr, adv_uhr, adv_vhr, fuhr, fvhr, detagdxhr, detagdyhr, bottomdraguhr, bottomdragvhr, viscuhr, viscvhr, Fxhr = compute_mom_budget(Seulerhr, Sadvechr, thr, n);

        duhrcg, dvhrcg = filter_hr(duhr, dvhr);
        dutransfer, dvtransfer = compute_transfer(uhrcgall[:,:,n], vhrcgall[:,:,n], 2 * duhrcg ./ Shr.grid.dtint .- 2 * dulr ./ Slr.grid.dtint, 2 * dvhrcg ./ Shr.grid.dtint .- 2 * dvlr ./ Slr.grid.dtint);

        tendency_transfer += (dutransfer + dvtransfer)

        advu_hrcg, advv_hrcg = filter_hr(adv_uhr, adv_vhr);
        advu_transfer, advv_transfer = compute_transfer(uhrcgall[:,:,n], vhrcgall[:,:,n], adv_ulr .- advu_hrcg, adv_vlr .- advv_hrcg);

        advection_transfer += advu_transfer + advv_transfer

        fvhrcg, fuhrcg = filter_hr(fvhr, fuhr);
        coriolisu_transfer, coriolisv_transfer = compute_transfer(uhrcgall[:,:,n], vhrcgall[:,:,n], -fvlr .+ fvhrcg, fulr .- fuhrcg);

        coriolis_transfer += coriolisu_transfer + coriolisv_transfer

        detagdxhrcg, detagdyhrcg = filter_hr(detagdxhr, detagdyhr);
        pressureu_transfer, pressurev_transfer = compute_transfer(uhrcgall[:,:,n], vhrcgall[:,:,n], detagdxlr .- detagdxhrcg, detagdylr .- detagdyhrcg);

        pressure_transfer += pressureu_transfer + pressurev_transfer

        bottomdraguhrcg, bottomdragvhrcg = filter_hr(bottomdraguhr, bottomdragvhr);
        bottomdragu_transfer, bottomdragv_transfer = compute_transfer(uhrcgall[:,:,n], vhrcgall[:,:,n], bottomdraguhrcg .- bottomdragulr, bottomdragvhrcg .- bottomdragvlr);

        bottomdrag_transfer += bottomdragu_transfer + bottomdragv_transfer

        viscuhrcg, viscvhrcg = filter_hr(viscuhr, viscvhr);
        viscu_transfer, viscv_transfer = compute_transfer(uhrcgall[:,:,n], vhrcgall[:,:,n],viscuhrcg .- visculr, viscvhrcg .- viscvlr);

        viscosity_transfer += viscu_transfer + viscv_transfer

    end

    return tendency_transfer, advection_transfer, pressure_transfer, viscosity_transfer, coriolis_transfer, bottomdrag_transfer

end

function energy_budget_offline()

    Plr = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        zb_forcing_momentum=false,
        zb_forcing_dissipation=false,
        zb_filtered=true,
        # nn_forcing_momentum=false,
        # nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=2,
        initial_cond="rest",
        # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
    );
    Slr = ShallowWaters.model_setup(Plr);

    Phr = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        zb_forcing_momentum=false,
        zb_forcing_dissipation=false,
        zb_filtered=true,
        # nn_forcing_momentum=false,
        # nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=1024,
        Ndays=2,
        initial_cond="rest",
        # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
    );
    Shr = ShallowWaters.model_setup(Phr);

    Padveclr = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=1e-12,
        H=500,
        cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        ω=1e-16,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        zb_forcing_momentum=false,
        zb_forcing_dissipation=false,
        zb_filtered=true,
        # nn_forcing_momentum=false,
        # nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=2,
        initial_cond="rest",
        # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
    );
    Sadveclr = ShallowWaters.model_setup(Padveclr);

    Padvechr = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=1e-12,
        H=500,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        ω=1e-16,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        zb_forcing_momentum=false,
        zb_forcing_dissipation=false,
        zb_filtered=true,
        # nn_forcing_momentum=false,
        # nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=1024,
        Ndays=2,
        initial_cond="rest",
        # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
    );
    Sadvechr = ShallowWaters.model_setup(Padvechr);

    Seulerlr = deepcopy(Slr);
    Seulerhr = deepcopy(Shr);

    tlr = 225 * Slr.grid.dtint
    thr = 1800 * Shr.grid.dtint

    n = 1096

    ulr_, vlr_, etalr_ = ShallowWaters.add_halo(uhrcgall[:,:,n], vhrcgall[:,:,n], etahrcgall[:,:,n], Slr)
    uhr_, vhr_, etahr_ = ShallowWaters.add_halo(uhrall[:,:,n], vhrall[:,:,n], etahrall[:,:,n], Shr)

    Seulerlr.Prog.u = ulr_
    Seulerlr.Prog.v = vlr_
    Seulerlr.Prog.η = etalr_

    Seulerhr.Prog.u = uhr_
    Seulerhr.Prog.v = vhr_
    Seulerhr.Prog.η = etahr_

    Sadveclr.Prog.u = ulr_
    Sadveclr.Prog.v = vlr_
    Sadveclr.Prog.η = etalr_

    Sadvechr.Prog.u = uhr_
    Sadvechr.Prog.v = vhr_
    Sadvechr.Prog.η = etahr_

    dulr, dvlr, adv_ulr, adv_vlr, fulr, fvlr, detagdxlr, detagdylr, bottomdragulr, bottomdragvlr, visculr, viscvlr, Fxlr = compute_mom_budget(Seulerlr, Sadveclr, tlr, n);
    duhr, dvhr, adv_uhr, adv_vhr, fuhr, fvhr, detagdxhr, detagdyhr, bottomdraguhr, bottomdragvhr, viscuhr, viscvhr, Fxhr = compute_mom_budget(Seulerhr, Sadvechr, thr, n);

    duhrcg, dvhrcg = filter_hr(duhr, dvhr);
    dutransfer, dvtransfer = compute_transfer(uhrcgall[:,:,n], vhrcgall[:,:,n], 2 * duhrcg ./ Shr.grid.dtint .- 2 * dulr ./ Slr.grid.dtint, 2 * dvhrcg ./ Shr.grid.dtint .- 2 * dvlr ./ Slr.grid.dtint);

    advu_hrcg, advv_hrcg = filter_hr(adv_uhr, adv_vhr);
    advu_transfer, advv_transfer = compute_transfer(uhrcgall[:,:,n], vhrcgall[:,:,n], adv_ulr .- advu_hrcg, adv_vlr .- advv_hrcg);

    fvhrcg, fuhrcg = filter_hr(fvhr, fuhr);
    coriolisu_transfer, coriolisv_transfer = compute_transfer(uhrcgall[:,:,n], vhrcgall[:,:,n], -fvlr .+ fvhrcg, fulr .- fuhrcg);

    detagdxhrcg, detagdyhrcg = filter_hr(detagdxhr, detagdyhr);
    pressureu_transfer, pressurev_transfer = compute_transfer(uhrcgall[:,:,n], vhrcgall[:,:,n], detagdxlr .- detagdxhrcg, detagdylr .- detagdyhrcg);

    bottomdraguhrcg, bottomdragvhrcg = filter_hr(bottomdraguhr, bottomdragvhr);
    bottomdragu_transfer, bottomdragv_transfer = compute_transfer(uhrcgall[:,:,n], vhrcgall[:,:,n], bottomdraguhrcg .- bottomdragulr, bottomdragvhrcg .- bottomdragvlr);

    viscuhrcg, viscvhrcg = filter_hr(viscuhr, viscvhr);
    viscu_transfer, viscv_transfer = compute_transfer(uhrcgall[:,:,n], vhrcgall[:,:,n],viscuhrcg .- visculr, viscvhrcg .- viscvlr);

    # detagdx2 = zeros(127,128)
    # detagdy2 = zeros(128,127)

    # ShallowWaters.∂x!(detagdx2, (Seuler.constants.g .* etahrcgall[:,:,n])./ Seuler.grid.Δ)
    # ShallowWaters.∂y!(detagdy2, (Seuler.constants.g .* etahrcgall[:,:,n])./ Seuler.grid.Δ)

    # check that budget was computed correctly

    tendencyhr = ((2 .* duhr) ./ Shr.grid.dtint)
    budgetsumhr = (adv_uhr + fvhr - detagdxhr + Fxhr + viscuhr + bottomdraguhr)

    fig = Figure(size=(1000, 300), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 1024),
    LinRange(0, 3840, 1024),
    tendencyhr,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Tendency"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,2], hm1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 1024),
    LinRange(0, 3840, 1024),
    budgetsumhr,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Budget sum"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,4], hm2)

    ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 1024),
    LinRange(0, 3840, 1024),
    tendencyhr .- budgetsumhr,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="HR tendency - budget sum"),
    colorrange=(-1.5e-15, 1.5e-15)
    );
    Colorbar(fig[1,6], hm3)

    # subgrid scale forcing contributions

    hrcgtendency = 2 * duhrcg ./ Shr.grid.dtint
    lrtendency = 2 * dulr ./ Slr.grid.dtint

    sgs_advection = -(adv_ulr .- advu_hrcg)
    sgs_coriolis = -fvlr .+ fvhrcg
    sgs_pressure = detagdxlr .- detagdxhrcg
    sgs_viscosity = viscuhrcg .- visculr
    sgs_bottomdrag = bottomdraguhrcg .- bottomdragulr

    Fxhrcg, _ = filter_hr(Fxhr, viscvhr);
    sgs_wind = Fxlr - Fxhrcg

    sgs_sum = sgs_advection + sgs_coriolis + sgs_pressure + sgs_viscosity + sgs_bottomdrag - sgs_wind

    fig = Figure(size=(1050, 700), fontsize=15);

    Label(
        fig[0, 3],
        "SGS contributions in coarse-grained high-resolution model, u equation",
        fontsize = 25,
        tellwidth = false
    )

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    hrcgtendency .- lrtendency,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Euler tendency"),
    colorrange=(-2e-5, 2e-5)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    ax12, hm12 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sgs_sum,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Sum of SGS terms"),
    colorrange=(-2e-5, 2e-5)
    );
    Colorbar(fig[1,4], hm12)
    hidexdecorations!(ax12)

    ax13, hm13 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (hrcgtendency .- lrtendency) .- sgs_advection,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Tendency - Sum"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,6], hm13)
    hidexdecorations!(ax13)

    ax2, hm2 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sgs_advection,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Nonlinear advection"),
    colorrange=(-2e-5, 2e-5)
    );
    hidedecorations!(ax2)
    Colorbar(fig[2,2], hm2)

    ax3, hm3 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sgs_coriolis,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Coriolis"),
    colorrange=(-2e-5, 2e-5)
    );
    hidedecorations!(ax3)
    Colorbar(fig[2,4], hm3)

    ax4, hm4 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sgs_pressure,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Pressure gradient"),
    colorrange=(-1e-5, 1e-5)
    );
    Colorbar(fig[2,6], hm4)
    hidexdecorations!(ax4)

    ax5, hm5 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sgs_viscosity,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Viscosity"),
    colorrange=(-2e-7, 2e-7)
    );
    hidedecorations!(ax5)
    Colorbar(fig[3,2], hm5)

    ax6, hm6 = heatmap(fig[3,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sgs_bottomdrag,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Bottom drag"),
    colorrange=(-1e-8, 1e-8)
    );
    hidedecorations!(ax6)
    Colorbar(fig[3,4], hm6)

    lr_freq = 1/30 .* freq(periodogram(uhrcgall[:,:,10]; radialavg=true, radialsum=false));

    transfersum = lr_freq.* ( -(advu_transfer + advv_transfer) +
        (coriolisu_transfer + coriolisv_transfer) +
        (pressureu_transfer + pressurev_transfer) + 
        (viscu_transfer + viscv_transfer) +
        (bottomdragu_transfer + bottomdragv_transfer)
    )

    fig = Figure(size=(1050, 700), fontsize=15);

    ax = Axis(fig[1,1],
        xscale=log10,
        xlabel="Wavelength (km)",
        ylabel=L"k T(k)",
        title="Kinetic energy transfer",
        xreversed=true,
        xticks=[700, 100, 30, 10, 2]
    )
    lines!(ax, 1 ./ lr_freq[2:end], (lr_freq.*(dutransfer + dvtransfer))[2:end], label="Tendency, Euler", color=:black)
    lines!(ax, 1 ./ lr_freq[2:end], (-lr_freq.*(advu_transfer + advv_transfer))[2:end], label="Advection")
    lines!(ax, 1 ./ lr_freq[2:end], (lr_freq.*(coriolisu_transfer + coriolisv_transfer))[2:end], label="Coriolis")
    lines!(ax, 1 ./ lr_freq[2:end], (lr_freq.*(pressureu_transfer + pressurev_transfer))[2:end], label="Pressure gradient")
    lines!(ax, 1 ./ lr_freq[2:end], (lr_freq.*(viscu_transfer + viscv_transfer))[2:end], label="Viscosity")
    lines!(ax, 1 ./ lr_freq[2:end], (lr_freq.*(bottomdragu_transfer + bottomdragv_transfer))[2:end], label="Bottom drag")
    lines!(ax, 1 ./ lr_freq[2:end], transfersum[2:end], label="Sum of terms", linestyle=:dash)

    axislegend(ax, position=:rt)

    # only re run if necessary, this computes the sum of 3 years of transfers to then look at the average. Still need to divide by 1096 when plotting
    # tendency_transfer, advection_transfer, pressure_transfer, viscosity_transfer, coriolis_transfer, bottomdrag_transfer = compute_averagedtransfer()

    transfers = load_object("./dissipation_constant/computing_trueS/energybudget_transfers_threeyearsums_computedwhrcgstates_tend_advec_pressure_coriolis_visc_bottomdrag.jld2");

    tendency_transfer = transfers[1];
    advection_transfer = transfers[2];
    pressure_transfer = transfers[3];
    coriolis_transfer = transfers[4];
    viscosity_transfer = transfers[5];
    bottomdrag_transfer = transfers[6];

    transfersum = -advection_transfer + pressure_transfer + viscosity_transfer + coriolis_transfer + bottomdrag_transfer

    fig = Figure(size=(1050, 700), fontsize=15);
    colors=Makie.resample_cmap(:lighttest, 9)
    ax = Axis(fig[1,1],
        xscale=log10,
        xlabel="Wavelength (km)",
        ylabel=L"k T(k)",
        title="Kinetic energy transfer",
        xreversed=true,
        xticks=[700, 100, 30, 10, 2]
    )
    lines!(ax, 1 ./ lr_freq[2:end], (lr_freq.*(tendency_transfer)./1096)[2:end], label="Tendency, Euler", color=:aqua)
    lines!(ax, 1 ./ lr_freq[2:end], (-lr_freq.*(advection_transfer)./1096)[2:end], label="Advection",color=colors[1])
    lines!(ax, 1 ./ lr_freq[2:end], (lr_freq.*(coriolis_transfer)./1096)[2:end], label="Coriolis",color=colors[2])
    lines!(ax, 1 ./ lr_freq[2:end], (lr_freq.*(pressure_transfer)./1096)[2:end], label="Pressure gradient",color=colors[4])
    lines!(ax, 1 ./ lr_freq[2:end], (lr_freq.*(viscosity_transfer)./1096)[2:end], label="Viscosity",color=colors[5])
    lines!(ax, 1 ./ lr_freq[2:end], (lr_freq.*(bottomdrag_transfer)./1096)[2:end], label="Bottom drag",color=colors[8])
    lines!(ax, 1 ./ lr_freq[2:end], (lr_freq .* transfersum ./ 1096)[2:end], label="Sum of terms", linestyle=:dash)

    axislegend(ax, position=:rt)

end