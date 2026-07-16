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

    # 30 km relative vorticity over the first three years
    ζzb = zeros(129,129,1096);
    ζnoparam = zeros(129,129,1096);
    ζmulti3more = zeros(129,129,1096);
    ζmulti2 = zeros(129,129,1096);
    ζ20s = zeros(129,129,1096);
    ζ30s = zeros(129,129,1096);
    ζ10s = zeros(129,129,1096);
    ζ5s = zeros(129,129,1096)
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
    five = kde(vec(ζ5s));
    ten = kde(vec(ζ10s));
    twenty = kde(vec(ζ20s));
    thirty = kde(vec(ζ30s));

    fig = Figure(size = (650, 420));
    ax = Axis(
        fig[1, 1],
        xlabel = "Relative Vorticity",
        ylabel = "Probability Density",
        title = "Three-year relative vorticity probability density"
        # yscale = log10
    )
    lines!(ax, hr.x, hr.density, label="3.75 km", color=colors[1])
    lines!(ax, hrcg.x, hrcg.density, label="Filtered, coarse-grained 3.75 km",color=:black)
    lines!(ax, zb.x, zb.density, label="ZB20", color=:red)
    lines!(ax, noparam.x, noparam.density, label="No closure, 30 km",color=:gray)
    # lines!(ax, multi2.x, multi2.density, label="Ensemble 2 day", color=colors[2])
    # lines!(ax, multi3more.x, multi3more.density, label="Ensemble 3 day", color=colors[3])
    # lines!(ax, multi10.x, multi10.density, label="Ensemble 10 day", color=colors[4])
    lines!(ax, ten.x, ten.density, label="10 day", color=colors[2])#, linestyle=:dashdot)
    lines!(ax, twenty.x, twenty.density, label="20 day", color=colors[3])#, linestyle=:dashdot)
    lines!(ax, thirty.x, thirty.density, label="30 day", color=colors[4])#, linestyle=:dashdot)
    lines!(ax, multi20.x, multi20.density, label="Ensemble 20 day", color=colors[6])#,linestyle=:dash)

    # lines!(ax, multi20.x, multi20.density, label="Ensemble 20 day")#,linestyle=:dash)
    # lines!(ax, thirty.x, thirty.density, label="30 day")#, linestyle=:dashdot)

    Legend(fig[1,2], ax)
    xlims!(ax, -0.000007, 0.000007)
    # ylims!(ax, 1, 10^4)

    # 10 year figure

    hr10 = kde(vec(cat(ζhr[:,:,1:7:1096], ζhr[:,:,1097:end]; dims=3)); npoints=16384);
    hrcg10 = kde(vec(cat(ζhrcg[:,:,1:7:1096], ζhrcg[:,:,1097:end]; dims=3)));
    zb10 = kde(vec(ζzb10));
    ten10 = kde(vec(ζ10s10));
    twenty10 = kde(vec(ζ20s10));
    thirty10 = kde(vec(ζ30s10));
    noparam10 = kde(vec(ζnoparam10));
    multi210 = kde(vec(ζmulti210));
    multi3more10 = kde(vec(ζmulti3more10));
    multi2010 = kde(vec(ζmulti2010));
    multi1010 = kde(vec(ζmulti1010));

    fig = Figure(size = (650, 420));
    ax = Axis(
        fig[1, 1],
        xlabel = "Relative Vorticity",
        ylabel = "Probability Density",
        title = "Ten-year relative vorticity probability density"
        # yscale = log10
    )
    # lines!(ax, hr10.x, hr10.density, label="3.75 km", color=colors[1])
    lines!(ax, hrcg10.x, hrcg10.density, label="Filtered, coarse-grained 3.75 km",color=:black)
    lines!(ax, zb10.x, zb10.density, label="ZB20", color=:red)
    lines!(ax, noparam10.x, noparam10.density, label="No closure, 30 km",color=:gray)
    lines!(ax, multi210.x, multi210.density, label="Ensemble 2 day", color=colors[1])
    lines!(ax, multi3more10.x, multi3more10.density, label="Ensemble 3 day", color=colors[2])
    lines!(ax, multi1010.x, multi1010.density, label="Ensemble 10 day", color=colors[3])
    # lines!(ax, ten10.x, ten10.density, label="10 day", color=colors[2])#, linestyle=:dashdot)
    # lines!(ax, twenty10.x, twenty10.density, label="20 day", color=colors[3])#, linestyle=:dashdot)
    # lines!(ax, thirty10.x, thirty10.density, label="30 day", color=colors[4])#, linestyle=:dashdot)
    # lines!(ax, multi2010.x, multi2010.density, label="Ensemble 20 day", color=colors[6])#,linestyle=:dash)

    Legend(fig[1,2], ax)
    xlims!(ax, -0.000007, 0.000007)
    # ylims!(ax, 1, 10^4)

    # 3 and 10 year

    fig = Figure(size = (1100, 500));
    ax = Axis(
        fig[1, 1],
        xlabel = "Relative Vorticity",
        ylabel = "Probability Density",
        title = "Three-year relative vorticity probability density"
        # yscale = log10
    )

    ax2 = Axis(
        fig[1, 2],
        xlabel = "Relative Vorticity",
        # ylabel = "Probability Density",
        title = "Ten-year relative vorticity probability density"
        # yscale = log10
    )

    lines!(ax, hr.x, hr.density, label="3.75 km", color=colors[1])
    lines!(ax, hrcg.x, hrcg.density, label="Filtered, coarse-grained 3.75 km",color=:black)
    lines!(ax, zb.x, zb.density, label="ZB20", color=:red)
    lines!(ax, noparam.x, noparam.density, label="No closure, 30 km",color=:gray)
    # lines!(ax, multi2.x, multi2.density, label="Ensemble 2 day", color=colors[2])
    # lines!(ax, multi3more.x, multi3more.density, label="Ensemble 3 day", color=colors[3])
    # lines!(ax, multi10.x, multi10.density, label="Ensemble 10 day", color=colors[4])
    lines!(ax, ten.x, ten.density, label="10 day", color=colors[2])#, linestyle=:dashdot)
    lines!(ax, twenty.x, twenty.density, label="20 day", color=colors[3])#, linestyle=:dashdot)
    lines!(ax, thirty.x, thirty.density, label="30 day", color=colors[4])#, linestyle=:dashdot)
    lines!(ax, multi20.x, multi20.density, label="Ensemble 20 day", color=colors[6])#,linestyle=:dash)

    xlims!(ax, -0.000007, 0.000007)
    ylims!(ax, 0, 2e5)

    lines!(ax2, hr10.x, hr10.density, label="3.75 km", color=colors[1])
    lines!(ax2, hrcg10.x, hrcg10.density, label="Filtered, coarse-grained 3.75 km",color=:black)
    lines!(ax2, zb10.x, zb10.density, label="ZB20", color=:red)
    lines!(ax2, noparam10.x, noparam10.density, label="No closure, 30 km",color=:gray)
    # lines!(ax, multi210.x, multi210.density, label="Ensemble 2 day", color=colors[2])
    # lines!(ax, multi3more10.x, multi3more10.density, label="Ensemble 3 day", color=colors[3])
    # lines!(ax, multi1010.x, multi1010.density, label="Ensemble 10 day", color=colors[4])
    lines!(ax2, ten10.x, ten10.density, label="10 day", color=colors[2])#, linestyle=:dashdot)
    lines!(ax2, twenty10.x, twenty10.density, label="20 day", color=colors[3])#, linestyle=:dashdot)
    lines!(ax2, thirty10.x, thirty10.density, label="30 day", color=colors[4])#, linestyle=:dashdot)
    lines!(ax2, multi2010.x, multi2010.density, label="Ensemble 20 day", color=colors[6])#,linestyle=:dash)

    xlims!(ax2, -0.000007, 0.000007)
    ylims!(ax2, 0, 2e5)

    Legend(fig[2, 1:2], ax, orientation = :horizontal)

    # plot of the vorticity itself

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

    # plot of the vorticity itself, best template
    # stable results in first three years
    fig = Figure(size=(900, 520), fontsize=15);

    Label(
        fig[0, 2],
        L"\zeta(3 \text{ years}, x, y)",
        fontsize = 20,
        tellwidth = false
    )

    t = 1096
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 1024),
    LinRange(0, 3840, 1024),
    ζhrcg[:,:, t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Filtered, \n coarse-grained 3.75 km"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[1,2], hm1, label="1/s")
    hidexdecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζnoparam[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="No closure, 30 km"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[1,4], hm1, label="1/s")
    hidedecorations!(ax2)

    ax3, hm3 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζzb[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="ZB20"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[1,6], hm1, label="1/s")
    hidedecorations!(ax3)

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti2[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 2 day"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[2,2], hm1, label="1/s")

    ax5, hm5 = heatmap(fig[2,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti3more[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 3 day"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[2,4], hm1, label="1/s")
    hideydecorations!(ax5)

    ax6, hm6 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti10[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 10 day"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[2,6], hm1, label="1/s")
    hideydecorations!(ax6)

    Colorbar(fig[1:2,4], hm1, label="1/s")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 2] = GridLayout()
    gc = fig[1, 3] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 2] = GridLayout()
    gf = fig[2, 3] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # examining the vorticity of the results that diverge over the last seven years
    # results in first 3 years
    fig = Figure(size=(900, 520), fontsize=15);

    Label(
        fig[0, 2],
        L"\zeta(3 \text{ years}, x, y)",
        fontsize = 20,
        tellwidth = false
    )

    t = 1096
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 1024),
    LinRange(0, 3840, 1024),
    ζhrcg[:,:, t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Filtered, \n coarse-grained 3.75 km"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[1,2], hm1, label="1/s")

    ax2, hm2 = heatmap(fig[1,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζzb[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="ZB20"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[1,4], hm1, label="1/s")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζ10s[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="10 day"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[1,6], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζ20s[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="20 day"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[2,2], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζ30s[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="30 day"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[2,4], hm1, label="1/s")

    ax3, hm3 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti20[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 20 day"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[2,6], hm1, label="1/s")

    Colorbar(fig[1:2,4], hm1, label="1/s")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 2] = GridLayout()
    gc = fig[1, 3] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 2] = GridLayout()
    gf = fig[2, 3] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # diverging results in final 7 years
    fig = Figure(size=(900, 520), fontsize=15);

    Label(
        fig[0, 2],
        L"\zeta(10 \text{ years}, x, y)",
        fontsize = 20,
        tellwidth = false
    )

    t = 522
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 1024),
    LinRange(0, 3840, 1024),
    ζhrcg[:,:,end],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Filtered, \n coarse-grained 3.75 km"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[1,2], hm1, label="1/s")

    ax2, hm2 = heatmap(fig[1,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζzb10[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="ZB20"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[1,4], hm1, label="1/s")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζ10s10[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="10 day"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[1,6], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζ20s10[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="20 day"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[2,2], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζ30s10[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="30 day"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[2,4], hm1, label="1/s")

    ax3, hm3 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ζmulti2010[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 20 day"),
    colorrange=(-maximum(abs.(ζhrcg[:,:, t])),maximum(abs.(ζhrcg[:,:, t])))
    );
    # Colorbar(fig[2,6], hm1, label="1/s")

    Colorbar(fig[1:2,4], hm1, label="1/s")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 2] = GridLayout()
    gc = fig[1, 3] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 2] = GridLayout()
    gf = fig[2, 3] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

end
