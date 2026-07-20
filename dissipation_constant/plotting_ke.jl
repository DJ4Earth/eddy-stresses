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

    # appendix

    kespec = zeros(Float64, N)
    hybrid = zeros(Float64, N)
    fourier = zeros(Float64, N)
    kespecpd = zeros(Float64, N)

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
        # kespec[j] = sum(abs2, ukespec[:,:,j]) + sum(abs2, vkespec[:,:,j])
        # hybrid[j] = sum(abs2, uhybrid[:,:,j]) + sum(abs2, vhybrid[:,:,j])
        # fourier[j] = sum(abs2, ufourier[:,:,j]) + sum(abs2, vfourier[:,:,j])
        # kespecpd[j] = sum(abs2, ukespecpd[:,:,j]) + sum(abs2, vkespecpd[:,:,j])
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
    # multi510 = zeros(Float64, N)
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

    end

    hrcg = zeros(Float64, 1461)
    threedayavg = zeros(Float64, 1461)
    tendayavg = zeros(Float64, 1461)

    for j = 1:1461
            hrcg[j] = sum(abs2, uhrcgall[:,:,j]) + sum(abs2, vhrcgall[:,:,j])
            # threedayavg[j] = sum(abs2, u3avgall[:,:,j]) + sum(abs2, v3avgall[:,:,j])
            # tendayavg[j] = sum(abs2, u10avgall[:,:,j]) + sum(abs2, v10avgall[:,:,j])
    end

    multi210p = zeros(Float64, 522)
    multi310p = zeros(Float64, 522)
    multi1010p = zeros(Float64, 522)

    avg3single = zeros(Float64, 522)
    avg10single = zeros(Float64, 424)
    for j = 1:522
            multi210p[j] = sum(abs2, umulti210more[:,:,j]) + sum(abs2, vmulti210more[:,:,j])
            multi310p[j] = sum(abs2, umulti310more[:,:,j]) + sum(abs2, vmulti310more[:,:,j])
            multi1010p[j] = sum(abs2, umulti1010more[:,:,j]) + sum(abs2, vmulti1010more[:,:,j])
            # avg3single[j] = sum(abs2, u3single[:,:,j]) + sum(abs2, v3single[:,:,j])
    end

    for j = 1:424
            avg10single[j] = sum(abs2, u10single[:,:,j]) + sum(abs2, v10single[:,:,j])
    end

    single31 = zeros(Float64, 783)
    single33 = zeros(Float64, 783)
    single38 = zeros(Float64, 783)
    single313 = zeros(Float64, 783)
    single323 = zeros(Float64, 783)
    single328 = zeros(Float64, 783)
    single333 = zeros(Float64, 783)
    single338 = zeros(Float64, 783)
    single341 = zeros(Float64, 783)
    single344 = zeros(Float64, 783)
    single348 = zeros(Float64, 783)
    single353 = zeros(Float64, 783)
    single358 = zeros(Float64, 783)
    single363 = zeros(Float64, 783)
    single368 = zeros(Float64, 783)
    single373 = zeros(Float64, 783)
    single378 = zeros(Float64, 783)
    single383 = zeros(Float64, 783)
    single386 = zeros(Float64, 783)

    single102 = zeros(Float64, 783)
    single1020 = zeros(Float64, 783)
    single1030 = zeros(Float64, 783)
    single1050 = zeros(Float64, 783)
    single1060 = zeros(Float64, 783)
    single1070 = zeros(Float64, 783)
    single1080 = zeros(Float64, 783)

    single1030_50 = zeros(Float64, 783)
    single1040_50 = zeros(Float64, 783)
    single1050_50 = zeros(Float64, 783)

    for j = 1:783

        single31[j] = sum(abs2, u3single1[:,:,j]) + sum(abs2, v3single1[:,:,j])
        single33[j] = sum(abs2, u3single3[:,:,j]) + sum(abs2, v3single3[:,:,j])
        single38[j] = sum(abs2, u3single8[:,:,j]) + sum(abs2, v3single8[:,:,j])
        single313[j] = sum(abs2, u3single13[:,:,j]) + sum(abs2, v3single13[:,:,j])
        single323[j] = sum(abs2, u3single23[:,:,j]) + sum(abs2, v3single23[:,:,j])
        single328[j] = sum(abs2, u3single28[:,:,j]) + sum(abs2, v3single8[:,:,j])
        single333[j] = sum(abs2, u3single33[:,:,j]) + sum(abs2, v3single33[:,:,j])
        single338[j] = sum(abs2, u3single38[:,:,j]) + sum(abs2, v3single38[:,:,j])
        single341[j] = sum(abs2, u3single41[:,:,j]) + sum(abs2, v3single41[:,:,j])
        single344[j] = sum(abs2, u3single44[:,:,j]) + sum(abs2, v3single44[:,:,j])
        single348[j] = sum(abs2, u3single48[:,:,j]) + sum(abs2, v3single48[:,:,j])
        single353[j] = sum(abs2, u3single53[:,:,j]) + sum(abs2, v3single53[:,:,j])
        single358[j] = sum(abs2, u3single58[:,:,j]) + sum(abs2, v3single58[:,:,j])
        single363[j] = sum(abs2, u3single63[:,:,j]) + sum(abs2, v3single63[:,:,j])
        single368[j] = sum(abs2, u3single68[:,:,j]) + sum(abs2, v3single68[:,:,j])
        single373[j] = sum(abs2, u3single73[:,:,j]) + sum(abs2, v3single73[:,:,j])
        single378[j] = sum(abs2, u3single78[:,:,j]) + sum(abs2, v3single78[:,:,j])
        single383[j] = sum(abs2, u3single83[:,:,j]) + sum(abs2, v3single83[:,:,j])
        single386[j] = sum(abs2, u3single86[:,:,j]) + sum(abs2, v3single86[:,:,j])
    
        single102[j] = sum(abs2, u10single2[:,:,j]) + sum(abs2, v10single2[:,:,j])
        single1020[j] = sum(abs2, u10single20[:,:,j]) + sum(abs2, v10single20[:,:,j])
        single1030[j] = sum(abs2, u10single30[:,:,j]) + sum(abs2, v10single30[:,:,j])
        single1050[j] = sum(abs2, u10single50[:,:,j]) + sum(abs2, v10single50[:,:,j])
        single1060[j] = sum(abs2, u10single60[:,:,j]) + sum(abs2, v10single60[:,:,j])
        single1070[j] = sum(abs2, u10single70[:,:,j]) + sum(abs2, v10single70[:,:,j])
        single1080[j] = sum(abs2, u10single80[:,:,j]) + sum(abs2, v10single80[:,:,j])

        single1030_50[j] = sum(abs2, u10single30start50[:,:,j]) + sum(abs2, v10single30start50[:,:,j])
        single1040_50[j] = sum(abs2, u10single40start50[:,:,j]) + sum(abs2, v10single40start50[:,:,j])
        single1050_50[j] = sum(abs2, u10single50start50[:,:,j]) + sum(abs2, v10single50start50[:,:,j])
    end

    multi2all = cat(multi210, multi210p[2:end]; dims=1);
    multi3all = cat(multi3more10, multi310p[2:end]; dims=1);
    multi10all = cat(multi1010, multi1010p[2:end]; dims=1);

    # 3 year, 10 year, 20 year

    colors = Makie.wong_colors()

    fig = Figure(size=(1000, 600), fontsize=15);

    hrcg_10 = cat(hrcg[1:7:1096], hrcg[1097:end];dims=1)

    ax3 = Axis(fig[1,1],
        # xlabel="Day",
        ylabel="Energy",
        title="Spatially averaged KE over 10 years"
    )
    lines!(ax3, LinRange(0, 10, 522),  hrcg_10 ./ (128*127), label="Filtered, coarse-grained 3.75 km",color=:black)
    lines!(ax3, LinRange(0, 10, 522), noparam10 ./ (128*127), label="No closure, 30 km",color=:gray)
    lines!(ax3, LinRange(0, 10, 522),  zb10 ./ (128*127), label="ZB20",color=:red)
    # lines!(fig[1,1], LinRange(0, 10*365, 522), fiveday./ (128*127), label="Online closure, 5 day")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), tenday10./ (128*127), label="Online closure, 10 day")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), twentyday./ (128*127), label="Online closure, 20 day")
    # lines!(fig[1,1], LinRange(0, 10*365, 522), thirtyday./ (128*127), label="Online closure, 30 day")
    lines!(ax3, LinRange(0, 10, 522), multi210 ./ (128*127), label="Ensemble 2 day",color=colors[1])
    # lines!(fig[1,1], LinRange(0, 10*365, 522), multi3 ./ (128*127), label="Online closure, batched 3 day, fewer initial conditions")
    lines!(ax3, LinRange(0, 10, 522), multi3more10 ./ (128*127), label="Ensemble 3 day",color=colors[2])

    lines!(ax3, LinRange(0, 10, 522), multi1010 ./ (128*127), label="Ensemble 10 day",color=colors[3])
    # lines!(ax3, LinRange(0, 10*365, 522), avg3_10 ./ (128*127), label="Averaged single initial conditions, 3 day optim.", color=colors[4])
    # lines!(ax3, LinRange(0, 10*365, 522), avg10_10 ./ (128*127), label="Averaged single initial conditions, 10 day optim.", color=colors[5])
    # lines!(ax3, LinRange(0, 10*365, 424), avg10single ./ (128*127), label="3 day optim.", color=colors[6])

    # lines!(fig[1,1], LinRange(0, 10*365, 522), multi20 ./ (128*127), label="Online closure, batched 20 day", color=:red3)
    # Legend(fig[2, 2], ax3)

    ax4 = Axis(fig[2,1],
        xlabel="Year",
        ylabel="Energy",
        title="Spatially averaged KE over 20 years"
    )
    lines!(ax4, LinRange(0, 20, 1043), multi2all ./ (128*127), label="Online closure, ensemble 2 day",color=colors[1])
    lines!(ax4, LinRange(0, 20, 1043), multi3all ./ (128*127), label="Online closure, ensemble 3 day",color=colors[2])
    lines!(ax4, LinRange(0, 20, 1043), multi10all ./ (128*127), label="Online closure, ensemble 3 day",color=colors[3])

    # lines!(ax3, LinRange(0, 10*365, 522),  cghr10_forplotting ./ (128*127), label="Filtered, coarse-grained 3.75 km",color=:red)

    # Legend(fig[3, 2], ax4)

    Legend(fig[3, 1], ax3, orientation = :horizontal)

    ga = fig[1, 1] = GridLayout()
    gb = fig[2, 1] = GridLayout()
    gc = fig[3, 1] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)"], [ga, gb])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end


    # Diverging results

    fig = Figure(size=(1000, 600), fontsize=15);
    ax = Axis(fig[1,1],
            # xlabel="Day",
            ylabel="Energy",
            title="Spatially averaged KE over 3 years"
    )
    lines!(ax, LinRange(0, 3, 1096),  hrcg[1:1096] ./ (128*127), label="Filtered, coarse-grained 3.75 km", color=:black)
    lines!(ax, LinRange(0, 3, 1096), noparam./ (128*127)), label="No closure, 30 km",color=:gray)
    lines!(ax, LinRange(0, 3, 1096), zb./ (128*127)), label="ZB20", color=:red)
    # lines!(ax, LinRange(0, 3*365, 1096), fiveday./ (128*127)), label="5 day")
    lines!(ax, LinRange(0, 3, 1096), tenday./ (128*127)), label="10 day", color=colors[4])
    # lines!(ax, LinRange(0, 3, 1096), threedayavg[1:1096] ./ (128*127)), label="Averaged 3 day")
    # lines!(ax, LinRange(0, 3, 1096), tendayavg[1:1096] ./ (128*127)), label="Averaged 3 day")
    lines!(ax,LinRange(0, 3, 1096), twentyday ./ (128*127), label="20 day", color=colors[5])
    # lines!(ax,LinRange(0, 3*365, 1096), twentydaycD ./ (128*127)), label="Online closure, 20 day with BD coeff")
    lines!(ax, LinRange(0, 3, 1096), thirtyday ./ (128*127), label="30 day",color=:purple)
    # lines!(ax, LinRange(0, 3*365, 1096), multi1 ./ (128*127)), label="Online closure, batched 1 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi1_more ./ (128*127)), label="Online closure, batched 1 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi2 ./ (128*127)), label="Online closure, ensemble 2 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi3more ./ (128*127)), label="Online closure, ensemble 3 day")
    # lines!(ax, LinRange(0, 3*365, 1096), multi5 ./ (128*127)), label="Online closure, batched 5 day", color=:mediumorchid)
    # lines!(ax, LinRange(0, 3*365, 1096), multi10 ./ (128*127)), label="Online closure, ensemble 10 day")#, color=:teal)
    lines!(ax, LinRange(0, 3, 1096), multi20 ./ (128*127), label="Ensemble 20 day", color=colors[7])
    # Legend(fig[1, 2], ax)

    ax3 = Axis(fig[2,1],
        xlabel="Year",
        ylabel="Energy",
        title="Spatially averaged KE over ten years"
    )
    lines!(ax3, LinRange(0, 10, 522),  hrcg_10 ./ (128*127), label="Filtered, coarse-grained 3.75 km",color=:black)
    lines!(ax3, LinRange(0, 10, 522), noparam10 ./ (128*127), label="No closure, 30 km",color=:gray)
    lines!(ax3, LinRange(0, 10, 522),  zb10 ./ (128*127), label="ZB20",color=:red)
    # lines!(ax3, LinRange(0, 10*365, 522), fiveday10./ (128*127), label="5 day")
    lines!(ax3, LinRange(0, 10, 522), tenday10./ (128*127), label="10 day",color=colors[4])
    lines!(ax3, LinRange(0, 10, 522), twentyday10./ (128*127), label="20 day",color=colors[5])
    lines!(ax3, LinRange(0, 10, 522), thirtyday10./ (128*127), label="30 day",color=:purple)
    lines!(ax3, LinRange(0, 10, 522), multi2010 ./ (128*127), label="Ensemble 20 day",color=colors[7])

    # Legend(fig[2, 2], ax3)

    Legend(fig[3, 1], ax, orientation = :horizontal)

    ga = fig[1, 1] = GridLayout()
    gb = fig[2, 1] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)"], [ga, gb])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end


    ax3 = Axis(fig[2,1],
        xlabel="Day",
        ylabel="J",
        title="Spatially averaged KE over ten years"
    )
    lines!(ax3, LinRange(0, 10*365, 522),  hrcg_10 ./ (128*127), label="Filtered, coarse-grained 3.75 km",color=:black)
    lines!(ax3, LinRange(0, 10*365, 522), noparam10 ./ (128*127), label="No closure, 30 km",color=:gray)
    lines!(ax3, LinRange(0, 10*365, 522),  zb10 ./ (128*127), label="ZB20",color=:red)
    # lines!(ax3, LinRange(0, 10*365, 522), fiveday10./ (128*127), label="5 day")
    lines!(ax3, LinRange(0, 10*365, 522), tenday10./ (128*127), label="10 day")
    lines!(ax3, LinRange(0, 10*365, 522), twentyday10./ (128*127), label="20 day")
    lines!(ax3, LinRange(0, 10*365, 522), thirtyday10./ (128*127), label="30 day")
    lines!(ax3, LinRange(0, 10*365, 522), multi2010 ./ (128*127), label="Ensemble 20 day")

    # Legend(fig[2, 2], ax3)

    Legend(fig[3, 1], ax, orientation = :horizontal)

    ga = fig[1, 1] = GridLayout()
    gb = fig[2, 1] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)"], [ga, gb])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # single initial condition runs

    # 10 day optimizations

    avg3_10 = cat(threedayavg[1:7:1096], threedayavg[1097:end];dims=1)
    avg10_10 = cat(tendayavg[1:7:1096], tendayavg[1097:end];dims=1)

    fig = Figure(size=(1000, 600), fontsize=15);
    ax = Axis(fig[1,1],
            # xlabel="Day",
            ylabel="Energy",
            title="10 day single i.c. runs, day 0 initial condition"
    )
    lines!(ax, LinRange(0, 10*365, 522), hrcg_10 ./ (128^2), label="Filtered, coarse-grained 3.75 km", color=:black)
    lines!(ax, LinRange(0, 10*365, 522), noparam10./ (128^2), label="No closure, 30 km",color=:gray)
    # lines!(ax, LinRange(0, 10*365, 522), zb10./ (128^2), label="ZB20", color=:red)
    # lines!(ax, LinRange(0, 10*365, 522), tenday10./ (128^2), label="10 day")
    # lines!(ax, LinRange(0, 10*365, 522), avg3_10 ./ (128^2), label="Averaged 3 day")
    lines!(ax, LinRange(0, 10*365, 522), avg10_10 ./ (128^2), label="Averaged 10 day")
    # lines!(ax, LinRange(0, 10*365, 522), single102[1:522] ./ 128^2, label="2")
    # lines!(ax, LinRange(0, 10*365, 522), single1020[1:522] ./ 128^2, label="20")
    # lines!(ax, LinRange(0, 10*365, 522), single1030[1:522] ./ 128^2, label="30")
    lines!(ax, LinRange(0, 10*365, 522), single1050[1:522] ./ 128^2, label="50")
    # lines!(ax, LinRange(0, 10*365, 522), single1060[1:522] ./ 128^2, label="60")
    # lines!(ax, LinRange(0, 10*365, 522), single1070[1:522] ./ 128^2, label="70")
    # lines!(ax, LinRange(0, 10*365, 522), single1080[1:522] ./ 128^2, label="80")
    Legend(fig[1, 2], ax)

    ax2 = Axis(fig[2,1],
            # xlabel="Day",
            ylabel="Energy",
            title="10 day single runs, day 50 initial condition"
    )
    lines!(ax2, LinRange(0, 10*365, 522), hrcg_10 ./ (128^2), label="Filtered, coarse-grained 3.75 km", color=:black)
    lines!(ax2, LinRange(0, 10*365, 522), noparam10./ (128^2), label="No closure, 30 km",color=:gray)
    # lines!(ax, LinRange(0, 10*365, 522), zb10./ (128^2), label="ZB20", color=:red)
    # lines!(ax, LinRange(0, 10*365, 522), tenday10./ (128^2), label="10 day")
    # lines!(ax, LinRange(0, 10*365, 522), avg3_10 ./ (128^2), label="Averaged 3 day")
    lines!(ax2, LinRange(0, 10*365, 522), avg10_10 ./ (128^2), label="Averaged 10 day")
    # lines!(ax2, LinRange(0, 10*365, 522), single1030_50[1:522] ./ 128^2, label="30")
    # lines!(ax2, LinRange(0, 10*365, 522), single1040_50[1:522] ./ 128^2, label="40")
    lines!(ax2, LinRange(0, 10*365, 522), single1050_50[1:522] ./ 128^2, label="50")
    Legend(fig[2, 2], ax2)

    # 3 day optimizations

    single31 = zeros(Float64, 783)
    single33 = zeros(Float64, 783)
    single38 = zeros(Float64, 783)
    single313 = zeros(Float64, 783)
    single323 = zeros(Float64, 783)
    single328 = zeros(Float64, 783)
    single333 = zeros(Float64, 783)
    single338 = zeros(Float64, 783)
    single341 = zeros(Float64, 783)
    single344 = zeros(Float64, 783)
    single348 = zeros(Float64, 783)
    single353 = zeros(Float64, 783)
    single358 = zeros(Float64, 783)
    single363 = zeros(Float64, 783)
    single368 = zeros(Float64, 783)
    single373 = zeros(Float64, 783)
    single378 = zeros(Float64, 783)
    single383 = zeros(Float64, 783)
    single386 = zeros(Float64, 783)

    avg3_10 = cat(threedayavg[1:7:1096], threedayavg[1097:end];dims=1)

    fig = Figure(size=(1000, 600), fontsize=15);
    ax = Axis(fig[1,1],
            # xlabel="Day",
            ylabel="Energy",
            title="3 day single i.c. runs, day 0 initial condition"
    )
    lines!(ax, LinRange(0, 10*365, 522), hrcg_10 ./ (128^2), label="Filtered, coarse-grained 3.75 km", color=:black)
    lines!(ax, LinRange(0, 10*365, 522), noparam10./ (128^2), label="No closure, 30 km",color=:gray)
    # lines!(ax, LinRange(0, 10*365, 522), zb10./ (128^2), label="ZB20", color=:red)
    # lines!(ax, LinRange(0, 10*365, 522), avg3_10 ./ (128^2), label="Averaged 3 day", color=:red)
    # lines!(ax, LinRange(0, 10*365, 522), single31[1:522] ./ 128^2, label="1")
    # lines!(ax, LinRange(0, 10*365, 522), single33[1:522] ./ 128^2, label="3")
    # lines!(ax, LinRange(0, 10*365, 522), single38[1:522] ./ 128^2, label="8")
    # lines!(ax, LinRange(0, 10*365, 522), single313[1:522] ./ 128^2, label="13")
    # lines!(ax, LinRange(0, 10*365, 522), single323[1:522] ./ 128^2, label="23")
    # lines!(ax, LinRange(0, 10*365, 522), single328[1:522] ./ 128^2, label="28")
    # lines!(ax, LinRange(0, 10*365, 522), single333[1:522] ./ 128^2, label="33")
    # lines!(ax, LinRange(0, 10*365, 522), single338[1:522] ./ 128^2, label="38")
    lines!(ax, LinRange(0, 10*365, 522), single341[1:522] ./ 128^2, label="41")
    lines!(ax, LinRange(0, 10*365, 522), single344[1:522] ./ 128^2, label="44")
    lines!(ax, LinRange(0, 10*365, 522), single348[1:522] ./ 128^2, label="48")
    lines!(ax, LinRange(0, 10*365, 522), single353[1:522] ./ 128^2, label="53")
    lines!(ax, LinRange(0, 10*365, 522), single358[1:522] ./ 128^2, label="58")
    lines!(ax, LinRange(0, 10*365, 522), single363[1:522] ./ 128^2, label="63")
    lines!(ax, LinRange(0, 10*365, 522), single368[1:522] ./ 128^2, label="68")
    lines!(ax, LinRange(0, 10*365, 522), single373[1:522] ./ 128^2, label="73")
    lines!(ax, LinRange(0, 10*365, 522), single378[1:522] ./ 128^2, label="78")
    lines!(ax, LinRange(0, 10*365, 522), single383[1:522] ./ 128^2, label="83")
    lines!(ax, LinRange(0, 10*365, 522), single386[1:522] ./ 128^2, label="86")
    Legend(fig[1, 2], ax)

    ax2 = Axis(fig[2,1],
            # xlabel="Day",
            ylabel="Energy",
            title="10 day single runs, day 50 initial condition"
    )
    lines!(ax2, LinRange(0, 10*365, 522), hrcg_10 ./ (128^2), label="Filtered, coarse-grained 3.75 km", color=:black)
    lines!(ax2, LinRange(0, 10*365, 522), noparam10./ (128^2), label="No closure, 30 km",color=:gray)
    # lines!(ax, LinRange(0, 10*365, 522), zb10./ (128^2), label="ZB20", color=:red)
    # lines!(ax, LinRange(0, 10*365, 522), tenday10./ (128^2), label="10 day")
    # lines!(ax, LinRange(0, 10*365, 522), avg3_10 ./ (128^2), label="Averaged 3 day")
    lines!(ax2, LinRange(0, 10*365, 522), avg10_10 ./ (128^2), label="Averaged 10 day")
    # lines!(ax2, LinRange(0, 10*365, 522), single1030_50[1:522] ./ 128^2, label="30")
    # lines!(ax2, LinRange(0, 10*365, 522), single1040_50[1:522] ./ 128^2, label="40")
    lines!(ax2, LinRange(0, 10*365, 522), single1050_50[1:522] ./ 128^2, label="50")
    Legend(fig[2, 2], ax2)

end