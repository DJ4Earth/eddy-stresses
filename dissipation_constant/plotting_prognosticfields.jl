function prognostic_plots()

    # Prognostic variables #############################################################

    # high versus low resolution eta
    t = 1096
    fig = Figure(size=(775, 300), fontsize=15);
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 1024),
    LinRange(0, 3840, 1024),
    etahrall[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\text{3.75 km } \mathbf{\eta}(t, x, y)"),
    colorrange=(-3,3)
    );
    Colorbar(fig[1,2], hm1, label="m")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    etanoparam[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\text{30 km } \mathbf{\eta}(t, x, y)"),
    colorrange=(-3,3)
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

    # time-averaged eta fields

    # change this if you want ten year or 3 year
    index = vcat(1:7:1096, 1097:1461)
    # index = 1:1096
    # index = 1:366
    total = length(index)
    timeavg_hr = sum(etahrcgall[:,:,index],dims=3)[:,:,1] ./ total

    fig = Figure(size=(850, 520), fontsize=15);
    fig.layout.alignmode = Outside();
    Label(
        fig[0, 2],
        "10-year averaged sea-surface height",
        fontsize = 20,
        tellwidth = false
    )

    ax, hm = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    timeavg_hr,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Filtered, \n coarse-grained 3.75 km"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    # Colorbar(fig[1,2], hm, label="m")
    hidexdecorations!(ax)

    index=1:522
    ax1, hm1 = heatmap(fig[1,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(etanoparam10[:,:,index], dims=3)[:,:,1] ./ total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="No closure, 30 km"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    # Colorbar(fig[1,4], hm1, label="m")
    hidedecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(etazb10[:,:,index], dims=3)[:,:,1] ./ total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="ZB20"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    # Colorbar(fig[1,5], hm1, label="m")
    hidedecorations!(ax2)

    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(etamulti210[:,:,index], dims=3)[:,:,1] ./ total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 2 day"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    # Colorbar(fig[2,2], hm1, label="m")

    ax4, hm4 = heatmap(fig[2,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(etamulti3more10, dims=3)[:,:,1] ./ total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 3 day"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    # Colorbar(fig[2,4], hm4, label="m")
    hideydecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(etamulti1010, dims=3)[:,:,1] ./ total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 10 day"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    # Colorbar(fig[2,6], hm1, label="m")
    hideydecorations!(ax5)

    Colorbar(fig[1:2,4], hm5, label="m")

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

    # 20 year average ssh, multi2 and multi3 results

    fig = Figure(size=(600, 300), fontsize=15);
    fig.layout.alignmode = Outside();
    Label(
        fig[0, 1:2],
        "20-year averaged sea-surface height",
        fontsize = 15,
        tellwidth = false
    )

    index = vcat(1:7:1096, 1096:1983)
    total = length(index)

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(etamulti2all[:,:,index], dims=3)[:,:,1]/total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 2 day"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    # Colorbar(fig[1,2], hm1, label="m")

    index=1:522
    ax2, hm2 = heatmap(fig[1,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(etamulti3all[:,:,index], dims=3)[:,:,1] ./ total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 3 day"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    # Colorbar(fig[1,4], hm1, label="m")

    Colorbar(fig[1,3], hm, label="m")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 2] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)"], [ga, gb])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # absolute difference in time-averaged ssh fields
    index = vcat(1:7:1096, 1097:1461)
    # index = 1:1096
    # index = 1:366
    total = length(index)
    timeavg_hr = sum(etahrcgall[:,:,index],dims=3)[:,:,1] ./ total

    fig = Figure(size=(850, 520), fontsize=15);
    fig.layout.alignmode = Outside();
    Label(
        fig[0, 3],
        "Difference in 10-year averaged sea-surface height",
        fontsize = 20,
        tellwidth = false
    )

    ax, hm = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    timeavg_hr,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Filtered, \n coarse-grained 3.75 km"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[1,2], hm, label="m")
    hidexdecorations!(ax)

    index=1:522
    ax1, hm1 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (timeavg_hr .- sum(etanoparam10[:,:,index], dims=3)[:,:,1] ./ total),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="No closure, 30 km"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[1,4], hm1, label="m")
    hidedecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (timeavg_hr .- sum(etazb10[:,:,index], dims=3)[:,:,1] ./ total),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="ZB20"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[1,6], hm1, label="m")
    hidedecorations!(ax2)

    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (timeavg_hr .- sum(etamulti210[:,:,index], dims=3)[:,:,1] ./ total),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 2 day"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[2,2], hm1, label="m")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (timeavg_hr .- sum(etamulti3more10, dims=3)[:,:,1] ./ total),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 3 day"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[2,4], hm4, label="m")
    hideydecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (timeavg_hr .- sum(etamulti1010, dims=3)[:,:,1] ./ total),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 10 day"),
    colorrange=(-maximum(abs.(timeavg_hr)),maximum(abs.(timeavg_hr)))
    );
    Colorbar(fig[2,6], hm1, label="m")
    hideydecorations!(ax5)

    # Colorbar(fig[2:3,4], hm, label="m")

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

function generalizability()

    fig = Figure(size=(600, 900), fontsize=15);
    fig.layout.alignmode = Outside();
    # Label(
    #     fig[0, 1:2],
    #     "3-year averaged sea-surface height",
    #     fontsize = 15,
    #     tellwidth = false
    # )
    Label(fig[0, 1], "Coarse-grained, 3.75 km",tellwidth=false)
    Label(fig[0, 2], "Multi 3",tellwidth=false)

    index = vcat(1:1096)
    total = length(index)

    cghr = sum(etahrcgall[:,:,1:1096], dims=3)[:,:,1]/total
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    cghr,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="45 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );
    # Colorbar(fig[1,2], hm1, label="m")
    hidexdecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(etamulti3more[:,:,1:1096], dims=3)[:,:,1]/total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="45 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );
    # Colorbar(fig[1,2], hm1, label="m")
    hidedecorations!(ax2)

    cghr50 = sum(etahrcg_50lat[:,:,1:1096], dims=3)[:,:,1]/total
    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    cghr50,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="50 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );
    hidexdecorations!(ax3)

    ax2, hm2 = heatmap(fig[2,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(eta_multi3_50lat[:,:,1:24:(24*365*3)], dims=3)[:,:,1]/total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="50 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );
    hidedecorations!(ax2)

    cghr0 = sum(etahrcg_0lat[:,:,1:1096], dims=3)[:,:,1]/total
    ax33, hm33 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    cghr0,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="0 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );
    hidexdecorations!(ax33)

    ax22, hm22 = heatmap(fig[3,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(eta_multi3_0lat[:,:,1:1096], dims=3)[:,:,1]/total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="0 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );
    hidedecorations!(ax22)

    cghrneg60 = sum(etahrcg_neg60lat[:,:,1:1096], dims=3)[:,:,1]/total
    ax3, hm3 = heatmap(fig[4,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    cghrneg60,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="-60 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );

    ax2, hm2 = heatmap(fig[4,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    sum(eta_multi3_neg60lat[:,:,1:24:(24*365*3)], dims=3)[:,:,1]/total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="-60 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );
    hideydecorations!(ax2)
    Colorbar(fig[1:4, 3], hm2, label="m")

    # difference plots

    fig = Figure(size=(600, 900), fontsize=15);
    fig.layout.alignmode = Outside();
    # Label(
    #     fig[0, 1:2],
    #     "3-year averaged sea-surface height",
    #     fontsize = 15,
    #     tellwidth = false
    # )
    Label(fig[0, 1], "Coarse-grained, 3.75 km",tellwidth=false)
    Label(fig[0, 2], "Multi 3",tellwidth=false)

    index = vcat(1:1096)
    total = length(index)

    cghr = sum(etahrcgall[:,:,1:1096], dims=3)[:,:,1]/total

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    cghr,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="45 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );
    # Colorbar(fig[1,2], hm1, label="m")
    hidexdecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    cghr .- sum(etamulti3more[:,:,1:1096], dims=3)[:,:,1]/total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="45 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );
    # Colorbar(fig[1,2], hm1, label="m")
    hidedecorations!(ax2)

    cghr50 = sum(etahrcg_50lat[:,:,1:1096], dims=3)[:,:,1]/total
    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    cghr50,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="50 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );
    hidexdecorations!(ax3)

    ax2, hm2 = heatmap(fig[2,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    cghr50 .- sum(eta_multi3_50lat[:,:,1:24:(24*365*3)], dims=3)[:,:,1]/total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="50 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );
    hidedecorations!(ax2)

    cghr0 = sum(etahrcg_0lat[:,:,1:1096], dims=3)[:,:,1]/total
    ax33, hm33 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    cghr0,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="0 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );
    hidexdecorations!(ax33)

    ax22, hm22 = heatmap(fig[3,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    cghr0 .- sum(eta_multi3_0lat[:,:,1:1096], dims=3)[:,:,1]/total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="0 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );
    hidedecorations!(ax22)

    cghrneg60 = sum(etahrcg_neg60lat[:,:,1:1096], dims=3)[:,:,1]/total
    ax3, hm3 = heatmap(fig[4,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    cghrneg60,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="-60 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );

    ax2, hm2 = heatmap(fig[4,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    cghrneg60 .- sum(eta_multi3_neg60lat[:,:,1:24:(24*365*3)], dims=3)[:,:,1]/total,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="-60 degree latitude"),
    colorrange=(-maximum(abs.(cghr)),maximum(abs.(cghr)))
    );
    hideydecorations!(ax2)
    Colorbar(fig[1:4, 3], hm2, label="m")


end
