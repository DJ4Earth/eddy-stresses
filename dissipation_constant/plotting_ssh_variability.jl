function ssh_variability()

    # std in the spatial dimension (resulting in a two-dimensional plot of time versus std)

    colors = Makie.wong_colors()

    # change this if you want ten year or 3 year
    index = vcat(1:7:1096, 1097:1461)
    # index = 1:1096
    # index = 1:366
    total = length(index)

    fig = Figure(size=(1100, 575), fontsize=15);
    fig.layout.alignmode = Outside();
    Label(
        fig[0, 1],
        "Sea-surface height variability (spatial)",
        fontsize = 20,
        tellwidth = false
    )
    ax = Axis(fig[1,1],xlabel="Day",ylabel="Standard-deviation")

    lines!(ax, LinRange(0,10*365, 522), std(etahrcgall[:,:,index], dims=[1,2])[1:522],label="Filtered, coarse-grained 3.75 km",color=:black)
    index=1:522
    lines!(ax, LinRange(0,10*365, 522), std(etanoparam10[:,:,index], dims=[1,2])[index], label="No closure, 30 km",color=:gray)
    lines!(ax, LinRange(0,10*365, 522), std(etazb10[:,:,index], dims=[1,2])[index], label="ZB20", color=:red)
    lines!(ax, LinRange(0,10*365, 522), std(etamulti210[:,:,index], dims=[1,2])[index], label="Ensemble 2 day", color=colors[1])
    lines!(ax, LinRange(0,10*365, 522), std(etamulti3more10[:,:,index], dims=[1,2])[index], label="Ensemble 3 day", color=colors[2])
    lines!(ax, LinRange(0,10*365, 522), std(etamulti1010[:,:,index], dims=[1,2])[index], label="Ensemble 10 day", color=colors[3])
    Legend(fig[2, 1], ax, orientation=:horizontal)

    # std in time (resulting in a three-dimensional plot std)

    colors = Makie.wong_colors()

    # change this if you want ten year or 3 year
    index = vcat(1:7:1096, 1097:1461)
    # index = 1:1096
    # index = 1:366

    fig = Figure(size=(850, 520), fontsize=15);
    fig.layout.alignmode = Outside();
    Label(
        fig[0, 2],
        "10-year sea-surface height variability",
        fontsize = 20,
        tellwidth = false
    )

    temp = std(etahrcgall[:,:,index], dims=[3])[:,:,1]

    ax, hm = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    temp,
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Filtered, \n coarse-grained 3.75 km"),
    colorrange=(0,maximum(abs.(temp)))
    );
    # Colorbar(fig[1,2], hm1, label="m")

    index=1:522
    ax2, hm2 = heatmap(fig[1,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    std(etanoparam10[:,:,index], dims=3)[:,:,1],
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="No closure, 30 km"),
    colorrange=(0,maximum(abs.(temp)))
    );
    # Colorbar(fig[1,4], hm1, label="m")

    ax2, hm2 = heatmap(fig[1,3],LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    std(etazb10[:,:,index], dims=3)[:,:,1],
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="ZB20"),
    colorrange=(0,maximum(abs.(temp)))
    );
    # Colorbar(fig[1,6], hm1, label="m")

    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    std(etamulti210[:,:,index], dims=3)[:,:,1],
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Ensemble 2 day"),
    colorrange=(0,maximum(abs.(temp)))
    );
    # Colorbar(fig[2,2], hm1, label="m")

    ax4, hm4 = heatmap(fig[2,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    std(etamulti3more10, dims=3)[:,:,1],
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Ensemble 3 day"),
    colorrange=(0,maximum(abs.(temp)))
    );
    # Colorbar(fig[2,4], hm4, label="m")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    std(etamulti1010, dims=3)[:,:,1],
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Ensemble 10 day"),
    colorrange=(0,maximum(abs.(temp)))
    );
    # Colorbar(fig[2,6], hm1, label="m")

    Colorbar(fig[1:2,4], hm, label="m")

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

    # absolute different in variability
    fig = Figure(size=(850, 520), fontsize=15);
    fig.layout.alignmode = Outside();
    Label(
        fig[0, 2],
        "Absolute difference in sea-surface height variability",
        fontsize = 20,
        tellwidth = false
    )

    temp = std(etahrcgall[:,:,index], dims=[3])[:,:,1]

    ax, hm = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    temp,
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Filtered, \n coarse-grained 3.75 km"),
    colorrange=(0,maximum(abs.(temp)))
    );
    # Colorbar(fig[1,2], hm1, label="m")
    hidexdecorations!(ax)

    index=1:522
    ax1, hm1 = heatmap(fig[1,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    abs.(temp .- std(etanoparam10[:,:,index], dims=3)[:,:,1]),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="No closure, 30 km"),
    colorrange=(0,maximum(abs.(temp)))
    );
    # Colorbar(fig[1,4], hm1, label="m")
    hidedecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,3],LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    abs.(temp .- std(etazb10[:,:,index], dims=3)[:,:,1]),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="ZB20"),
    colorrange=(0,maximum(abs.(temp)))
    );
    # Colorbar(fig[1,6], hm1, label="m")
    hidedecorations!(ax2)

    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    abs.(temp .- std(etamulti210[:,:,index], dims=3)[:,:,1]),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Ensemble 2 day"),
    colorrange=(0,maximum(abs.(temp)))
    );
    # Colorbar(fig[2,2], hm1, label="m")

    ax4, hm4 = heatmap(fig[2,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    abs.(temp .- std(etamulti3more10, dims=3)[:,:,1]),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Ensemble 3 day"),
    colorrange=(0,maximum(abs.(temp)))
    );
    # Colorbar(fig[2,4], hm4, label="m")
    hideydecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    abs.(temp .- std(etamulti1010, dims=3)[:,:,1]),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="Ensemble 10 day"),
    colorrange=(0,maximum(abs.(temp)))
    );
    # Colorbar(fig[2,6], hm1, label="m")
    hideydecorations!(ax5)

    Colorbar(fig[1:2,4], hm, label="m")

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

    # lines!(ax, LinRange(0, 3*365, 1096), fiveday./ (128^2), label="5 day")
    # lines!(ax, LinRange(0, 3*365, 1096), tenday./ (128^2), label="10 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentyday ./ (128^2), label="20 day")
    # lines!(ax,LinRange(0, 3*365, 1096), twentydaycD ./ (128^2), label="Online closure, 20 day with BD coeff")
    # lines!(ax, LinRange(0, 3*365, 1096), thirtyday ./ (128^2), label="30 day")


end