
function appendix_plots()

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