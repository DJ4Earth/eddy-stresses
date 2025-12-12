using GLMakie

coarse_grained_hrstates = load_object("./spinup_files/1024_filtered_downsized_uveta_3years_postspinup_dailysaves.jld2");
uhrcg = coarse_grained_hrstates[1];
vhrcg = coarse_grained_hrstates[2];
etahrcg = coarse_grained_hrstates[3];

u5eta = ncread("./5daystatewitheta/u.nc", "u");
u3 = ncread("./multi3day/u.nc", "u");
u10daynoeta = ncread("./results/10daystate_noeta_oneyear/u.nc", "u")

v = ncread("./multi5day/v.nc", "v");
eta = ncread("./multi5day/eta.nc", "eta");

fig = Figure(fontsize=15);

framerate = 60
timestamps = range(1, 827, step=1)
ax = Axis(fig[1, 1], xlabel="km", ylabel="km", title = "v(x, y)")

record(fig, "output.mp4", 1:827) do t
    tempframe = v[:, :, t]
    hm =heatmap!(ax,
        LinRange(0, 3840, 128),
        LinRange(0, 3840, 128),
        tempframe,
        colormap = :balance,
        colorrange=(-maximum(abs.(vhrcg[:,:,10])),maximum(abs.(vhrcg[:,:,10])))
    )
    Colorbar(fig[1,2], hm, label="m/s")
end