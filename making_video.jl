using GLMakie

coarse_grained_hrstates = load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_3years_postspinup_dailysaves.jld2");
uhrcg = coarse_grained_hrstates[1];
vhrcg = coarse_grained_hrstates[2];
etahrcg = coarse_grained_hrstates[3];

u5eta = ncread("./5daystatewitheta/u.nc", "u");
u3 = ncread("./multi3day/u.nc", "u");
u10daynoeta = ncread("./dissipation_constant/results/128_online_stateweights_10dayoptimization_startfrom5daystate_noeta_oneyear/u.nc", "u")

v = ncread("./multi5day/v.nc", "v");
eta = ncread("./multi5day/eta.nc", "eta");

u20day = ncread("./dissipation_constant/results/128_online_gelu_stateweights_20dayoptimization_startfrom10daystate_3years_dailysaves/u.nc", "u");
v20day = ncread("./dissipation_constant/results/128_online_gelu_stateweights_20dayoptimization_startfrom10daystate_3years_dailysaves/v.nc", "v");
eta20day = ncread("./dissipation_constant/results/128_online_gelu_stateweights_20dayoptimization_startfrom10daystate_3years_dailysaves/eta.nc", "eta");


fig = Figure(fontsize=15);

framerate = 40
timestamps = range(1, 1098, step=1)
ax = Axis(fig[1, 1], xlabel="km", ylabel="km", title = "η(x, y)")

record(fig, "output.mp4", 1:1098) do t
    tempframe = eta20day[:, :, t]
    hm =heatmap!(ax,
        LinRange(0, 3840, 128),
        LinRange(0, 3840, 128),
        tempframe,
        colormap = :balance,
        colorrange=(-maximum(abs.(etahrcg[:,:,365])),maximum(abs.(etahrcg[:,:,365])))
    )
    Colorbar(fig[1,2], hm, label="m")
end