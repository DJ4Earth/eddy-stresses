using GLMakie

coarse_grained_hrstates = load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_3years_postspinup_dailysaves.jld2");
uhrcg = coarse_grained_hrstates[1];
vhrcg = coarse_grained_hrstates[2];
etahrcg = coarse_grained_hrstates[3];

u5eta = ncread("./5daystatewitheta/u.nc", "u");
u3 = ncread("./multi3day/u.nc", "u");
u10daynoeta = ncread("./dissipation_constant/results/128_online_stateweights_10dayoptimization_startfrom5daystate_noeta_oneyear/u.nc", "u")

u = ncread("./dissipation_constant/results/result_online_state_pluscD_weights_20dayoptimization_startfrom20day_3years_dailysaves/u.nc", "u")

umulti3 = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_10years_weeklysaves/u.nc", "u");

u20day = ncread("./dissipation_constant/results/result_online_stateweights_20dayoptimization_startfrom10day_fixedcfl_3years_dailysaves/u.nc", "u");
# v20day = ncread("./dissipation_constant/results/128_online_gelu_stateweights_20dayoptimization_startfrom10daystate_3years_dailysaves/v.nc", "v");
# eta20day = ncread("./dissipation_constant/results/128_online_gelu_stateweights_20dayoptimization_startfrom10daystate_3years_dailysaves/eta.nc", "eta");

u3smag = ncread("./dissipation_smagorinsky/results_with_parameterization/result_online_stateweights_3dayoptimization_smagorinskydiffusion_startfromoffline/u.nc", "u")

fig = Figure(fontsize=15);

framerate = 30
ax = Axis(fig[1, 1], xlabel="km", ylabel="km", title = "3.75 km x-velocity")

record(fig, "unoparam_zb20_multi2_3years_dailysaves.mp4", 1:1096) do t
    tempframe = uhr_all[:, :, t]
    hm =heatmap!(ax,
        LinRange(0, 3840, 1024),
        LinRange(0, 3840, 1024),
        tempframe,
        colormap = :balance,
        colorrange=(-maximum(abs.(uhr1[:,:,365])),maximum(abs.(uhr1[:,:,365])))
    )
    Colorbar(fig[1,2], hm, label="m/s")
end

fig = Figure(fontsize=15, size=(1000, 275));

Label(
    fig[0, 3],
    "u-velocity over three years",
    fontsize = 20,
    tellwidth = false
)
framerate = 30
ax1 = Axis(fig[1, 1], xlabel="km", ylabel="km", title = "No closure")
ax2 = Axis(fig[1, 3], xlabel="km", ylabel="km", title = "ZB20")
ax3 = Axis(fig[1, 5], xlabel="km", ylabel="km", title = "Ensemble 2 day")

record(fig, "unoparam_zb20_multi2_3years_dailysaves.mp4", 1:1096) do t
    hm = heatmap!(ax1,
        LinRange(0, 3840, 128),
        LinRange(0, 3840, 128),
        unoparam[:, :, t],
        colormap = :balance,
        colorrange=(-maximum(abs.(uhrcg1[:,:,365])),maximum(abs.(uhrcg1[:,:,365])))
    )
    Colorbar(fig[1,2], hm, label="m/s")
    hm = heatmap!(ax2,
        LinRange(0, 3840, 128),
        LinRange(0, 3840, 128),
        uzb[:, :, t],
        colormap = :balance,
        colorrange=(-maximum(abs.(uhrcg1[:,:,365])),maximum(abs.(uhrcg1[:,:,365])))
    )
    Colorbar(fig[1,4], hm, label="m/s")
    hm = heatmap!(ax3,
        LinRange(0, 3840, 128),
        LinRange(0, 3840, 128),
        umulti2[:, :, t],
        colormap = :balance,
        colorrange=(-maximum(abs.(uhrcg1[:,:,365])),maximum(abs.(uhrcg1[:,:,365])))
    )
    Colorbar(fig[1,6], hm, label="m/s")
end