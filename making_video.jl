using GLMakie

coarse_grained_hrstates = load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_3years_postspinup_dailysaves.jld2");
uhrcg = coarse_grained_hrstates[1];
vhrcg = coarse_grained_hrstates[2];
etahrcg = coarse_grained_hrstates[3];

u5eta = ncread("./5daystatewitheta/u.nc", "u");
u3 = ncread("./multi3day/u.nc", "u");
u10daynoeta = ncread("./dissipation_constant/results/128_online_stateweights_10dayoptimization_startfrom5daystate_noeta_oneyear/u.nc", "u")

u = ncread("./dissipation_constant/results/result_online_state_pluscD_weights_20dayoptimization_startfrom20day_3years_dailysaves/u.nc", "u")

umulti3 = ncread("./dissipation_constant/results/128_online_multistateweights_3dayoptimization_3-25-30-40-50-60-80initdays_startfrom20daystate_3years_dailysaves/u.nc", "u")

u20day = ncread("./dissipation_constant/results/result_online_stateweights_20dayoptimization_startfrom10day_fixedcfl_3years_dailysaves/u.nc", "u");
# v20day = ncread("./dissipation_constant/results/128_online_gelu_stateweights_20dayoptimization_startfrom10daystate_3years_dailysaves/v.nc", "v");
# eta20day = ncread("./dissipation_constant/results/128_online_gelu_stateweights_20dayoptimization_startfrom10daystate_3years_dailysaves/eta.nc", "eta");

u3smag = ncread("./dissipation_smagorinsky/results_with_parameterization/result_online_stateweights_3dayoptimization_smagorinskydiffusion_startfromoffline/u.nc", "u")

fig = Figure(fontsize=15);

framerate = 40
timestamps = range(1, 1096, step=1)
ax = Axis(fig[1, 1], xlabel="km", ylabel="km", title = "u(x, y)")

record(fig, "u20day.mp4", 1:1096) do t
    tempframe = u20day[:, :, t]
    hm =heatmap!(ax,
        LinRange(0, 3840, 128),
        LinRange(0, 3840, 128),
        tempframe,
        colormap = :balance,
        colorrange=(-maximum(abs.(uhrcg[:,:,365])),maximum(abs.(uhrcg[:,:,365])))
    )
    Colorbar(fig[1,2], hm, label="m/s")
end