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

fframerate = 15

xs = LinRange(0, 3840, 128)
ys = LinRange(0, 3840, 128)

inv_scale = 1 / Shr.constants.scale
# clim = maximum(abs.(ζhr[:,:,1000] .* inv_scale))  # compute ONCE
clim = .0001

ax = Axis(fig[1, 1])

hm = heatmap!(
    ax,
    xs, ys,
    etamulti2[:, :, 1],
    colormap = :balance,
    colorrange = (-3, 3)
)
hidedecorations!(ax) 

Colorbar(fig[1, 2], hm, label = "1/s")

record(fig, "ssh_threeyears_multi2.mp4", 1:1096; framerate = framerate) do t
    hm[3][] = etamulti2[:, :, t]
end


# double video 
fig = Figure(size=(900,350));
framerate = 15

xs = LinRange(0, 3840, 129)
ys = LinRange(0, 3840, 129)

xsl = LinRange(0, 3840, 129)
ysl = LinRange(0, 3840, 129)

clim = .0001

ax1 = Axis(fig[1, 1],
    xlabel = "km",
    ylabel = "km",
    title  = "30 km vorticity, with ensemble 2 day closure"
)

hm1 = heatmap!(
    ax1,
    xs, ys,
    ζmulti2[:, :, 1],
    colormap = :balance,
    colorrange = (-clim, clim)
)

ax2 = Axis(fig[1, 3],
    xlabel = "km",
    ylabel = "km",
    title  = "30 km vorticity, with ensemble 3 day closure"
)

hm2 = heatmap!(
    ax2,
    xs, ys,
    ζmulti3more[:, :, 1],
    colormap = :balance,
    colorrange = (-clim, clim)
)

Colorbar(fig[1, 2], hm1, label = "1/s")
Colorbar(fig[1, 4], hm2, label = "1/s")


record(fig, "multi2_multi3_vorticity.mp4", 1:1096; framerate = framerate) do t
    hm1[3][] = ζmulti2[:, :, t]
    hm2[3][] = ζmulti3more[:, :, t]
end


# triple video

fig = Figure(size=(900,350));
framerate = 15

xs = LinRange(0, 3840, 129)
ys = LinRange(0, 3840, 129)

xsl = LinRange(0, 3840, 129)
ysl = LinRange(0, 3840, 129)

clim = .0001

Label(
    fig[0, 2],
    "30 km vorticity",
    fontsize = 25,
    tellwidth = false
)

ax0 = Axis(fig[1, 1],
    xlabel = "km",
    ylabel = "km",
    title  = "Data"
)

hm0 = heatmap!(
    ax0,
    xs, ys,
    ζhrcg[:, :, 1],
    colormap = :balance,
    colorrange = (-clim, clim)
)

ax1 = Axis(fig[1, 2],
    xlabel = "km",
    ylabel = "km",
    title  = "Non-parameterized model"
)

hm1 = heatmap!(
    ax1,
    xs, ys,
    ζnoparam[:, :, 1],
    colormap = :balance,
    colorrange = (-clim, clim)
)
hideydecorations!(ax1)

ax2 = Axis(fig[1, 3],
    xlabel = "km",
    ylabel = "km",
    title  = "Parameterized model"
)

hm2 = heatmap!(
    ax2,
    xs, ys,
    ζmulti2[:, :, 1],
    colormap = :balance,
    colorrange = (-clim, clim)
)
hideydecorations!(ax2)

# Colorbar(fig[1, 2], hm1, label = "1/s")
# Colorbar(fig[1, 4], hm2, label = "1/s")
Colorbar(fig[1, 4], hm1, label = "1/s")

# Colorbar(fig[1, 4], hm2, label = "1/s")

record(fig, "data_noparam_multi2_vorticity.mp4", 1:1096; framerate = framerate) do t
    hm0[3][] = ζhrcg[:, :, t]
    hm1[3][] = ζnoparam[:, :, t]
    hm2[3][] = ζmulti2[:, :, t]
end