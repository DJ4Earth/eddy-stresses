using NetCDF, Measures
using Plots, LaTeXStrings

etanoslipdissfilter1 = ncread("./data_files_gamma0.3/128_spinup_wforcing_dissipation_wfilter_1pass/eta.nc", "eta")

etanoslip = ncread("./data_files_gamma0.3/1024_spinup_noslip/eta.nc", "eta")

eta_avg = zeros(1024,1024)

for t = 1:3651
    eta_avg += etanoslip[:, :, t]
end

eta_avg = eta_avg ./ 3651

etaaveraged = heatmap(1:3.75:3840,
    1:3.75:3840,
    eta_avg[:, :, end]',
    clim=(-1,1),
    c=:balance,
    xlabel="x (km)",
    xguidefontsize=13,
    xtickfontsize=11,
    ylabel="y (km)",
    yguidefontsize=13,
    ytickfontsize=11,
    title="Time averaged sea surface height",
    plot_titlefontsize=13,
    colorbar_title=L"m",
    colorbar_titlefontsize=13,
    dpi=300,
    size=(700,600),
    colorbar=:top
)

y_u = 3750*Array(1:1024) .- 3750/2
Fx(y) = (64*3750*0.12/1000/500)*(cos.(2π*(y/3840e3 .- 1/2)) + 2*sin.(π*(y/3840e3 .- 1/2)))

wind_stress = plot(Fx(y_u), 
    y_u,
    title="Wind stress",
    xlabel="Pa",
    xguidefontsize=13,
    yticks=[],
    legend=false,
    dpi=300
)

plot(etaaveraged,wind_stress, layout=grid(1,2,
    widths=(6/8,2/8)),
    size=(1000,600),
    margin=5mm
)

uhr = ncread("./data_files_gamma0.3/1024_spinup_noslip/u.nc", "u")
vhr = ncread("./data_files_gamma0.3/1024_spinup_noslip/v.nc", "v")

dudy = zeros(1023, 1023)
dvdx = zeros(1023, 1023)

@inbounds for j ∈ 1:1023, i ∈ 1:1023
    dvdx[i,j] = vhr[i+1,j, end] - vhr[i,j, end]
end

@inbounds for j ∈ 1:1023, i ∈ 1:1023
    dudy[i,j] = uhr[i,j+1,end] - uhr[i,j,end]
end

vorticity = (dvdx - dudy) / 3750

vorticity_plot = heatmap(LinRange(0, 3840, 1023),
    LinRange(0, 3840, 1023),
    vorticity',
    c=:balance,
    clim=(-.0001, .0001),
    xlabel="x (km)",
    # xguidefontsize=14,
    # xtickfontsize=11,
    ylabel="y (km)",
    # yguidefontsize=14,
    # ytickfontsize=11,
    title="Relative vorticity",
    # plot_titlefontsize=13,
    colorbar_title=L"1/s",
    # colorbar_titlefontsize=14,
    dpi=300,
    margin=3mm,
    size=(500,400)
)

energy = uhr[:, 1:1023, end].^2 + vhr[1:1023, :, end].^2
energy_plot = heatmap(LinRange(0, 3840, 1023),
    LinRange(0, 3840, 1023),
    energy',
    c=:amp,
    clim = (0,3),
    xlabel="x (km)",
    xguidefontsize=13,
    xtickfontsize=11,
    ylabel="y (km)",
    yguidefontsize=13,
    ytickfontsize=11,
    title="Kinetic energy",
    plot_titlefontsize=13,
    colorbar_title=L"J",
    colorbar_titlefontsize=13,
    dpi=300,
    size=(700, 600)
)

plot(energy_plot,vorticity_plot, layout=grid(1,2,
    widths=(4/8,4/8)),
    size=(1200,400),
    margin=5mm
)

plot(etahr',
    c=:balance,
    clim=(-.1,.1),
    xlabel=L"x",
    xguidefontsize=13,
    ylabel=L"y",
    yguidefontsize=13,
    title="Sea surface height",
    plot_titlefontsize=13,
    colorbar_title=L"1/s",
    colorbar_titlefontsize=13,
    dpi=300
)

eta_avg = zeros(1024,1024)

for t = 1:3651
    eta_avg += etahr[:, :, t]
end

eta_avg = eta_avg ./ 3651

etalr = ncread("/Users/swilliamson/Documents/GitHub/ShallowWaters.jl/my_scripts/data_files_gamma0.3/128_spinup_noforcing_noslipbc/eta.nc", "eta")
etalr_withparameterization = ncread("/Users/swilliamson/Documents/GitHub/ShallowWaters.jl/my_scripts/data_files_gamma0.3/128_spinup_wforcing_dissipation_wfilter_1pass_noslipbc/eta.nc", "eta")

heatmap(etalr[:, :, end]',
    clim=(-2,2),
    c=:balance,
    xlabel=L"x",
    xguidefontsize=13,
    ylabel=L"y",
    yguidefontsize=13,
    title="Sea surface height",
    plot_titlefontsize=13,
    colorbar_title=L"1/s",
    colorbar_titlefontsize=13,
    dpi=300
)

etafreeslip = etanoslip
y_u = 3750*Array(1:1024) .- 3750/2
Fx(y) = (64*3750*0.12/1000/500)*(cos.(2π*(y/3840e3 .- 1/2)) + 2*sin.(π*(y/3840e3 .- 1/2)))

etaplot = heatmap(etafreeslip[:, :, end]',
    clim=(-2,2),
    c=:balance,
    xlabel=L"x",
    xguidefontsize=13,
    ylabel=L"y",
    yguidefontsize=13,
    title="Sea surface height",
    plot_titlefontsize=13,
    colorbar_title=L"1/s",
    colorbar_titlefontsize=13,
    dpi=300
)

wind_stress = plot(Fx(y_u), 
    y_u,
    title="Wind stress",
    xlabel="Pa",
    yticks=[],
    legend=false
)

plot(etaplot,wind_stress, layout=grid(1,2,
    widths=(6/8,2/8)),
    size=(950,400),
    margin=5mm
)

uhr = ncread("./data_files_gamma0.3/1024_spinup_noslip/u.nc", "u")
vhr = ncread("./data_files_gamma0.3/1024_spinup_noslip/v.nc", "v")

ulr = ncread("./data_files_gamma0.3/128_spinup_noforcing_noslipbc/u.nc", "u")
vlr = ncread("./data_files_gamma0.3/128_spinup_noforcing_noslipbc/v.nc", "v")

energy = ulr[:, 1:127, end].^2 + vlr[1:127, :, end].^2
energy_plot = heatmap(LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    energy',
    c=:amp,
    clim = (0,3),
    xlabel="x (km)",
    xguidefontsize=13,
    xtickfontsize=11,
    ylabel="y (km)",
    yguidefontsize=13,
    ytickfontsize=11,
    title="Kinetic energy",
    plot_titlefontsize=13,
    colorbar_title=L"J",
    colorbar_titlefontsize=13,
    dpi=300,
    size=(700, 600)
)

ulr_param = ncread("./data_files_gamma0.3/128_spinup_wforcing_dissipation_wfilter_1pass_noslipbc/u.nc", "u")
vlr_param = ncread("./data_files_gamma0.3/128_spinup_wforcing_dissipation_wfilter_1pass_noslipbc/v.nc", "v")

energy_hr = zeros(3651)
energy_lr = zeros(3657)
energy_lr_param = zeros(3657)

for j = 1:3651

    energy_hr[j] = (sum(uhr[:,:,j].^2) + sum(vhr[:,:,j].^2)) / (1024 * 1024)

end

for j = 1:3657

    energy_lr[j] = (sum(ulr[:,:,j].^2) + sum(vlr[:,:,j].^2)) / (128 * 128)
    energy_lr_param[j] = (sum(ulr_param[:,:,j].^2) + sum(vlr_param[:,:,j].^2)) / (128 * 128)

end

xlr = LinRange(0, 10, 3657)
xhr = LinRange(0, 10, 3651)

plot(xlr, energy_lr, label=("30km resolution, without parameterization"),dpi=300)
plot!(xlr, energy_lr_param, label=("30km resolution, with parameterization"),dpi=300)
plot!(xhr, energy_hr, label=("3.75km resolution"),dpi=300)
xlabel!("Years")
ylabel!("Spatially averaged energy")
title!("Energy during spinup")

eta_hr = ncread("./data_files_gamma0.3/1024_spinup_noslip/eta.nc", "eta");
eta_lr = ncread("./data_files_gamma0.3/128_spinup_noforcing_noslipbc/eta.nc", "eta");

one = heatmap(eta_hr[:, :, end]',
    clim=(-2,2),
    c=:balance,
    xlabel=L"x",
    xguidefontsize=13,
    ylabel=L"y",
    yguidefontsize=13,
    title="3.75km resolution ssh",
    plot_titlefontsize=13,
    colorbar_title=L"m",
    colorbar_titlefontsize=13,
    dpi=300
)

heatmap(ulr[:, :, end]',
    clim=(-2,2),
    c=:balance,
    xlabel=L"x",
    xguidefontsize=13,
    ylabel=L"y",
    yguidefontsize=13,
    title="3.75km resolution x-velocity",
    plot_titlefontsize=13,
    colorbar_title=L"m",
    colorbar_titlefontsize=13,
    dpi=300
)

two = heatmap(eta_lr[:, :, end]',
clim=(-2,2),
c=:balance,
xlabel=L"x",
xguidefontsize=13,
ylabel=L"y",
yguidefontsize=13,
title="30km resolution ssh",
plot_titlefontsize=13,
colorbar_title=L"m",
colorbar_titlefontsize=13,
dpi=300
)

plot(one, two, layout=grid(1,2,
    widths=(4/8,4/8)),
    size=(950,400),
    margin=5mm
)