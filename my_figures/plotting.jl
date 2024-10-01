using NetCDF, Measures
using Plots, LaTeXStrings

##### NEED TO LOAD DATA ######


upadded = copy(uhr[:, :, end])
upadded = cat(zeros(2,1024),upadded,zeros(2,1024),dims=1)
upadded = cat(zeros(1023+4,2),upadded,zeros(1023+4,2),dims=2)

vpadded = copy(vhr[:, :, end])
v = cat(zeros(T,halo,nvy),v,zeros(T,halo,nvy),dims=1)
v = cat(zeros(T,nvx+2*halo,halo),v,zeros(T,nvx+2*halo,halo),dims=2)


dudy = zeros(1023, 1023)
dvdx = zeros(1023, 1023)

@inbounds for j ∈ 1:1023, i ∈ 1:1023
    dvdx[i,j] = vhr[i+1,j, end] - vhr[i,j, end]
end

@inbounds for j ∈ 1:1023, i ∈ 1:1023
    dudy[i,j] = uhr[i,j+1,end] - uhr[i,j,end]
end

vorticity = (dvdx - dudy) / 3750

vorticity_plot = heatmap(2:3.75:3840,
    2:3.75:3840,
    vorticity',
    c=:balance,
    clim=(-.0001, .0001),
    xlabel="x (km)",
    xguidefontsize=13,
    ylabel="y (km)",
    yguidefontsize=13,
    title="Relative vorticity",
    plot_titlefontsize=13,
    colorbar_title=L"1/s",
    colorbar_titlefontsize=13,
    dpi=300,
    size=(700, 600),
    margin=5mm
)

savefig("hr_vorticity_withscaling_050624.png")

energy = uhr[:, 1:1023, end].^2 + vhr[1:1023, :, end].^2 ./ 3750^2
energy_plot = heatmap(2:3.75:3840,
    2:3.75:3840,energy',
    clim = (0,3),
    c=:amp,
    xlabel="x (km)",
    xguidefontsize=13,
    ylabel="y (km)",
    yguidefontsize=13,
    title="Kinetic energy",
    plot_titlefontsize=13,
    colorbar_title=L"J",
    colorbar_titlefontsize=13,
    dpi=300
)

plot(energy_plot,vorticity_plot, layout=grid(1,2,
    widths=(4/8,4/8)),
    size=(1200,400),
    margin=5mm
)

savefig("energy_vorticity_noslip_withscaling_050624.png")

etaplot = heatmap(1:3.75:3840,
    1:3.75:3840,
    eta[:, :, end]',
    clim=(-4,4),
    c=:balance,
    xlabel="x (km)",
    xguidefontsize=13,
    ylabel="y (km)",
    yguidefontsize=13,
    title="Sea surface height",
    plot_titlefontsize=13,
    colorbar_title=L"m",
    colorbar_titlefontsize=13,
    dpi=300,
    size=(700, 600),
    colorbar=:top
)

y_u = 3750*Array(1:1024) .- 3750/2
Fx(y) = (64*3750*0.12/1000/500)*(cos.(2π*(y/3840e3 .- 1/2)) + 2*sin.(π*(y/3840e3 .- 1/2)))

wind_stress = plot(Fx(y_u), 
    y_u,
    title="Wind stress",
    xlabel="Pa",
    yticks=[],
    legend=false,
    dpi=300
)

plot(etaplot,wind_stress, layout=grid(1,2,
    widths=(6/8,2/8)),
    size=(1000,600),
    margin=5mm
)

savefig("ssh_windstress_newaxes_050624.png")

eta_avg = zeros(1024,1024)

for t = 1:3651
    eta_avg += etahr[:, :, t]
end

eta_avg = eta_avg ./ 3651

etaaveraged = heatmap(1:3.75:3840,
    1:3.75:3840,
    eta_avg[:, :, end]',
    clim=(-1,1),
    c=:balance,
    xlabel="x (km)",
    xguidefontsize=13,
    ylabel="y (km)",
    yguidefontsize=13,
    title="Time averaged sea surface height",
    plot_titlefontsize=13,
    colorbar_title=L"m",
    colorbar_titlefontsize=13,
    dpi=300,
    size=(700,600)
)

plot(etaaveraged,wind_stress, layout=grid(1,2,
    widths=(6/8,2/8)),
    size=(1000,600),
    margin=5mm
)

savefig("ssh_averaged_windstress.png")

savefig("time_averaged_ssh_050624.png")

using JLD2

energy_hr = []
energy_lr = []
energy_lr_param = []

