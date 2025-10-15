S = ShallowWaters.model_setup(T=Float32;
    output=false,
    L_ratio=1,
    g=9.81,
    H=500,
    wind_forcing_x="double_gyre",
    Lx=3840e3,
    seasonal_wind_x=false,
    topography="flat",
    bc="nonperiodic",
    bottom_drag="quadratic",
    tracer_advection=false,
    tracer_relaxation=false,
    zb_forcing_momentum=false,
    zb_forcing_dissipation=false,
    zb_filtered=true,
    nn_forcing_momentum=false,
    nn_forcing_dissipation=true,
    N=1,
    α=2,
    nx=128,
    Ndays=Ndays,
    initial_cond="ncfile",
    initpath="./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825",
    init_starti=1
)

T11 = SNN.Diag.CNNVars.T11
T12 = SNN.Diag.CNNVars.T12
T22 = SNN.Diag.CNNVars.T22

fig = Figure(fontsize = 15);
ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T11,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"T_{11}"),
colorrange=(-maximum(T11),
maximum(T11))
);
Colorbar(fig[1,2], hm1)

ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T12,
colormap=:balance,
title="All inputs (derivaties + velocities)",
axis=(xlabel="km", ylabel="km", title=L"T_{12}"),
colorrange=(-maximum(T12),
maximum(T12))
);
Colorbar(fig[1,4], hm2)

ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T22,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"T_{22}"),
colorrange=(-maximum(T22),
maximum(T22))
);
Colorbar(fig[2,2], hm3)

ζD_filtered = SZB.Diag.ZBVars.ζD_filtered
ζDhat_filtered = SZB.Diag.ZBVars.ζDhat_filtered
trace_filtered = SZB.Diag.ZBVars.trace_filtered

ζsqT = SZB.Diag.ZBVars.ζsqT
ζDT = SZB.Diag.ZBVars.ζDT
ζDhat = SZB.Diag.ZBVars.ζDhat

denom = SZB.grid.Δ^2 * SZB.grid.scale

fig = Figure(fontsize = 15);
ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
(ζsqT - ζDT)./denom,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\zeta^2 - \zeta D"),
colorrange=((-maximum((ζsqT - ζDT)./denom)),
maximum((ζsqT - ζDT)./denom))
);
Colorbar(fig[1,2], hm1)

ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
ζDhat ./ denom,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\zeta \hat{D}"),
colorrange=(-maximum(ζDhat ./ denom),
maximum(ζDhat ./ denom))
);
Colorbar(fig[1,4], hm2)

ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
(ζsqT + ζDT)./denom,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\zeta^2 + \zeta D"),
colorrange=(-maximum((ζsqT + ζDT)./denom),
maximum((ζsqT + ζDT)./denom))
);
Colorbar(fig[2,2], hm3)


S_u = SZB.Diag.ZBVars.S_u
S_v = SZB.Diag.ZBVars.S_v

fig = Figure(fontsize = 15);
ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
S_u,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"S_u"),
colorrange=((-maximum(S_u)),
maximum(S_u))
);
Colorbar(fig[1,2], hm1)

ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
S_v,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"S_v"),
colorrange=(-maximum(S_v),
maximum(S_v))
);
Colorbar(fig[1,4], hm2)