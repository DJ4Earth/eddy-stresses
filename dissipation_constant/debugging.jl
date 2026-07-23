Plr = ShallowWaters.Parameter(T=Float64,
    output=false,
    L_ratio=1,
    g=9.81,
    H=500,
    cfl=.898,
    wind_forcing_x="double_gyre",
    Lx=3840e3,
    RKo=2,
    seasonal_wind_x=false,
    topography="flat",
    bc="nonperiodic",
    bottom_drag="quadratic",
    tracer_advection=false,
    tracer_relaxation=false,
    zb_forcing_momentum=false,
    zb_forcing_dissipation=false,
    zb_filtered=true,
    # nn_forcing_momentum=false,
    # nn_forcing_dissipation=true,
    N=1,
    α=2,
    nx=128,
    Ndays=2,
    initial_cond="rest",
    # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
);
Slr = ShallowWaters.model_setup(Plr);

Phr = ShallowWaters.Parameter(T=Float64,
    output=false,
    L_ratio=1,
    g=9.81,
    H=500,
    RKo=2,
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
    # nn_forcing_momentum=false,
    # nn_forcing_dissipation=true,
    N=1,
    α=2,
    nx=1024,
    Ndays=2,
    initial_cond="rest",
    # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
);
Shr = ShallowWaters.model_setup(Phr);

Seulerlr = deepcopy(Slr);
Seulerhr = deepcopy(Shr);
# Seulerhr = deepcopy(Shr);
Smom = deepcopy(Slr);

Plr2 = ShallowWaters.Parameter(T=Float64,
    output=false,
    L_ratio=1,
    g=1e-12,
    H=500,
    cfl=.898,
    wind_forcing_x="double_gyre",
    Lx=3840e3,
    RKo=2,
    ω=1e-16,
    seasonal_wind_x=false,
    topography="flat",
    bc="nonperiodic",
    bottom_drag="quadratic",
    tracer_advection=false,
    tracer_relaxation=false,
    zb_forcing_momentum=false,
    zb_forcing_dissipation=false,
    zb_filtered=true,
    # nn_forcing_momentum=false,
    # nn_forcing_dissipation=true,
    N=1,
    α=2,
    nx=128,
    Ndays=2,
    initial_cond="rest",
    # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
);
Slr2 = ShallowWaters.model_setup(Plr2);

Phr2 = ShallowWaters.Parameter(T=Float64,
    output=false,
    L_ratio=1,
    g=1e-12,
    H=500,
    wind_forcing_x="double_gyre",
    Lx=3840e3,
    RKo=2,
    ω=1e-12,
    seasonal_wind_x=false,
    topography="flat",
    bc="nonperiodic",
    bottom_drag="quadratic",
    tracer_advection=false,
    tracer_relaxation=false,
    zb_forcing_momentum=false,
    zb_forcing_dissipation=false,
    zb_filtered=true,
    # nn_forcing_momentum=false,
    # nn_forcing_dissipation=true,
    N=1,
    α=2,
    nx=1024,
    Ndays=2,
    initial_cond="rest",
    # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
);
Shr2 = ShallowWaters.model_setup(Phr2);

Sadvlr = deepcopy(Slr2);
Sadvhr = deepcopy(Shr2);

mom_u = zeros(Slr.grid.nux, Slr.grid.nuy);
mom_v = zeros(Slr.grid.nvx, Slr.grid.nvy);

adv_ulr = zeros(Slr.grid.nux, Slr.grid.nuy);
adv_vlr = zeros(Slr.grid.nvx, Slr.grid.nvy);

adv_uhr = zeros(Shr.grid.nux, Shr.grid.nuy);
adv_vhr = zeros(Shr.grid.nvx, Shr.grid.nvy);

du = zeros(Slr.grid.nux, Slr.grid.nuy);
dv = zeros(Slr.grid.nvx, Slr.grid.nvy);
deta = zeros(Slr.grid.nx, Slr.grid.ny);

duhr = zeros(Shr.grid.nux, Shr.grid.nuy);
dvhr = zeros(Shr.grid.nvx, Shr.grid.nvy);
detahr = zeros(Shr.grid.nx, Shr.grid.ny);

t = 225 * Slr.grid.dtint

n = 1096
u_, v_, eta_ = ShallowWaters.add_halo(uhrcgall[:,:,n], vhrcgall[:,:,n], etahrcgall[:,:,n], Slr)
uhr_, vhr_, etahr_ = ShallowWaters.add_halo(uhrall[:,:,n], vhrall[:,:,n], etahrall[:,:,n], Shr)

Seulerlr.Prog.u = u_
Seulerlr.Prog.v = v_
Seulerlr.Prog.η = eta_

Seulerhr.Prog.u = uhr_
Seulerhr.Prog.v = vhr_
Seulerhr.Prog.η = etahr_

Smom.Prog.u = u_
Smom.Prog.v = v_
Smom.Prog.η = eta_

Sadvlr.Prog.u = u_
Sadvlr.Prog.v = v_
Sadvlr.Prog.η = eta_

Sadvhr.Prog.u = uhr_
Sadvhr.Prog.v = vhr_
Sadvhr.Prog.η = etahr_

compute_tendencies_witheuler!(du, dv, deta, Seulerlr, n*t)
compute_tendencies_witheuler!(duhr, dvhr, detahr, Seulerhr, 8*n*t)
compute_momentum_new!(mom_u, mom_v, Smom, n*t)
compute_momentum_new!(adv_ulr, adv_vlr, Sadvlr, n*t)
compute_momentum_new!(adv_uhr, adv_vhr, Sadvhr, n*t*8)

ker = ImageFiltering.Kernel.gaussian((30e3/3750))

filteredu = imfilter(adv_uhr, reflect(ker))
filteredv = imfilter(adv_vhr, reflect(ker))

filtereddu = imfilter(duhr, reflect(ker))
filtereddv = imfilter(dvhr, reflect(ker))

adv_uhrcg = (filteredu[8:8:end, 4:8:end] .+ filteredu[8:8:end, 5:8:end]) .* 0.5
adv_vhrcg = (filteredv[4:8:end, 8:8:end] .+ filteredv[5:8:end, 8:8:end]) .* 0.5

du_hrcg = (filtereddu[8:8:end, 4:8:end] .+ filtereddu[8:8:end, 5:8:end]) .* 0.5
dv_hrcg = (filtereddv[4:8:end, 8:8:end] .+ filtereddv[5:8:end, 8:8:end]) .* 0.5

fig = Figure(size=(950, 500), fontsize=15);

Label(
    fig[0, 3],
    fontsize = 20,
    tellwidth = false
)

t = 1096

ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
((2 .* du) ./ 384),
colormap=:balance,
axis=(xlabel="km", ylabel="km", title="Euler"),
colorrange=(-1.5e-5, 1.5e-5)
);
Colorbar(fig[1,2], hm1)
hidexdecorations!(ax1)

ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
(mom_u),
colormap=:balance,
axis=(xlabel="km", ylabel="km", title="Momentum"),
colorrange=(-1.5e-5, 1.5e-5)
);
hidedecorations!(ax2)
Colorbar(fig[1,4], hm2)

ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
((2 .* du) ./ 384) .- mom_u,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title="Euler - Momentum"),
# colorrange=(-.005, .005)
);
hidedecorations!(ax3)
Colorbar(fig[1,6], hm3)

ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
(adv_ulr),
colormap=:balance,
axis=(xlabel="km", ylabel="km", title="Advection"),
colorrange=(-1.5e-5, 1.5e-5)
);
hidedecorations!(ax3)
Colorbar(fig[2,2], hm4)

ax5, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
mom_u .- adv_ulr,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title="Momentum - advec"),
colorrange=(-1.5e-5, 1.5e-5)
);
hidedecorations!(ax3)
Colorbar(fig[2,4], hm5)

ax6, hm6 = heatmap(fig[2,5], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
Slr.grid.Δ .* Suapprox[:,:,t],
colormap=:balance,
axis=(xlabel="km", ylabel="km", title="Nonlinear advection approximation"),
colorrange=(-1.5e-5, 1.5e-5)
);
hidedecorations!(ax3)
Colorbar(fig[2,6], hm6)

fig = Figure(size=(950, 500), fontsize=15);

Label(
    fig[0, 3],
    fontsize = 20,
    tellwidth = false
)

t = 1096

ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
((2 .* du) ./ 384),
colormap=:balance,
axis=(xlabel="km", ylabel="km", title="Euler"),
colorrange=(-1.5e-5, 1.5e-5)
);
Colorbar(fig[1,2], hm1)
hidexdecorations!(ax1)

ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
(adv_ulr),
colormap=:balance,
axis=(xlabel="km", ylabel="km", title="Advec"),
colorrange=(-2.5e-4, 2.5e-4)
);
hidedecorations!(ax2)
Colorbar(fig[1,4], hm2)

ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
((2 .* du) ./ 384) .- (adv_ulr),
colormap=:balance,
axis=(xlabel="km", ylabel="km", title="Approx"),
# colorrange=(-1.5e-5, 1.5e-5)
);
hidedecorations!(ax3)
Colorbar(fig[1,6], hm3)

ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
((2 .* du_hrcg) ./ 48 .- (2 .* du) ./ 384) .- (adv_uhrcg .- adv_ulr),
colormap=:balance,
axis=(xlabel="km", ylabel="km", title="Euler - Advec"),
colorrange=(-1.5e-5, 1.5e-5)
);
hidedecorations!(ax3)
Colorbar(fig[2,2], hm4)

ax5, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
((2 .* du_hrcg) ./ 48 .- (2 .* du) ./ 384) .- Slr.grid.Δ .* Suapprox[:,:,t],
colormap=:balance,
axis=(xlabel="km", ylabel="km", title="Euler - approx"),
colorrange=(-1.5e-5, 1.5e-5)
);
hidedecorations!(ax3)
Colorbar(fig[2,4], hm5)

ax6, hm6 = heatmap(fig[2,5], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
Slr.grid.Δ .* Suapprox[:,:,t],
colormap=:balance,
axis=(xlabel="km", ylabel="km", title="Nonlinear advection approximation"),
colorrange=(-1.5e-5, 1.5e-5)
);
hidedecorations!(ax3)
Colorbar(fig[2,6], hm6)