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

Seuler = deepcopy(Slr);
Sadv = deepcopy(Slr2);

t = 225 * Slr.grid.dtint
n = 1096

# offline check to make sure I'm computing everything correctly. The momentum budget terms should
# sum up to the euler tendencies I compute, du and dv
# check passed
function momentum_budget_offlinecheck()

    u_, v_, eta_ = ShallowWaters.add_halo(uhrcgall[:,:,n], vhrcgall[:,:,n], etahrcgall[:,:,n], Slr)

    Seuler.Prog.u = u_
    Seuler.Prog.v = v_
    Seuler.Prog.η = eta_

    Sadv.Prog.u = u_
    Sadv.Prog.v = v_
    Sadv.Prog.η = eta_

    du, dv, adv_u, adv_v, fu, fv, detagdx, detagdy = compute_mom_budget(Seuler, Sadv, t, n)

    detagdx2 = zeros(127,128)
    detagdy2 = zeros(128,127)

    ShallowWaters.∂x!(detagdx2, (Seuler.constants.g .* etahrcgall[:,:,n])./ Seuler.grid.Δ)
    ShallowWaters.∂y!(detagdy2, (Seuler.constants.g .* etahrcgall[:,:,n])./ Seuler.grid.Δ)

    # all of the rhs terms are computed with a *scaled* prognostic variable, so to actually 
    # see the momentum budget terms we need to remove that scaling. The division by Delta is because
    # all of the timestepping is done with dt / Delta, and we want the term that gets multiplied by dt
    # I should have corrected for the scalings in the values returned by compute_mom_budget already,
    # but for the stuff I'm pulling after (see below) I need to do it manually

    bottomdragu = (Seuler.Diag.Bottomdrag.Bu[2:end-1,2:end-1]) ./ (Seuler.grid.scale * Seuler.grid.Δ)
    bottomdragv = (Seuler.Diag.Bottomdrag.Bv[2:end-1,2:end-1]) ./ (Seuler.grid.scale * Seuler.grid.Δ)
    viscu = (Seuler.Diag.Smagorinsky.LLu1[:,2:end-1] + Seuler.Diag.Smagorinsky.LLu2[2:end-1,:]) ./ (Seuler.grid.scale * Seuler.grid.Δ);
    viscv = (Seuler.Diag.Smagorinsky.LLv1[:, 2:end-1] + Seuler.Diag.Smagorinsky.LLv2[2:end-1,:]) ./ (Seuler.grid.scale * Seuler.grid.Δ)

    Fx = Seuler.forcing.Fx ./ (Seuler.grid.scale * Seuler.grid.Δ)

    fig = Figure(size=(1050, 600), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ((2 .* du) ./ 384),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Euler tendency"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (adv_u + fv - detagdx + Fx + viscu + bottomdragu),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Sum of budget terms"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax2)
    Colorbar(fig[1,4], hm2)

    ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (adv_u),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Advection"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax2)
    Colorbar(fig[1,6], hm3)

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    fv,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Coriolis"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    hidedecorations!(ax4)
    Colorbar(fig[2,2], hm4)

    ax5, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -detagdx,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Gravity"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    hidedecorations!(ax5)
    Colorbar(fig[2,4], hm5)

    ax6, hm6 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Fx,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Wind stress"),
    # colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax6)
    Colorbar(fig[2,6], hm6)

    ax7, hm7 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    viscu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Viscosity"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    hidedecorations!(ax7)
    Colorbar(fig[3,2], hm7)

    ax8, hm8 = heatmap(fig[3,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    bottomdragu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Bottom drag"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    hidedecorations!(ax8)
    Colorbar(fig[3,4], hm8)

end
