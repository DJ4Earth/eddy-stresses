# offline check to make sure I'm computing everything correctly. The momentum budget terms should
# sum up to the euler tendencies I compute, du and dv
# check passed
function compute_timeaveraged_momentumbudget(u, v, eta, Seuler, Sadv)

    dutotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)
    adv_utotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)
    fvtotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)
    detagdxtotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)
    bottomdragutotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)
    viscutotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)
    Fxtotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)
    Sutotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)

    dvtotal = zeros(Seuler.grid.nvx, Seuler.grid.nvy)
    adv_vtotal = zeros(Seuler.grid.nvx, Seuler.grid.nvy)
    futotal = zeros(Seuler.grid.nvx, Seuler.grid.nvy)
    detagdytotal = zeros(Seuler.grid.nvx, Seuler.grid.nvy)
    bottomdragvtotal = zeros(Seuler.grid.nvx, Seuler.grid.nvy)
    viscvtotal = zeros(Seuler.grid.nvx, Seuler.grid.nvy)
    Svtotal = zeros(Seuler.grid.nvx, Seuler.grid.nvy)

    for n = 1:1096

        t = 225 * Seuler.grid.dtint * n

        S = deepcopy(Seuler)
        S2 = deepcopy(Sadv)

        u_, v_, eta_ = ShallowWaters.add_halo(u[:,:,n], v[:,:,n], eta[:,:,n], Seuler)

        S.Prog.u = u_
        S.Prog.v = v_
        S.Prog.η = eta_

        S2.Prog.u = u_
        S2.Prog.v = v_
        S2.Prog.η = eta_

        du, dv, adv_u, adv_v, fu, fv, detagdx, detagdy, bottomdragu, bottomdragv, viscu, viscv, Fx = compute_mom_budget(S, S2, t, n);

        if S.parameters.nn_forcing_dissipation
            Sutotal += S.Diag.CNNVars.S_u ./ (S.grid.Δ * S.constants.scale)
            Svtotal += S.Diag.CNNVars.S_v ./ (S.grid.Δ * S.constants.scale)
        end

        dutotal += du
        adv_utotal += adv_u
        fvtotal += fv
        detagdxtotal += detagdx
        bottomdragutotal += bottomdragu
        viscutotal += viscu
        Fxtotal += Fx

        dvtotal += dv
        adv_vtotal += adv_v
        futotal += fu
        detagdytotal += detagdy
        bottomdragvtotal += bottomdragv
        viscvtotal += viscv

    end

    scale = 1 / 1096
    return dutotal .* scale, dvtotal .* scale, adv_utotal .* scale, adv_vtotal .* scale, futotal .* scale, fvtotal .* scale, detagdxtotal .* scale, detagdytotal .* scale, bottomdragutotal .* scale, bottomdragvtotal .* scale, viscutotal .* scale, viscvtotal .* scale, Fxtotal .* scale, Sutotal .* scale, Svtotal .* scale

end

function compute_timeaveraged_energybudget(u, v, eta, Seuler, Sadv)

    udutotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)
    uadv_utotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)
    ufvtotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)
    udetagdxtotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)
    ubottomdragutotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)
    uviscutotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)
    uFxtotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)
    Sutotal = zeros(Seuler.grid.nux, Seuler.grid.nuy)

    vdvtotal = zeros(Seuler.grid.nvx, Seuler.grid.nvy)
    vadv_vtotal = zeros(Seuler.grid.nvx, Seuler.grid.nvy)
    vfutotal = zeros(Seuler.grid.nvx, Seuler.grid.nvy)
    vdetagdytotal = zeros(Seuler.grid.nvx, Seuler.grid.nvy)
    vbottomdragvtotal = zeros(Seuler.grid.nvx, Seuler.grid.nvy)
    vviscvtotal = zeros(Seuler.grid.nvx, Seuler.grid.nvy)
    Svtotal = zeros(Seuler.grid.nvx, Seuler.grid.nvy)

    for n = 1:1096

        t = 225 * Seuler.grid.dtint * n

        S = deepcopy(Seuler)
        S2 = deepcopy(Sadv)

        u_, v_, eta_ = ShallowWaters.add_halo(u[:,:,n], v[:,:,n], eta[:,:,n], Seuler)

        S.Prog.u = u_
        S.Prog.v = v_
        S.Prog.η = eta_

        S2.Prog.u = u_
        S2.Prog.v = v_
        S2.Prog.η = eta_

        du, dv, adv_u, adv_v, fu, fv, detagdx, detagdy, bottomdragu, bottomdragv, viscu, viscv, Fx = compute_mom_budget(S, S2, t, n);


        if S.parameters.nn_forcing_dissipation
            Sutotal += u[:,:,n] .* S.Diag.CNNVars.S_u ./ (S.grid.Δ * S.constants.scale)
            Svtotal += v[:,:,n] .* S.Diag.CNNVars.S_v ./ (S.grid.Δ * S.constants.scale)
        end

        udutotal += u[:,:,n] .* du
        uadv_utotal += u[:,:,n] .* adv_u
        ufvtotal += u[:,:,n] .* fv
        udetagdxtotal += u[:,:,n] .* detagdx
        ubottomdragutotal += u[:,:,n] .* bottomdragu
        uviscutotal += u[:,:,n] .* viscu
        uFxtotal += u[:,:,n] .* Fx

        vdvtotal += v[:,:,n] .* dv
        vadv_vtotal += v[:,:,n] .* adv_v
        vfutotal += v[:,:,n] .* fu
        vdetagdytotal += v[:,:,n] .* detagdy
        vbottomdragvtotal += v[:,:,n] .* bottomdragv
        vviscvtotal += v[:,:,n] .* viscv

    end

    scale = 1 / 1096
    return udutotal .* scale, vdvtotal .* scale, uadv_utotal .* scale, vadv_vtotal .* scale, vfutotal .* scale, ufvtotal .* scale, udetagdxtotal .* scale, vdetagdytotal .* scale, ubottomdragutotal .* scale, vbottomdragvtotal .* scale, uviscutotal .* scale, vviscvtotal .* scale, uFxtotal .* scale,  Sutotal .* scale, Svtotal .* scale

end

function momentum_budget_offline()

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
        g=1e-20,
        H=500,
        cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        ω=1e-22,
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

    u_, v_, eta_ = ShallowWaters.add_halo(uhrcgall[:,:,n], vhrcgall[:,:,n], etahrcgall[:,:,n], Slr)

    Seuler.Prog.u = u_
    Seuler.Prog.v = v_
    Seuler.Prog.η = eta_

    Sadv.Prog.u = u_
    Sadv.Prog.v = v_
    Sadv.Prog.η = eta_

    du, dv, adv_u, adv_v, fu, fv, detagdx, detagdy, bottomdragu, bottomdragv, viscu, viscv, Fx = compute_mom_budget(Seuler, Sadv, t, n);
    dutotal, dvtotal, adv_utotal, adv_vtotal, futotal, fvtotal, detagdxtotal, detagdytotal, bottomdragutotal, bottomdragvtotal, viscutotal, viscvtotal, Fxtotal, _, _ = compute_timeaveraged_momentumbudget(uhrcgall, vhrcgall, etahrcgall, Seuler, Sadv);
    udutotal, vdvtotal, uadv_utotal, vadv_vtotal, vfutotal, ufvtotal, udetagdxtotal, vdetagdytotal, ubottomdragutotal, vbottomdragvtotal, uviscutotal, vviscvtotal, uFxtotal, _, _ = compute_timeaveraged_energybudget(uhrcgall, vhrcgall, etahrcgall, Seuler, Sadv);

    fig = Figure(fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ((2 .* dutotal) ./ 384) .- (adv_utotal + fvtotal - detagdxtotal + Fxtotal + viscutotal + bottomdragutotal),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Euler tendency - sum of budget terms"),
    colorrange=(-2.5e-15, 2.5e-15)
    );
    Colorbar(fig[1,2], hm1)

    # fig = Figure(size=(1050, 700), fontsize=15);

    fig = Figure(size=(1050, 650),fontsize=15);
    Label(
        fig[0, 3],
        "Momentum budget snapshot, CGHR model",
        fontsize = 20,
        tellwidth = false
    )

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ((2 .* du) ./ 384),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\partial u/ \partial t_{\text{Euler}}"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    # ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # (adv_u + fv - detagdx + Fx + viscu + bottomdragu),
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title="Sum of budget terms"),
    # colorrange=(-1.5e-5, 1.5e-5)
    # );
    # hidedecorations!(ax2)
    # Colorbar(fig[1,4], hm2)

    # ax13, hm13 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # ((2 .* du) ./ 384) .- (adv_u + fv - detagdx + Fx + viscu + bottomdragu),
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title="Tendency - sum"),
    # colorrange=(-1.5e-13, 1.5e-13)
    # );
    # hidedecorations!(ax13)
    # Colorbar(fig[1,6], hm13)

    ax3, hm3 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (adv_u),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{adv}}^{u}"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax3)
    Colorbar(fig[1,4], hm3)

    ax4, hm4 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    fv,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{cor}}^{u}"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    Colorbar(fig[1,6], hm4)
    hidedecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -detagdx,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{press}}^{u}"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    # hidedecorations!(ax5)
    Colorbar(fig[2,2], hm5)

    # ax6, hm6 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # Fx,
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title="Wind stress"),
    # # colorrange=(-1.5e-5, 1.5e-5)
    # );
    # hideydecorations!(ax6)
    # Colorbar(fig[2,6], hm6)

    ax7, hm7 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    viscu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{visc}}^{u}"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    hideydecorations!(ax7)
    Colorbar(fig[2,4], hm7)

    ax8, hm8 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    bottomdragu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{BD}}^{u}"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    hideydecorations!(ax8)
    Colorbar(fig[2,6], hm8)

    ax8, hm8 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -detagdx + fv,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{press}}^{u} + G_{\text{cor}}^{u}"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hideydecorations!(ax8)
    Colorbar(fig[3,2], hm8)

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    # gg = fig[3, 1] = GridLayout()
    # gh = fig[3, 3] = GridLayout()

    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # same as above but time-averaged

    fig = Figure(size=(1200, 650),fontsize=15);
    Label(
        fig[0, 3],
        "3-year averaged momentum budget components, HRCG model",
        fontsize = 20,
        tellwidth = false
    )

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ((2 .* dutotal) ./ 384),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\partial u/ \partial t_{\text{Euler}}"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    ax3, hm3 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (adv_utotal),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{adv}}^{u}"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    hidedecorations!(ax3)
    Colorbar(fig[1,4], hm3)

    ax4, hm4 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    fvtotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{cor}}^{u}"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,6], hm4)
    hidedecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -detagdxtotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{press}}^{u}"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    # hidedecorations!(ax5)
    Colorbar(fig[2,2], hm5)

    ax7, hm7 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    viscutotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{visc}}^{u}"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    hideydecorations!(ax7)
    Colorbar(fig[2,4], hm7)

    ax8, hm8 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    bottomdragutotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{BD}}^{u}"),
    colorrange=(-1.5e-8, 1.5e-8)
    );
    hideydecorations!(ax8)
    Colorbar(fig[2,6], hm8)

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    # gg = fig[3, 1] = GridLayout()
    # gh = fig[3, 3] = GridLayout()

    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    fig = Figure(size=(1050, 700), fontsize=15);

    Label(
        fig[0, 3],
        "Momentum budget, v-components",
        fontsize = 20,
        tellwidth = false
    )

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ((2 .* dv) ./ 384),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Euler tendency"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (adv_v - fu - detagdy + viscv + bottomdragv),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Sum of budget terms"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax2)
    Colorbar(fig[1,4], hm2)

    ax13, hm13 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ((2 .* dv) ./ 384) .- (adv_v - fu - detagdy + viscv + bottomdragv),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Tendency - sum"),
    colorrange=(-1.5e-13, 1.5e-13)
    );
    hidedecorations!(ax13)
    Colorbar(fig[1,6], hm13)

    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (adv_v),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Advection"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax3)
    Colorbar(fig[2,2], hm3)

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -fu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Coriolis"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    Colorbar(fig[2,4], hm4)
    hidexdecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -detagdy,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Pressure gradient"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    hidedecorations!(ax5)
    Colorbar(fig[2,6], hm5)

    # ax6, hm6 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # Fx,
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title="Wind stress"),
    # # colorrange=(-1.5e-5, 1.5e-5)
    # );
    # hideydecorations!(ax6)
    # Colorbar(fig[2,6], hm6)

    ax7, hm7 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    viscv,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Viscosity"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    Colorbar(fig[3,2], hm7)

    ax8, hm8 = heatmap(fig[3,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    bottomdragv,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Bottom drag"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    hideydecorations!(ax8)
    Colorbar(fig[3,4], hm8)

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    gg = fig[3, 1] = GridLayout()
    gh = fig[3, 3] = GridLayout()

    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"], [ga, gb, gc, gd, ge, gf, gg, gh])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # KE budget components

    fig = Figure(size=(1200, 650),fontsize=15);
    Label(
        fig[0, 3],
        "Momentum budget, u-components",
        fontsize = 20,
        tellwidth = false
    )

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uhrcgall[:,:,n] .* ((2 .* du) ./ 384),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Euler tendency"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    # ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # (adv_u + fv - detagdx + Fx + viscu + bottomdragu),
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title="Sum of budget terms"),
    # colorrange=(-1.5e-5, 1.5e-5)
    # );
    # hidedecorations!(ax2)
    # Colorbar(fig[1,4], hm2)

    # ax13, hm13 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # ((2 .* du) ./ 384) .- (adv_u + fv - detagdx + Fx + viscu + bottomdragu),
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title="Tendency - sum"),
    # colorrange=(-1.5e-13, 1.5e-13)
    # );
    # hidedecorations!(ax13)
    # Colorbar(fig[1,6], hm13)

    ax3, hm3 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uhrcgall[:,:,n] .* (adv_u),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Advection"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax3)
    Colorbar(fig[1,4], hm3)

    ax4, hm4 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uhrcgall[:,:,n] .* fv,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Coriolis"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    Colorbar(fig[1,6], hm4)
    hidedecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uhrcgall[:,:,n] .* (-detagdx),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Pressure gradient"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    # hidedecorations!(ax5)
    Colorbar(fig[2,2], hm5)

    # ax6, hm6 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # Fx,
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title="Wind stress"),
    # # colorrange=(-1.5e-5, 1.5e-5)
    # );
    # hideydecorations!(ax6)
    # Colorbar(fig[2,6], hm6)

    ax7, hm7 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uhrcgall[:,:,n] .* viscu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Viscosity"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    Colorbar(fig[2,4], hm7)

    ax8, hm8 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uhrcgall[:,:,n] .* bottomdragu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Bottom drag"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    # hideydecorations!(ax8)
    Colorbar(fig[2,6], hm8)

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    # gg = fig[3, 1] = GridLayout()
    # gh = fig[3, 3] = GridLayout()

    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # time-averaged ke budget components

    fig = Figure(size=(1050, 650),fontsize=15);
    Label(
        fig[0, 3],
        "Time-averaged KE budget components",
        fontsize = 20,
        tellwidth = false
    )

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ((2 .* udutotal) ./ 384),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u\partial u/ \partial t_{\text{Euler}}"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    ax3, hm3 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uadv_utotal),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"uG_{\text{adv}}^{u}"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    hidedecorations!(ax3)
    Colorbar(fig[1,4], hm3)

    ax4, hm4 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ufvtotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"uG_{\text{cor}}^{u}"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,6], hm4)
    hidedecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -udetagdxtotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"uG_{\text{press}}^{u}"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    # hidedecorations!(ax5)
    Colorbar(fig[2,2], hm5)

    ax7, hm7 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uviscutotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"uG_{\text{visc}}^{u}"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    hideydecorations!(ax7)
    Colorbar(fig[2,4], hm7)

    ax8, hm8 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ubottomdragutotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"uG_{\text{BD}}^{u}"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    hideydecorations!(ax8)
    Colorbar(fig[2,6], hm8)

    ax8, hm8 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -udetagdxtotal + ufvtotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"uG_{\text{press}}^{u} + uG_{\text{cor}}^{u}"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    hideydecorations!(ax8)
    Colorbar(fig[3,2], hm8)

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()


    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], [ga, gb, gc, gd, ge, gf])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # rk4 computed tendency
    average_tend = sum(((2 .* durk4) ./ 384), dims=3)[:,:,1] ./ 1096
    average_energy = sum(uhrcgall[:,:,1:1096] .* ((2 .* durk4) ./ 384), dims=3)[:,:,1] ./ 1096

    fig = Figure(size=(900, 450),fontsize=15);
    Label(
        fig[0, 1:2],
        "Time-averaged RK4 tendency and KE budget",
        fontsize = 20,
        tellwidth = false
    )
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    average_tend,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\partial u / \partial t_{\text{RK4}}"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    # Colorbar(fig[1,2], hm1)

    ax1, hm1 = heatmap(fig[1,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    average_energy,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u\partial u / \partial t_{\text{RK4}}"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    Colorbar(fig[1,3], hm1)
    hideydecorations!(ax1)

end

function budget_term_checks()

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
        g=1e-20,
        H=500,
        cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
        ω=1e-22,
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

    u_, v_, eta_ = ShallowWaters.add_halo(uhrcgall[:,:,n], vhrcgall[:,:,n], etahrcgall[:,:,n], Slr)

    Seuler.Prog.u = u_
    Seuler.Prog.v = v_
    Seuler.Prog.η = eta_

    Sadv.Prog.u = u_
    Sadv.Prog.v = v_
    Sadv.Prog.η = eta_

    du, dv, adv_u, adv_v, fu, fv, detagdx, detagdy, bottomdragu, bottomdragv, viscu, viscv, Fx = compute_mom_budget(Seuler, Sadv, t, n);

    println("Check that initial condition doesn't propagate (should be zero): ", norm(Seuler.Prog.η .- eta_))

    # check volume flux components
    # these are used in the advection scheme to determine fu and fv
    # my thought process says that if coriolisoverh is correct, and the mass fluxes are correct,
    # then fu and fv should also be correct

    ShallowWaters.thickness!(Seuler.Diag.VolumeFluxes.h, eta_, Seuler.forcing.H)
    ShallowWaters.Ix!(Seuler.Diag.VolumeFluxes.h_u, Seuler.Diag.VolumeFluxes.h)
    ShallowWaters.Iy!(Seuler.Diag.VolumeFluxes.h_v, Seuler.Diag.VolumeFluxes.h)
    ShallowWaters.Ixy!(Seuler.Diag.Vorticity.h_q, Seuler.Diag.VolumeFluxes.h)

    println("Check that the sea-surface height h is correct (should be zero): ", norm(Seuler.Diag.VolumeFluxes.h .- (eta_ .+ Seuler.forcing.H)))

    U2 = (u_[2:end-1,2:end-1].*Seuler.Diag.VolumeFluxes.h_u).*Seuler.constants.scale_inv

    println("Check that the volume flux U is correct (should be zero): ", norm(Seuler.Diag.VolumeFluxes.U .- U2))

    # check the computed pressure gradient, should be identical to this
    detagdx2 = zeros(127,128)
    detagdy2 = zeros(128,127)

    ShallowWaters.∂x!(detagdx2, (Seuler.constants.g .* etahrcgall[:,:,n])./ Seuler.grid.Δ)
    ShallowWaters.∂y!(detagdy2, (Seuler.constants.g .* etahrcgall[:,:,n])./ Seuler.grid.Δ)

    println("Check that the pressure gradient is correct (should be zero): ", norm(detagdx .- detagdx2))

    fig = Figure(fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    detagdx .- detagdx2,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Pressure gradient check"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    fig

    # check the coriolis force computation -- from this we compute the terms fu and fv in the budget
    # theoretically, if this computation is correct then all computations to get fu, fv should also be correct
    m,n = size(Seuler.Diag.Vorticity.q)
    coriolisoverh = zeros(m,n)
    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            coriolisoverh[i,j] = Seuler.Diag.Vorticity.q[i,j] - (Seuler.Diag.Vorticity.dvdx[i+1,j+1] - Seuler.Diag.Vorticity.dudy[i+1,j+1]) / Seuler.Diag.Vorticity.h_q[i,j]
        end
    end
    coriolisoverh2 = Seuler.grid.f_q ./ S.Diag.Vorticity.h_q

    fig = Figure(fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    coriolisoverh .- coriolisoverh2,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="f / h check"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    fig = Figure(fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ((2 .* du) ./ 384) .- (adv_u + fv - detagdx2 + Fx + viscu + bottomdragu),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Euler tendency"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

end

function momentum_budget_multi2()

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
        nn_forcing_momentum=false,
        nn_forcing_dissipation=true,
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
        nn_forcing_momentum=false,
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=2,
        initial_cond="rest",
        # initpath = "./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves"
    );
    Slr2 = ShallowWaters.model_setup(Plr2);

    Seuler = deepcopy(Slr);
    current = 1
    for m in (Seuler.Diag.CNNVars.model_Su, Seuler.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-6-8-10-13-15-18-23-28-33-38-41-44-48-51-53-58-63-65-68-73-78-83-86-88daystart_2dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    Sadv = deepcopy(Slr2);
    current = 1
    for m in (Sadv.Diag.CNNVars.model_Su, Sadv.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-6-8-10-13-15-18-23-28-33-38-41-44-48-51-53-58-63-65-68-73-78-83-86-88daystart_2dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    t = 225 * Slr.grid.dtint
    n = 1096

    u_, v_, eta_ = ShallowWaters.add_halo(umulti2[:,:,n], vmulti2[:,:,n], etamulti2[:,:,n], Slr)

    Seuler.Prog.u = u_
    Seuler.Prog.v = v_
    Seuler.Prog.η = eta_

    Sadv.Prog.u = u_
    Sadv.Prog.v = v_
    Sadv.Prog.η = eta_

    du, dv, adv_u, adv_v, fu, fv, detagdx, detagdy, bottomdragu, bottomdragv, viscu, viscv, Fx = compute_mom_budget(Seuler, Sadv, t, n);
    dutotal, dvtotal, adv_utotal, adv_vtotal, futotal, fvtotal, detagdxtotal, detagdytotal, bottomdragutotal, bottomdragvtotal, viscutotal, viscvtotal, Fxtotal, Sutotal, Svtotal = compute_timeaveraged_momentumbudget(umulti2, vmulti2, etamulti2, Seuler, Sadv);
    udutotal, vdvtotal, uadv_utotal, vadv_vtotal, vfutotal, ufvtotal, udetagdxtotal, vdetagdytotal, ubottomdragutotal, vbottomdragvtotal, uviscutotal, vviscvtotal, uFxtotal, uSutotal, vSvtotal = compute_timeaveraged_energybudget(umulti2, vmulti2, etamulti2, Seuler, Sadv);

    # all of the rhs terms are computed with a *scaled* prognostic variable, so to actually 
    # see the momentum budget terms we need to remove that scaling. The division by Delta is because
    # all of the timestepping is done with dt / Delta, and we want the term that gets multiplied by dt
    # I should have corrected for the scalings in the values returned by compute_mom_budget already,
    # but for the stuff I'm pulling after (see below) I need to do it manually

    utendency = ((2 .* dutotal) ./ Seuler.grid.dtint)
    parameterization = Seuler.Diag.CNNVars.S_u ./ (Seuler.grid.Δ * Seuler.constants.scale);
    budgetsum = (adv_utotal + fvtotal - detagdxtotal + Fxtotal + viscutotal + bottomdragutotal + Sutotal);

    fig = Figure(fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    utendency .- budgetsum,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Euler tendency - budget sum"),
    colorrange=(-5e-14, 5e-14)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    fig = Figure(size=(1050, 700), fontsize=15);

    Label(
        fig[0, 3],
        "Momentum budget snapshot, ensemble 2 day parameterization",
        fontsize = 25,
        tellwidth = false
    )

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    2 .* du ./ Seuler.grid.dtint,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\partial u / \partial t_{\text{Euler}}"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    ax3, hm3 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (adv_u),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{adv}}^{u}"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax3)
    Colorbar(fig[1,4], hm3)

    ax4, hm4 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    fv,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{cor}}^{u}"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    Colorbar(fig[1,6], hm4)
    hidedecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -detagdx,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{press}}^{u}"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    hidexdecorations!(ax5)
    Colorbar(fig[2,2], hm5)

    # ax6, hm6 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # Fx,
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title="Wind stress"),
    # # colorrange=(-1.5e-5, 1.5e-5)
    # );
    # hideydecorations!(ax6)
    # Colorbar(fig[2,6], hm6)

    ax7, hm7 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    viscu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{visc}}^{u}"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    hideydecorations!(ax7)
    Colorbar(fig[2,4], hm7)

    ax8, hm8 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    bottomdragu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{BD}}^{u}"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    hideydecorations!(ax8)
    Colorbar(fig[2,6], hm8)

    ax9, hm9 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Seuler.Diag.CNNVars.S_u ./ (Seuler.grid.Δ * Seuler.constants.scale),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"S_u"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[3,2], hm9)

    ax9, hm9 = heatmap(fig[3,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -detagdxtotal + fvtotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{press}}^{u} + G_{\text{cor}}^{u}"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[3,4], hm9)

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    gg = fig[3, 1] = GridLayout()
    # gh = fig[3, 3] = GridLayout()

    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"], [ga, gb, gc, gd, ge, gf, gg])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # same as above but time-averaged

    fig = Figure(size=(1050, 650),fontsize=15);
    Label(
        fig[0, 3],
        "3-year averaged momentum budget components, ensemble 2 day model",
        fontsize = 20,
        tellwidth = false
    )

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ((2 .* dutotal) ./ 384),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\partial u/ \partial t_{\text{Euler}}"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    ax3, hm3 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (adv_utotal),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{adv}}^{u}"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    hidedecorations!(ax3)
    Colorbar(fig[1,4], hm3)

    ax4, hm4 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    fvtotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{cor}}^{u}"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,6], hm4)
    hidedecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -detagdxtotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{press}}^{u}"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    # hidedecorations!(ax5)
    Colorbar(fig[2,2], hm5)

    ax7, hm7 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    viscutotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{visc}}^{u}"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    hideydecorations!(ax7)
    Colorbar(fig[2,4], hm7)

    ax8, hm8 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    bottomdragutotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"G_{\text{BD}}^{u}"),
    colorrange=(-1.5e-8, 1.5e-8)
    );
    hideydecorations!(ax8)
    Colorbar(fig[2,6], hm8)

    ax9, hm9 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Sutotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"S_u"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[3,2], hm9)

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    gg = fig[3, 1] = GridLayout()
    # gh = fig[3, 3] = GridLayout()

    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"], [ga, gb, gc, gd, ge, gf, gg])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # time-averaged ke budget components

    fig = Figure(size=(1000, 650),fontsize=15);
    Label(
        fig[0, 3],
        "Time-averaged KE budget components, ensemble 2 day model",
        fontsize = 20,
        tellwidth = false
    )

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ((2 .* udutotal) ./ 384),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"u\partial u/ \partial t_{\text{Euler}}"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    ax3, hm3 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (uadv_utotal),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"uG_{\text{adv}}^{u}"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    hidedecorations!(ax3)
    Colorbar(fig[1,4], hm3)

    ax4, hm4 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ufvtotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"uG_{\text{cor}}^{u}"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,6], hm4)
    hidedecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -udetagdxtotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"uG_{\text{press}}^{u}"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    # hidedecorations!(ax5)
    Colorbar(fig[2,2], hm5)

    ax7, hm7 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uviscutotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"uG_{\text{visc}}^{u}"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    hideydecorations!(ax7)
    Colorbar(fig[2,4], hm7)

    ax8, hm8 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    ubottomdragutotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"uG_{\text{BD}}^{u}"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    hideydecorations!(ax8)
    Colorbar(fig[2,6], hm8)

    ax9, hm9 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    uSutotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"S_u"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    Colorbar(fig[3,2], hm9)

    ax9, hm9 = heatmap(fig[3,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -udetagdxtotal + ufvtotal,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"uG_{\text{press}}^{u} + uG_{\text{cor}}^{u}"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    Colorbar(fig[3,4], hm9)

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    gc = fig[1, 5] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 3] = GridLayout()
    gf = fig[2, 5] = GridLayout()
    gg = fig[3, 1] = GridLayout()
    # gh = fig[3, 3] = GridLayout()

    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"], [ga, gb, gc, gd, ge, gf, gg])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

end

function momentum_budget_zannabolton()
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
        zb_forcing_dissipation=true,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=false,
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
        zb_forcing_dissipation=true,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=false,
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

    u_, v_, eta_ = ShallowWaters.add_halo(uzb[:,:,n], vzb[:,:,n], etazb[:,:,n], Slr)

    Seuler.Prog.u = u_
    Seuler.Prog.v = v_
    Seuler.Prog.η = eta_

    Sadv.Prog.u = u_
    Sadv.Prog.v = v_
    Sadv.Prog.η = eta_

    du, dv, adv_u, adv_v, fu, fv, detagdx, detagdy, bottomdragu, bottomdragv, viscu, viscv, Fx = compute_mom_budget(Seuler, Sadv, t, n);

    # detagdx2 = zeros(127,128)
    # detagdy2 = zeros(128,127)

    # ShallowWaters.∂x!(detagdx2, (Seuler.constants.g .* etahrcgall[:,:,n])./ Seuler.grid.Δ)
    # ShallowWaters.∂y!(detagdy2, (Seuler.constants.g .* etahrcgall[:,:,n])./ Seuler.grid.Δ)

    # all of the rhs terms are computed with a *scaled* prognostic variable, so to actually 
    # see the momentum budget terms we need to remove that scaling. The division by Delta is because
    # all of the timestepping is done with dt / Delta, and we want the term that gets multiplied by dt
    # I should have corrected for the scalings in the values returned by compute_mom_budget already,
    # but for the stuff I'm pulling after (see below) I need to do it manually

    eulertendency = ((2 .* du) ./ Seuler.grid.dtint);
    parameterization = Seuler.Diag.ZBVars.S_u ./ (Seuler.grid.Δ * Seuler.constants.scale);
    budgetsum = (adv_u + fv - detagdx + Fx + viscu + bottomdragu + 2 .* parameterization);

    fig = Figure(size=(1050, 700), fontsize=15);

    Label(
        fig[0, 3],
        "Online u budget, ZB20 parameterization",
        fontsize = 20,
        tellwidth = false
    )

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    eulertendency,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Euler tendency"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    budgetsum,
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
    hidedecorations!(ax3)
    Colorbar(fig[1,6], hm3)

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    fv,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Coriolis"),
    colorrange=(-1.5e-4, 1.5e-4)
    );
    Colorbar(fig[2,2], hm4)
    hidexdecorations!(ax4)

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
    Colorbar(fig[3,2], hm7)

    ax8, hm8 = heatmap(fig[3,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    bottomdragu,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Bottom drag"),
    colorrange=(-1.5e-7, 1.5e-7)
    );
    hideydecorations!(ax8)
    Colorbar(fig[3,4], hm8)

    ax9, hm9 = heatmap(fig[3,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    parameterization,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Parameterization"),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hideydecorations!(ax9)
    Colorbar(fig[3,6], hm9)

end