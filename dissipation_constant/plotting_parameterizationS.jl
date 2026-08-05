function create_models()

    Scghr = deepcopy(S);
    S5 = deepcopy(S);
    S10 = deepcopy(S);
    current = 1
    for m in (S10.Diag.CNNVars.model_Su, S10.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_states_10dayoptimization_startfrom5daystate_noeta_gelu_30iterations.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end
    S20 = deepcopy(S);
    current = 1
    for m in (S20.Diag.CNNVars.model_Su, S20.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_states_20dayoptimization_startfrom10daystate_noeta_10iterations.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end
    S30 = deepcopy(S);
    current = 1
    for m in (S30.Diag.CNNVars.model_Su, S30.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/states_noetainloss/result_online_state_30dayoptimzation_startfrom30day6iterations_constantdissipation_15iterations_21totaliterations_8hourdata_200maxhistory_fixedcfl.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    Smulti2 = deepcopy(S);
    current = 1
    for m in (Smulti2.Diag.CNNVars.model_Su, Smulti2.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-6-8-10-13-15-18-23-28-33-38-41-44-48-51-53-58-63-65-68-73-78-83-86-88daystart_2dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    Smulti3 = deepcopy(S);
    current = 1
    for m in (Smulti3.Diag.CNNVars.model_Su, Smulti3.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86daystart_3dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    Smulti10 = deepcopy(S);
    current = 1
    for m in (Smulti10.Diag.CNNVars.model_Su, Smulti10.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/result_multistate_5-20-35-50-65-75daystart_10dayoptimization_initialweights20daystate_fixedcfl_15iterations_constdissipation.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    Smulti20 = deepcopy(S);
    current = 1
    for m in (Smulti20.Diag.CNNVars.model_Su, Smulti20.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/result_multistate_5-25-45-65daystart_20dayoptimization_initialweights20daystate_fixedcfl_15iterations_constdissipation.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    Soffline = deepcopy(S);
    current = 1
    for m in (Soffline.Diag.CNNVars.model_Su, Soffline.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(load_object("./dissipation_constant/tuned_weights/result_offline_150iterations_geluactivation_111925.jld2").solution[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    Szb = deepcopy(S);
    Scghr = deepcopy(S);


end

function online_S_first3years()

    # online plot within the first three years
    t = 1096

    uzb_, vzb_, etazb_ = ShallowWaters.add_halo(Float64.(uzb[:,:,t]), Float64.(vzb[:,:,t]), Float64.(etazb[:,:,t]), zeros(128,128), S);
    u10s_, v10s_, eta10s_ = ShallowWaters.add_halo(Float64.(u10s[:,:,t]), Float64.(v10s[:,:,t]), Float64.(eta10s[:,:,t]), zeros(128,128), S);
    u20s_, v20s_, eta20s_ = ShallowWaters.add_halo(Float64.(u20s[:,:,t]), Float64.(v20s[:,:,t]), Float64.(eta20s[:,:,t]), zeros(128,128), S);
    u30s_, v30s_, eta30s_ = ShallowWaters.add_halo(Float64.(u30s[:,:,t]), Float64.(v30s[:,:,t]), Float64.(eta30s[:,:,t]), zeros(128,128), S);
    # umulti1_, vmulti1_, etamulti1_ = ShallowWaters.add_halo(Float64.(umulti1[:,:,t]), Float64.(vmulti1[:,:,t]), Float64.(etamulti1[:,:,t]), zeros(128,128), S);
    # umulti1more_, vmulti1more_, etamulti1more_ = ShallowWaters.add_halo(Float64.(umulti1more[:,:,t]), Float64.(vmulti1more[:,:,t]), Float64.(etamulti1more[:,:,t]), zeros(128,128), S);
    umulti2_, vmulti2_, etamulti2_ = ShallowWaters.add_halo(Float64.(umulti2[:,:,t]), Float64.(vmulti2[:,:,t]), Float64.(etamulti2[:,:,t]), zeros(128,128), S);
    # umulti3_, vmulti3_, etamulti3_ = ShallowWaters.add_halo(Float64.(umulti3[:,:,t]), Float64.(vmulti3[:,:,t]), Float64.(etamulti3[:,:,t]), zeros(128,128), S);
    umulti3more_, vmulti3more_, etamulti3more_ = ShallowWaters.add_halo(Float64.(umulti3more[:,:,t]), Float64.(vmulti3more[:,:,t]), Float64.(etamulti3more[:,:,t]), zeros(128,128), S);
    umulti10_, vmulti10_, etamulti10_ = ShallowWaters.add_halo(Float64.(umulti10[:,:,t]), Float64.(vmulti10[:,:,t]), Float64.(etamulti10[:,:,t]), zeros(128,128), S);
    umulti20_, vmulti20_, etamulti20_ = ShallowWaters.add_halo(Float64.(umulti20[:,:,t]), Float64.(vmulti20[:,:,t]), Float64.(etamulti20[:,:,t]), zeros(128,128), S);

    ShallowWaters.ZB_momentum(Float64.(uzb_), Float64.(vzb_), Szb, Szb.Diag);
    ShallowWaters.CNN_momentum(u20s_, v20s_, S20);
    ShallowWaters.CNN_momentum(u30s_, v30s_, S30);
    ShallowWaters.CNN_momentum(u10s_, v10s_, S10);

    ShallowWaters.CNN_momentum(umulti2_, vmulti2_, Smulti2);
    ShallowWaters.CNN_momentum(umulti3more_, vmulti3more_, Smulti3);
    ShallowWaters.CNN_momentum(umulti10_, vmulti10_, Smulti10);
    ShallowWaters.CNN_momentum(umulti20_, vmulti20_, Smulti20);

    fig = Figure(size=(1040, 520), fontsize=15);

    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Suhr[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\overline{S_u}(3 \; \text{years}, x, y)"),
    colorrange=(-maximum(abs.(Suhr[:,:,t])),maximum(abs.(Suhr[:,:,t])))
    );
    Colorbar(fig[1,2], hm1, label="1/s")
 
    ax1, hm1 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Szb.Diag.ZBVars.S_u ./ Szb.grid.Δ,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="ZB20"),
    colorrange=(-maximum(abs.(Suhr[:,:,t])),maximum(abs.(Suhr[:,:,t])))
    );
    Colorbar(fig[1,4], hm1, label="1/s")

    ax2, hm2 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti2.Diag.CNNVars.S_u ./ Smulti2.grid.Δ,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 2 day"),
    colorrange=(-maximum(abs.(Suhr[:,:,t])),maximum(abs.(Suhr[:,:,t])))
    );
    Colorbar(fig[1,6], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti3.Diag.CNNVars.S_u ./ Smulti3.grid.Δ,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 3 day"),
    colorrange=(-maximum(abs.(Suhr[:,:,t])),maximum(abs.(Suhr[:,:,t])))
    );
    Colorbar(fig[2,2], hm1, label="1/s")

    ax4, hm4 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti10.Diag.CNNVars.S_u./ Smulti10.grid.Δ,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 10 day"),
    colorrange=(-maximum(abs.(Suhr[:,:,t])),maximum(abs.(Suhr[:,:,t])))
    );
    Colorbar(fig[2,4], hm1, label="1/s")

    ax3, hm3 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    S30.Diag.CNNVars.S_u ./ S30.grid.Δ,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="30 day"),
    colorrange=(-maximum(abs.(Suhr[:,:,t])),maximum(abs.(Suhr[:,:,t])))
    );
    Colorbar(fig[2,6], hm1, label="1/s")

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

end

function offline_S_firstthreeyears()

    # offline plot within the first three years

    t = 1096
    uhrcg_, vhrcg_, etahrcg_ = ShallowWaters.add_halo(Float64.(uhrcgall[:,:,t]), Float64.(vhrcgall[:,:,t]), Float64.(etahrcgall[:,:,t]), zeros(128,128), S);

    ShallowWaters.ZB_momentum(uhrcg_, vhrcg_, Szb, Szb.Diag);
    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, S20);
    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, S30);
    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, S10);

    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, Smulti2);
    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, Smulti3);
    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, Smulti10);
    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, Smulti20);
    ShallowWaters.CNN_momentum(uhrcg_, vhrcg_, Soffline);

    # quick comparison

    fig = Figure(fontsize=15, size = (900, 300));

    ax00, hm00 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (advec_hrcg[1][:,:,t] .- advec_cg[1][:,:,t]),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"(S_{\text{adv}})_u"),
    colorrange=(-1e-5, 1e-5)
    );

    ax0, hm0 = heatmap(fig[1,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Szb.Diag.ZBVars.S_u ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"(S_{\text{ZB}})_u"),
    colorrange=(-1e-5, 1e-5)
    );
    hideydecorations!(ax0)

    ax1, hm1 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (advec_hrcg[1][:,:,t] .- advec_cg[1][:,:,t]) .- Szb.Diag.ZBVars.S_u ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"(S_{\text{adv}})_u - (S_{\text{ZB}})_u"),
    colorrange=(-1e-5, 1e-5)
    );
    Colorbar(fig[1,4], hm1)
    hideydecorations!(ax1)

    # S_u
    fig = Figure(size=(900, 780), fontsize=15);

    Label(
        fig[0, 2],
        L"S_u(3 \text{ years}, x, y)",
        fontsize = 20,
        tellwidth = false
    )

    t = 1096
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (tend_rk4_hrcg[1][:,:,t]./48 .- tend_rk4_cg[1][:,:,t]./384),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Total tendencies, RK4"),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j]))),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidexdecorations!(ax1)

    ax0, hm0 = heatmap(fig[1,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    2 .* (tend_euler_hrcg[1][:,:,t]./48 .- tend_euler_cg[1][:,:,t]./384),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Total tendencies, Euler"),
    # colorrange=(-maximum(abs.(Suadvec[:,:,t])),maximum(abs.(Suadvec[:,:,t]))),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax0)

    ax00, hm00 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (advec_hrcg[1][:,:,t] .- advec_cg[1][:,:,t]),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Nonlinear advection, Euler"),
    # colorrange=(-maximum(abs.(Suadvec[:,:,t])),maximum(abs.(Suadvec[:,:,t]))),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax00)
    # Colorbar(fig[1,2], hm0)

    ax2, hm2 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    S.grid.Δ .* Suapprox[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Nonlinear advection approx."),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j])))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidexdecorations!(ax2)


    s = Szb.grid.Δ * Szb.grid.scale
    ax3, hm3 = heatmap(fig[2,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Szb.Diag.ZBVars.S_u ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="ZB20"),
    # colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u./ s)),maximum(abs.(Szb.Diag.ZBVars.S_u./ s)))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    # Colorbar(fig[1,2], hm2, label=L"m/s^2")
    hidedecorations!(ax3)

    # ax4, hm4 = heatmap(fig[2,2], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # Soffline.Diag.CNNVars.S_u ./ s,
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title="Offline-learned NN"),
    # # colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u./ s)),maximum(abs.(Szb.Diag.ZBVars.S_u./ s)))
    # colorrange=(-1.5e-5, 1.5e-5)
    # );
    # # Colorbar(fig[1,4], hm2, label=L"m/s^2")
    # hidedecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti2.Diag.CNNVars.S_u ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 2 day"),
    # colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u./ s)),maximum(abs.(Szb.Diag.ZBVars.S_u./ s)))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    # Colorbar(fig[1,6], hm2, label=L"m/s^2")
    hidedecorations!(ax5)

    ax6, hm6 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti3.Diag.CNNVars.S_u ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 3 day"),
    # colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u./ s)),maximum(abs.(Szb.Diag.ZBVars.S_u./ s)))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    # Colorbar(fig[2,2], hm2, label=L"m/s^2")

    ax7, hm7 = heatmap(fig[3,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti10.Diag.CNNVars.S_u ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 10 day"),
    # colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u./ s)),maximum(abs.(Szb.Diag.ZBVars.S_u./ s)))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    # Colorbar(fig[2,4], hm2, label=L"m/s^2")
    hideydecorations!(ax7)

    ax8, hm8 = heatmap(fig[3,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    S30.Diag.CNNVars.S_u ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="30 day"),
    # colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u./ s)),maximum(abs.(Szb.Diag.ZBVars.S_u./ s)))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    # Colorbar(fig[2,6], hm2, label=L"m/s^2")
    hideydecorations!(ax8)

    Colorbar(fig[1:3,4], hm1, label=L"m/s^2")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 2] = GridLayout()
    gc = fig[1, 3] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 2] = GridLayout()
    gf = fig[2, 3] = GridLayout()
    gg = fig[3, 1] = GridLayout()
    gh = fig[3, 2] = GridLayout()
    gi = fig[3, 3] = GridLayout()

    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"], [ga, gb, gc, gd, ge, gf, gg, gh, gi])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # S_v
    fig = Figure(size=(900, 780), fontsize=15);

    Label(
        fig[0, 2],
        L"S_v(3 \text{ years}, x, y)",
        fontsize = 20,
        tellwidth = false
    )

    t = 1096
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (tend_rk4_hrcg[2][:,:,t]./48 .- tend_rk4_cg[2][:,:,t]./384),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Total tendencies, RK4"),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j]))),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidexdecorations!(ax1)

    ax0, hm0 = heatmap(fig[1,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    2 .* (tend_euler_hrcg[2][:,:,t]./48 .- tend_euler_cg[2][:,:,t]./384),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Total tendencies, Euler"),
    # colorrange=(-maximum(abs.(Suadvec[:,:,t])),maximum(abs.(Suadvec[:,:,t]))),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax0)

    ax00, hm00 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (advec_hrcg[2][:,:,t] .- advec_cg[2][:,:,t]),
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Nonlinear advection, Euler"),
    # colorrange=(-maximum(abs.(Suadvec[:,:,t])),maximum(abs.(Suadvec[:,:,t]))),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax00)
    # Colorbar(fig[1,2], hm0)

    ax2, hm2 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    S.grid.Δ .* Svapprox[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Nonlinear advection approx."),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j])))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidexdecorations!(ax2)


    s = Szb.grid.Δ * Szb.grid.scale
    ax3, hm3 = heatmap(fig[2,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Szb.Diag.ZBVars.S_v ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="ZB20"),
    # colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u./ s)),maximum(abs.(Szb.Diag.ZBVars.S_u./ s)))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    # Colorbar(fig[1,2], hm2, label=L"m/s^2")
    hidedecorations!(ax3)

    # ax4, hm4 = heatmap(fig[2,2], LinRange(0, 3840, 128),
    # LinRange(0, 3840, 128),
    # Soffline.Diag.CNNVars.S_u ./ s,
    # colormap=:balance,
    # axis=(xlabel="km", ylabel="km", title="Offline-learned NN"),
    # # colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u./ s)),maximum(abs.(Szb.Diag.ZBVars.S_u./ s)))
    # colorrange=(-1.5e-5, 1.5e-5)
    # );
    # # Colorbar(fig[1,4], hm2, label=L"m/s^2")
    # hidedecorations!(ax4)

    ax5, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti2.Diag.CNNVars.S_v ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 2 day"),
    # colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u./ s)),maximum(abs.(Szb.Diag.ZBVars.S_u./ s)))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    # Colorbar(fig[1,6], hm2, label=L"m/s^2")
    hidedecorations!(ax5)

    ax6, hm6 = heatmap(fig[3,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti3.Diag.CNNVars.S_v ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 3 day"),
    # colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u./ s)),maximum(abs.(Szb.Diag.ZBVars.S_u./ s)))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    # Colorbar(fig[2,2], hm2, label=L"m/s^2")

    ax7, hm7 = heatmap(fig[3,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Smulti10.Diag.CNNVars.S_v ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Ensemble 10 day"),
    # colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u./ s)),maximum(abs.(Szb.Diag.ZBVars.S_u./ s)))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    # Colorbar(fig[2,4], hm2, label=L"m/s^2")
    hideydecorations!(ax7)

    ax8, hm8 = heatmap(fig[3,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    S30.Diag.CNNVars.S_v ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="30 day"),
    # colorrange=(-maximum(abs.(Szb.Diag.ZBVars.S_u./ s)),maximum(abs.(Szb.Diag.ZBVars.S_u./ s)))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    # Colorbar(fig[2,6], hm2, label=L"m/s^2")
    hideydecorations!(ax8)

    Colorbar(fig[1:3,4], hm1, label=L"m/s^2")

    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 2] = GridLayout()
    gc = fig[1, 3] = GridLayout()
    gd = fig[2, 1] = GridLayout()
    ge = fig[2, 2] = GridLayout()
    gf = fig[2, 3] = GridLayout()
    gg = fig[3, 1] = GridLayout()
    gh = fig[3, 2] = GridLayout()
    gi = fig[3, 3] = GridLayout()

    for (label, layout) in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"], [ga, gb, gc, gd, ge, gf, gg, gh, gi])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

end

function extras()

    # viscosity and bottom drag SGS contributions to S_tot, also looking at the nonlinear advection

    fig = Figure(size=(750, 325), fontsize=15);
    j = 1096
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (Mu_hrdownsized[:,:,j] .- visc_cg[1][:,:,j]) ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\overline{M_{MM}(3 \text{ years}, u)} - M_{MM}(3 \text{ years}, \overline{u})"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    Colorbar(fig[1,2], hm1, label=L"m/s^2")

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (Bu_hrdownsized[:,:,j] .- bd_cg[1][:,:,j]) ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\overline{M_{BF}(3 \text{ years}, u)} - M_{BF}(3 \text{ years}, \overline{u})"),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j])))
    colorrange=(-1.5e-7, 1.5e-7)
    );
    Colorbar(fig[1,4], hm2, label=L"m/s^2")
    hideydecorations!(ax2)
    
    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 3] = GridLayout()
    for (label, layout) in zip(["(a)", "(b)"], [ga, gb])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = 15,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right)
    end

    # comparing "true" Su fields
    s = S.grid.Δ * S.grid.scale

    fig = Figure(size=(950, 500), fontsize=15);

    Label(
        fig[0, 3],
        L"S_u(3 \text{ years}, x, y)",
        fontsize = 20,
        tellwidth = false
    )

    t = 1096
    euler =  2. * (tendu_hrdownsized[:,:,t]./48 - tend_euler_cg[1][:,:,t]./384)
    advec = (Advecu_hrdownsized[:,:,t]./ (3.75e3) - advec_cg[1][:,:,t]./(30e3))
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -Suhr[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Total tendencies, RK4"),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j])))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (2 .* tend_euler_hrcg[1][:,:,t]) ./ 48,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Euler tendency"),
    # colorrange=(-maximum(abs.((Advecu_hrdownsized[:,:,j] .- advec_cghr[1][:,:,j]) ./ S.grid.Δ)),maximum(abs.((Advecu_hrdownsized[:,:,j] .- advec_cghr[1][:,:,j]) ./ S.grid.Δ)))
    # colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax2)
    Colorbar(fig[1,4], hm2)

    ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    Advecu_hrdownsized[:,:,t] ./ 3.75e3,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Advection term"),
    # colorrange=(-maximum(abs.(- Advecu_hrdownsized[:,:,t]./(64)  + advec_cg[1][:,:,t]./(30e3 * 64))), maximum(abs.(- Advecu_hrdownsized[:,:,t]./(3.75e3 * 64)  + advec_cg[1][:,:,t]./(64))))
    colorrange=(-2e-5, 2e-5)
    );
    hidedecorations!(ax3)
    Colorbar(fig[1,6], hm3)

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    advec,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Nonlinear advection"),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j])))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[2,2], hm4)

    visc = visc_hrcg[1][:,:,t] ./ (64 * 3.75e3) - visc_cg[1][:,:,t] ./ (64 * 30e3)
    ax5, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    visc,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Viscosity"),
    # colorrange=(-maximum(abs.(visc)),maximum(abs.(visc)))
    colorrange=(-1e-6, 1e-6)
    );
    Colorbar(fig[2,4], hm5)

    bd = bd_hrcg[1][:,:,t] ./ (64 * 3.75e3) -  bd_cg[1][:,:,t]./ (64 * 30e3)
    ax6, hm6 = heatmap(fig[2,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    bd,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Bottom drag"),
    colorrange=(-maximum(abs.(bd)),maximum(abs.(bd)))
    # colorrange=(-1.5e-7, 1.5e-7)
    );
    Colorbar(fig[2,6], hm6)

    # Su difference fields
    fig = Figure(size=(950, 500), fontsize=15);

    Label(
        fig[0, 3],
        L"S_u(3 \text{ years}, x, y)",
        fontsize = 20,
        tellwidth = false
    )

    t = 1096
    euler = 2 .* (tend_euler_hrcg[1][:,:,t]./48 - tend_euler_cg[1][:,:,t]./384)
    advec = adv_uhrcg .- adv_ulr
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -Suhr[:,:,j] .- euler,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Total tendencies, |RK4 - Euler|"),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j])))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    abs.(-Suhr[:,:,j] - advec),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="|RK4 - nonlinear advection|"),
    # colorrange=(-maximum(abs.((Advecu_hrdownsized[:,:,j] .- advec_cghr[1][:,:,j]) ./ S.grid.Δ)),maximum(abs.((Advecu_hrdownsized[:,:,j] .- advec_cghr[1][:,:,j]) ./ S.grid.Δ)))
    colorrange=(0, 4e-5)
    );
    hidedecorations!(ax2)
    Colorbar(fig[1,4], hm2)

    ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    euler .- advec,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="|Euler - advec|"),
    # colorrange=(-maximum(abs.((Advecu_hrdownsized[:,:,j] .- advec_cg[1][:,:,j]) ./ S.grid.Δ)),maximum(abs.((Advecu_hrdownsized[:,:,j] .- advec_cg[1][:,:,j]) ./ S.grid.Δ)))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax3)
    Colorbar(fig[1,6], hm3)

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    euler .- 2 .* advec,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="|Euler - 2 * nonlinear advection|"),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j])))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[2,2], hm4)

    ax5, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    abs.(euler .- S.grid.Δ .* Suapprox[:,:,j]),
    colormap=:amp,
    axis=(xlabel="km", ylabel="km", title="|Euler - nonlinear advection approx.|"),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j])))
    colorrange=(0, 4e-5)
    );
    hideydecorations!(ax5)
    Colorbar(fig[2,4], hm5)

    # comparing "true" Sv fields
    s = S.grid.Δ * S.grid.scale

    fig = Figure(size=(800, 400), fontsize=15);

    Label(
        fig[0, 3],
        L"S_v(3 \text{ years}, x, y)",
        fontsize = 20,
        tellwidth = false
    )

    j = 1096
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -Svhr[:,:,j],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Total tendencies"),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j]))),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    Colorbar(fig[1,2], hm1)
    hidexdecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -(Advecv_hrdownsized[:,:,j] .- advec_cghr[2][:,:,j]) ./ S.grid.Δ,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Nonlinear advection"),
    # colorrange=(-maximum(abs.((Advecu_hrdownsized[:,:,j] .- advec_cghr[1][:,:,j]) ./ S.grid.Δ)),maximum(abs.((Advecu_hrdownsized[:,:,j] .- advec_cghr[1][:,:,j]) ./ S.grid.Δ)))
    colorrange=(-6e-5, 6e-5)
    );
    hidedecorations!(ax2)
    Colorbar(fig[1,4], hm2)

    ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    S.grid.Δ .* Svapprox[:,:,j],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Nonlinear advection approx."),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j])))
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidedecorations!(ax3)
    Colorbar(fig[1,6], hm3)

    ax4, hm4 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (Bv_hrdownsized[:,:,j] .- bd_cghr[2][:,:,j]) ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\overline{M_{BF}(3 \text{ years}, u)} - M_{BF}(3 \text{ years}, \overline{u})"),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j])))
    colorrange=(-1.5e-7, 1.5e-7)
    );
    Colorbar(fig[2,2], hm4)

    ax5, hm5 = heatmap(fig[2,3], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    (Mv_hrdownsized[:,:,j] .- visc_hrcg[2][:,:,j]) ./ s,
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title=L"\overline{M_{MM}(3 \text{ years}, u)} - M_{MM}(3 \text{ years}, \overline{u})"),
    colorrange=(-1.5e-6, 1.5e-6)
    );
    Colorbar(fig[2,4], hm5)
    hideydecorations!(ax5)

    # windowed velocity stuff


    Subefore = load_object("./dissipation_constant/ke_transfers/trueS_fromwindowedtendencies_withRK4_SuSv_first3years_dailysaves_S(overline(window(u))).jld2")[1];
    Subefore_noeta = load_object("./dissipation_constant/ke_transfers/trueS_fromwindowedtendencies_withRK4_SuSv_first3years_dailysaves_S(overline(window(u)))_nowindowedeta.jld2")[1];
    
    fig = Figure(size=(800, 800), fontsize=15);
    Label(
        fig[0, 2],
        L"S_u(3 \text{ years}, x, y)",
        fontsize = 20,
        tellwidth = false
    )

    t = 1096
    ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -Suhr[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Total tendencies, no window"),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j]))),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hidexdecorations!(ax1)

    ax2, hm2 = heatmap(fig[1,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -winubefore .* Suhr[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Total tendencies, window after"),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j]))),
    colorrange=(-1.5e-5, 1.5e-5)
    );
    hideydecorations!(ax1)

    ax3, hm3 = heatmap(fig[2,1], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -Subefore[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Total tendencies, window before"),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j]))),
    colorrange=(-1.5e-5, 1.5e-5)
    );

    ax4, hm4 = heatmap(fig[2,2], LinRange(0, 3840, 128),
    LinRange(0, 3840, 128),
    -Subefore_noeta[:,:,t],
    colormap=:balance,
    axis=(xlabel="km", ylabel="km", title="Total tendencies, window before, did not window eta"),
    # colorrange=(-maximum(abs.(Suhr[:,:,j])),maximum(abs.(Suhr[:,:,j]))),
    colorrange=(-1.5e-5, 1.5e-5)
    );

    Colorbar(fig[1:2,3], hm1)

end