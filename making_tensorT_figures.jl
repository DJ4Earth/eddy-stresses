function compute_Ts(j, Ndays)

    ulr = ncread("./spinup_files/128_cginitcond_1year_postspinup_noforcing/u.nc", "u")
    vlr = ncread("./spinup_files/128_cginitcond_1year_postspinup_noforcing/v.nc", "v")
    etalr = ncread("./spinup_files/128_cginitcond_1year_postspinup_noforcing/eta.nc", "eta")

    Ppred = ShallowWaters.Parameter(T=Float64;
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
        initpath="./spinup_files/128_cginitcond_1year_postspinup_noforcing/",
        init_starti=1
    );
    Spred = ShallowWaters.model_setup(Ppred);

    PNN = ShallowWaters.Parameter(T=Float64;
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
        initpath="./spinup_files/128_cginitcond_1year_postspinup_noforcing/",
        init_starti=1
    );
    SNN = ShallowWaters.model_setup(PNN);

    cgstates = load_object("./offline_files/1024_filtered_downsized_uveta_10days_postspinup_hourlysaves_111925.jld2")
    ucg = cgstates[1]
    vcg = cgstates[2]
    etacg = cgstates[3]

    u, v, eta = ShallowWaters.add_halo(ucg[:,:,j], vcg[:,:,j], etacg[:,:,j], zeros(128,128), SNN)
    snapshot = [u, v]

    param_guess = load_object("./tuned_weights/result_offline_150iterations_geluactivation_111925.jld2").solution;
    current = 1
    for model in (SNN.Diag.CNNVars.model_Su, SNN.Diag.CNNVars.model_Sv)
        for layers in model[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
                    # println(norm((reshape(param_guess[current:(current + sz - 1)], size(array)...) .- array)))
                    current += sz
            end
        end
    end

    ShallowWaters.CNN_momentum(snapshot[1], snapshot[2], Spred);
    ShallowWaters.CNN_momentum(snapshot[1], snapshot[2], SNN);

    T11_pred = Spred.Diag.CNNVars.T11;
    T12_pred = Spred.Diag.CNNVars.T12;
    T22_pred = Spred.Diag.CNNVars.T22;

    T11_NN = SNN.Diag.CNNVars.T11;
    T12_NN = SNN.Diag.CNNVars.T12;
    T22_NN = SNN.Diag.CNNVars.T22;

    PZB = ShallowWaters.Parameter(T=Float64;
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
        zb_forcing_dissipation=true,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=false,
        N=1,
        α=2,
        nx=128,
        Ndays=1,
        initial_cond="ncfile",
        initpath="./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825",
        init_starti=1
    );
    SZB = ShallowWaters.model_setup(PZB);
    ShallowWaters.ZB_momentum(snapshot[1], snapshot[2], SZB, SZB.Diag);

    # high-resolution T's

    nx = SZB.grid.nx
    ny = SZB.grid.ny
    nux = SZB.grid.nux
    nuy = SZB.grid.nuy
    nvx = SZB.grid.nvx
    nvy = SZB.grid.nvy
    halo = SZB.grid.halo
    haloη = SZB.grid.haloη

    true_Ts = load_object("./offline_files/trueTs_filtered_downsized_T11T22T12_hourlysaves_111925.jld2")

    T11_true = true_Ts[1][:,:,j];
    T22_true = true_Ts[2][:,:,j];
    T12_true = true_Ts[3][:,:,j];

    ζD_filtered = SZB.Diag.ZBVars.ζD_filtered;
    ζDhat_filtered = SZB.Diag.ZBVars.ζDhat_filtered;
    trace_filtered = SZB.Diag.ZBVars.trace_filtered;

    return T11_true, T12_true, T22_true, ζD_filtered, ζDhat_filtered, trace_filtered, T11_NN, T12_NN, T22_NN, T11_pred, T12_pred, T22_pred

end

function plots()

PZB = ShallowWaters.Parameter(T=Float64;
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
    zb_forcing_dissipation=true,
    zb_filtered=true,
    nn_forcing_momentum=false,
    nn_forcing_dissipation=false,
    N=1,
    α=2,
    nx=128,
    Ndays=1,
    initial_cond="ncfile",
    initpath="./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825",
    init_starti=1
);
SZB = ShallowWaters.model_setup(PZB);

# snapshot to use, can be anything within j = 3:3:241
# these are the values the offline weights were trained on
Ndays = 1
T11_true, T12_true, T22_true, ζD_filtered, ζDhat_filtered, trace_filtered, T11_NN, T12_NN, T22_NN, T11_pred, T12_pred, T22_pred = compute_Ts(j, Ndays);
denom = SZB.grid.Δ^2 * SZB.grid.scale

## True versus NN T's

fig = Figure(fontsize = 15,size=(900,450));

ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T11_true,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"T_{11}"),
colorrange=(-maximum(abs.(T11_true)),
maximum(abs.(T11_true)))
);
Colorbar(fig[1,2], hm1)

ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T12_true,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"T_{12}"),
colorrange=(-maximum(abs.(T12_true)),
maximum(abs.(T12_true)))
);
Colorbar(fig[1,4], hm2)

ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T22_true,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"T_{22}"),
colorrange=(-maximum(abs.(T22_true)),
maximum(abs.(T22_true)))
);
Colorbar(fig[1,6], hm3)

ax1, hm1 = heatmap(fig[2,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T11_NN,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\tilde{T}_{11}"),
colorrange=(-maximum(abs.(T11_true)),
maximum(abs.(T11_true)))
);
Colorbar(fig[2,2], hm1)

ax2, hm2 = heatmap(fig[2,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T12_NN,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\tilde{T}_{12}"),
colorrange=(-maximum(abs.(T12_true)),
maximum(abs.(T12_true)))
);
Colorbar(fig[2,4], hm2)

ax3, hm3 = heatmap(fig[2,5], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T22_NN,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\tilde{T}_{22}"),
colorrange=(-maximum(abs.(T22_true)),
maximum(abs.(T22_true)))
);
Colorbar(fig[2,6], hm3)

# Zanna-Bolton output versus the NN output

fig = Figure(fontsize = 15, size=(900,450));
ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
(trace_filtered - ζD_filtered)./denom,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\zeta^2 - \zeta D"),
colorrange=((-maximum(abs.((trace_filtered - ζD_filtered)./denom))),
maximum(abs.((trace_filtered - ζD_filtered)./denom)))
);
Colorbar(fig[1,2], hm1)

ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
ζDhat_filtered ./ denom,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\zeta \hat{D}"),
colorrange=(-maximum(abs.(ζDhat_filtered ./ denom)),
maximum(abs.(ζDhat_filtered ./ denom)))
);
Colorbar(fig[1,4], hm2)

ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
(trace_filtered + ζD_filtered)./denom,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\zeta^2 + \zeta D"),
colorrange=(-maximum(abs.((trace_filtered + ζD_filtered)./denom)),
maximum(abs.((trace_filtered + ζD_filtered)./denom)))
);
Colorbar(fig[1,6], hm3)

ax1, hm1 = heatmap(fig[2,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T11_NN,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\tilde{T}_{11}"),
colorrange=((-maximum(abs.(T11_NN))),
maximum(abs.(T11_NN)))
);
Colorbar(fig[2,2], hm1)

ax2, hm2 = heatmap(fig[2,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T12_NN,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\tilde{T}_{12}"),
colorrange=(-maximum(abs.(T12_NN)),
maximum(abs.(T12_NN)))
);
Colorbar(fig[2,4], hm2)

ax3, hm3 = heatmap(fig[2,5], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T22_NN,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\tilde{T}_{22}"),
colorrange=(-maximum(abs.(T22_NN)),
maximum(abs.(T22_NN)))
);
Colorbar(fig[2,6], hm3)

# Prediction versus the NN output

fig = Figure(fontsize = 15, size=(900,450));
ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T11_pred,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\tilde{T}_{11}(-)"),
colorrange=((-maximum(abs.(T11_pred))),
maximum(abs.(T11_pred)))
);
Colorbar(fig[1,2], hm1)

ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T12_pred,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\tilde{T}_{12}(-)"),
colorrange=(-maximum(abs.(T12_pred)),
maximum(abs.(T12_pred)))
);
Colorbar(fig[1,4], hm2)

ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T22_pred,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\tilde{T}_{22}(-)"),
colorrange=(-maximum(abs.(T22_pred)),
maximum(abs.(T22_pred)))
);
Colorbar(fig[1,6], hm3)

ax1, hm1 = heatmap(fig[2,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T11_NN,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\tilde{T}_{11}"),
colorrange=((-maximum(abs.(T11_NN))),
maximum(abs.(T11_NN)))
);
Colorbar(fig[2,2], hm1)

ax2, hm2 = heatmap(fig[2,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T12_NN,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\tilde{T}_{12}"),
colorrange=(-maximum(abs.(T12_NN)),
maximum(abs.(T12_NN)))
);
Colorbar(fig[2,4], hm2)

ax3, hm3 = heatmap(fig[2,5], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T22_NN,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\tilde{T}_{22}"),
colorrange=(-maximum(abs.(T22_NN)),
maximum(abs.(T22_NN)))
);
Colorbar(fig[2,6], hm3)

# ZB
fig = Figure(fontsize = 15, size=(1500,400));
ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
(trace_filtered - ζD_filtered)./denom,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\zeta^2 - \zeta D"),
colorrange=((-maximum(abs.((trace_filtered - ζD_filtered)./denom))),
maximum(abs.((trace_filtered - ζD_filtered)./denom)))
);
Colorbar(fig[1,2], hm1)

ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
ζDhat_filtered ./ denom,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\zeta \hat{D}"),
colorrange=(-maximum(abs.(ζDhat_filtered ./ denom)),
maximum(abs.(ζDhat_filtered ./ denom)))
);
Colorbar(fig[1,4], hm2)

ax3, hm3 = heatmap(fig[1,5], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
(trace_filtered + ζD_filtered)./denom,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\zeta^2 + \zeta D"),
colorrange=(-maximum(abs.((trace_filtered + ζD_filtered)./denom)),
maximum(abs.((trace_filtered + ζD_filtered)./denom)))
);
Colorbar(fig[1,6], hm3)

# true?

model = InitWeightsModel{Float64}();

fig = Figure(fontsize = 15, size=(1500,400));
ax1, hm1 = heatmap(fig[1,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
model.T11,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"T_{11}"),
colorrange=((-maximum(abs.(model.T11))),
maximum(abs.(model.T11)))
);
Colorbar(fig[1,2], hm1)

ax2, hm2 = heatmap(fig[1,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
(model.T22),
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"T_{22}"),
colorrange=(-maximum(abs.((model.T22))),
maximum(abs.((model.T22))))
);
Colorbar(fig[1,4], hm2)

ax2, hm2 = heatmap(fig[1,5], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
(model.T12),
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"T_{22}"),
colorrange=(-maximum(abs.((model.T12))),
maximum(abs.((model.T12))))
);
Colorbar(fig[1,6], hm2)


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

end