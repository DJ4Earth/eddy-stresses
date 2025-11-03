SNN = ShallowWaters.model_setup(T=Float32;
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
    Ndays=1,
    initial_cond="ncfile",
    initpath="./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825",
    init_starti=1
)

param_guess = load_object("./offline_files/offlineresult_3-25-25_1e-3objective_relu_activation.jld2").solution
# param_guess = zeros(Lux.parameterlength(SNN.Diag.CNNVars.model_Su) + Lux.parameterlength(SNN.Diag.CNNVars.model_Sv))
current = 1
for model in (SNN.Diag.CNNVars.model_Su, SNN.Diag.CNNVars.model_Sv)
    for layers in model[1]
        for array in layers
                sz = prod(size(array))
                param_guess[current:(current + sz - 1)] .= vec(array)
                current += sz
        end
    end
end

ShallowWaters.CNN_momentum(SNN.Prog.u, SNN.Prog.v, SNN)

T11 = SNN.Diag.CNNVars.T11
T12 = SNN.Diag.CNNVars.T12
T22 = SNN.Diag.CNNVars.T22

SZB = ShallowWaters.model_setup(T=Float32;
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
    Ndays=Ndays,
    initial_cond="ncfile",
    initpath="./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825",
    init_starti=1
)
ShallowWaters.ZB_momentum(SZB.Prog.u, SZB.Prog.v, SZB, SZB.Diag)

# high-resolution T's
# ubar = ShallowWaters.coarse_grain_u(hru[:,:,1], 1024, SZB)
# vbar = ShallowWaters.coarse_grain_v(hrv[:,:,1], 1024, SZB)

# ubarsq = ShallowWaters.coarse_grain_u(hru[:,:,1].^2, 1024, SZB)
# vbarsq = ShallowWaters.coarse_grain_v(hrv[:,:,1].^2, 1024, SZB)

nx = SZB.grid.nx
ny = SZB.grid.ny
nux = SZB.grid.nux
nuy = SZB.grid.nuy
nvx = SZB.grid.nvx
nvy = SZB.grid.nvy
halo = SZB.grid.halo
haloη = SZB.grid.haloη

ucg = load_object("./offline_files/coarsegrained_hr_ubar_ubarsq_t1_foroffline_101525.jld2")
vcg = load_object("./offline_files/coarsegrained_hr_vbar_vbarsq_t1_foroffline_101525.jld2")
uvbar = load_object("./offline_files/coarsegrained_hr_uvbar_t1_foroffline_101525.jld2")

ubar = ucg[1]
ubarsq = ucg[2]

vbar = vcg[1]
vbarsq = vcg[2]

ubarh = cat(zeros(T,nux+2*halo,halo),cat(zeros(T,halo,nuy),ubar,zeros(T,halo,nuy),dims=1),zeros(T,nux+2*halo,halo),dims=2)
ubarhsq = cat(zeros(T,nux+2*halo,halo),cat(zeros(T,halo,nuy),ubarsq,zeros(T,halo,nuy),dims=1),zeros(T,nux+2*halo,halo),dims=2)

vbarh = cat(zeros(T,nvx+2*halo,halo),cat(zeros(T,halo,nvy),vbar,zeros(T,halo,nvy),dims=1),zeros(T,nvx+2*halo,halo),dims=2)
vbarhsq = cat(zeros(T,nvx+2*halo,halo),cat(zeros(T,halo,nvy),vbarsq,zeros(T,halo,nvy),dims=1),zeros(T,nvx+2*halo,halo),dims=2)

uvbarh = cat(zeros(T,nx+2*haloη,haloη),cat(zeros(T,haloη,ny),uvbar,zeros(T,haloη,ny),dims=1),zeros(T,nx+2*haloη,haloη),dims=2)

T12_true = ShallowWaters.Iy(ubarh)[2:end-1,2:end-1] .* ShallowWaters.Ix(vbarh)[2:end-1,2:end-1] - ShallowWaters.Ixy(uvbarh)

T11_true = (ShallowWaters.Ixy(ShallowWaters.Iy(ubarh)[2:end-1,2:end-1])).^2 - ShallowWaters.Ixy(ShallowWaters.Iy(ubarhsq)[2:end-1,2:end-1])
T22_true = ShallowWaters.Ixy((ShallowWaters.Ix(vbarh)[2:end-1,2:end-1])).^2 - ShallowWaters.Ixy(ShallowWaters.Ix(vbarhsq)[2:end-1,2:end-1])

ζD_filtered = SZB.Diag.ZBVars.ζD_filtered
ζDhat_filtered = SZB.Diag.ZBVars.ζDhat_filtered
trace_filtered = SZB.Diag.ZBVars.trace_filtered

ζsqT = SZB.Diag.ZBVars.ζsqT
ζDT = SZB.Diag.ZBVars.ζDT
ζDhat = SZB.Diag.ZBVars.ζDhat

denom = SZB.grid.Δ^2 * SZB.grid.scale

## "True"

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

## NN

ax1, hm1 = heatmap(fig[2,1], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T11,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\tilde{T}_{11}"),
colorrange=(-maximum(abs.(T11)),
maximum(abs.(T11)))
);
Colorbar(fig[2,2], hm1)

ax2, hm2 = heatmap(fig[2,3], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T12,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\tilde{T}_{12}"),
colorrange=(-maximum(abs.(T12)),
maximum(abs.(T12)))
);
Colorbar(fig[2,4], hm2)

ax3, hm3 = heatmap(fig[2,5], LinRange(0, 3840, 128),
LinRange(0, 3840, 128),
T22,
colormap=:balance,
axis=(xlabel="km", ylabel="km", title=L"\tilde{T}_{22}"),
colorrange=(-maximum(abs.(T22)),
maximum(abs.(T22)))
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


# the following is to create the coarse grained uv term for computing the off-diagonal entries in T

uhrh = cat(zeros(T,1023+2*halo,halo),cat(zeros(T,halo,1024),hru[:,:,1],zeros(T,halo,1024),dims=1),zeros(T,1023+2*halo,halo),dims=2)
vhrh = cat(zeros(T,1024+2*halo,halo),cat(zeros(T,halo,1023),hrv[:,:,1],zeros(T,halo,1023),dims=1),zeros(T,1024+2*halo,halo),dims=2)

# moving to hr corner grid and cut off the halo

uhrq = ShallowWaters.Iy(uhrh)[2:end-1,2:end-1]
vhrq = ShallowWaters.Ix(vhrh)[2:end-1,2:end-1]

uhrT = zeros(1024,1024)
vhrT = zeros(1024,1024)

ShallowWaters.Ixy!(uhrT,uhrq)
ShallowWaters.Ixy!(vhrT,vhrq)

uvbar = ShallowWaters.coarse_grain_eta(uhrT .* vhrT, 1024, SZB)