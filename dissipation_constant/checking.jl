T = Float64
Ndays = 90

P1 = ShallowWaters.Parameter(T=T,
    output=true,
    # output_dt=168,
    output_dt=1,
    L_ratio=1,
    g=9.81,
    H=500,
    cfl=.898,
    wind_forcing_x="double_gyre",
    Lx=3840e3,
    seasonal_wind_x=false,
    adv_scheme="ArakawaHsu",
    topography="flat",
    bc="nonperiodic",
    bottom_drag="quadratic",
    diffusion="constant",        # this is the only new parameter to be adjusted in the new spinups
    tracer_advection=false,
    tracer_relaxation=false,
    zb_forcing_momentum=false,
    zb_forcing_dissipation=false,
    zb_filtered=true,
    nn_forcing_momentum=false,
    nn_forcing_dissipation=false,
    N=1,
    α=2,
    nx=128,
    Ndays=Ndays,
    initial_cond="ncfile",
    initpath="./dissipation_constant/spinup_files/10yearspinup_128_noslipbc_noforcing_float64prog/"
);

S1 = ShallowWaters.model_setup(P1);

ShallowWaters.time_integration(S1)

P2 = ShallowWaters.Parameter(T=T,
    output=true,
    # output_dt=168,
    output_dt=1,
    L_ratio=1,
    g=9.81,
    H=500,
    cfl=.898,
    wind_forcing_x="double_gyre",
    Lx=3840e3,
    seasonal_wind_x=false,
    adv_scheme="Sadourny",
    topography="flat",
    bc="nonperiodic",
    bottom_drag="quadratic",
    diffusion="constant",        # this is the only new parameter to be adjusted in the new spinups
    tracer_advection=false,
    tracer_relaxation=false,
    zb_forcing_momentum=false,
    zb_forcing_dissipation=false,
    zb_filtered=true,
    nn_forcing_momentum=false,
    nn_forcing_dissipation=false,
    N=1,
    α=2,
    nx=128,
    Ndays=Ndays,
    initial_cond="ncfile",
    initpath="./dissipation_constant/spinup_files/10yearspinup_128_noslipbc_noforcing_float64prog/"
);

S2 = ShallowWaters.model_setup(P2);

ShallowWaters.time_integration(S2)

ah = []
s = []

uah = ncread("./run_0002/u.nc", "u");
vah = ncread("./run_0002/v.nc", "v");

us = ncread("./run_0003/u.nc", "u");
vs = ncread("./run_0003/v.nc", "v");

for j = 1:2251
    push!(ah, sum(uah[:,1:end-1,j].^2 .+ vah[1:end-1,:,j].^2))
    push!(s, sum(us[:,1:end-1,j].^2 .+ vs[1:end-1,:,j].^2))
end

fig = Figure();
lines(fig[1,1], ah[1:24:end] ./ 128^2)
lines!(fig[1,1], s[1:24:end] ./ 128^2)