include("eddy_paper.jl")

# Shr = ShallowWaters.model_setup(output=false,
#         L_ratio=1,
#         g=9.81,
#         H=500,
#         wind_forcing_x="double_gyre",
#         Lx=3840e3,
#         seasonal_wind_x=false,
#         topography="flat",
#         bc="nonperiodic",
#         bottom_drag="quadratic",
#         tracer_advection=false,
#         tracer_relaxation=false,
#         α=2,
#         nx=1024,
#         Ndays=30,
#         initial_cond="ncfile",
#         initpath="/scratch/swilliamson/eddy-stresses/spinup_files/1024_spinup_noslip_10years_050624"
# )

# hrstates = save_states(Shr)

# jldsave("1024_postspinup_thirtydays_hourly_071425.jld2", hrstates=hrstates)

u_hr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/u.nc", "u")
v_hr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/v.nc", "v")
eta_hr = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/eta.nc", "eta")

Slr = ShallowWaters.model_setup(output=false,
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
        Ndays=1
)
coarse_grain(u_hr[:,:,1], v_hr[:,:,1], η_hr[:,:,1], 1024, S_lr)
