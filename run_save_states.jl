include("eddy_paper.jl")

# Shr = ShallowWaters.model_setup(output=true,
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
#         initpath="./spinup_files/1024_spinup_noslip_10years_050624"
# )

# hrstates = save_states(Shr)

# jldsave("1024_postspinup_thirtydays_hourly_071625.jld2", hrstates=hrstates)

u_hr = ncread("./spinup_files/1024_30days_postspinup_noslip_071625/u.nc", "u")
v_hr = ncread("./spinup_files/1024_30days_postspinup_noslip_071625/v.nc", "v")
eta_hr = ncread("./spinup_files/1024_30days_postspinup_noslip_071625/eta.nc", "eta")

# Slr = ShallowWaters.model_setup(output=false,
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
#         nx=128,
#         Ndays=20
# )

# u = []
# v = []
# eta = []

# for j = 1:20

#         push!(u, ShallowWaters.coarse_grain_u(u_hr[:,:,j], 1024, Slr))
#         push!(v, ShallowWaters.coarse_grain_v(v_hr[:,:,j], 1024, Slr))
#         push!(eta, ShallowWaters.coarse_grain_eta(eta_hr[:,:,j], 1024, Slr))

# end

# jldsave("1024_coarsegrainedu_20days_hourlysaves_071525.jld2", u=u)
# jldsave("1024_coarsegrainedv_20days_hourlysaves_071525.jld2", v=v)
# jldsave("1024_coarsegrainedeta_20days_hourlysaves_071525.jld2", eta=eta)

Shr = ShallowWaters.model_setup(output=true,
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
        α=2,
        nx=1024,
        Ndays=30,
        initial_cond="ncfile",
        initpath="./spinup_files/1024_spinup_noslip_10years_050624"
)
