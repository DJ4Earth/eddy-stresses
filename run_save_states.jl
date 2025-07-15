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

u_hr = ncread("./run_0001/u.nc", "u")
v_hr = ncread("./run_0001/v.nc", "v")
eta_hr = ncread("./run_0001/eta.nc", "eta")

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

u = []
v = []
eta = []
for j = 1:30

        push!(u, ShallowWaters.coarse_grain_u(u_hr[:,:,j], 1024, Slr))
        push!(v, ShallowWaters.coarse_grain_v(v_hr[:,:,j], 1024, Slr))
        push!(eta, ShallowWaters.coarse_grain_eta(eta_hr[:,:,j], 1024, Slr))

end

jldsave("hru_coarsegrained_30days_dailysaves.jld2", u=u)
jldsave("hrv_coarsegrained_30days_dailysaves.jld2", v=v)
jldsave("hreta_coarsegrained_30days_dailysaves.jld2", eta=eta)

