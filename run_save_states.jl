include("eddy_paper.jl")

Shr = ShallowWaters.model_setup(output=false,
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
        initpath="/scratch/swilliamson/eddy-stresses/spinup_files/1024_spinup_noslip_10years_050624"
)

hrstates = save_states(Shr)

jldsave("1024_postspinup_thirtydays_hourly_071425.jld2", hrstates=hrstates)