"""
runs just the adjoint computation
"""
S30, dS30 = run_adjoint(
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
    nn_forcing_dissipation=true,
    handwritten=false,
    α=2,
    # νB=1000,
    nx=128,
    Ndays=1,
    # initial_cond="ncfile",
    # initpath="./run_0001/"
)

"""
runs the optim experiment
"""
obj_fg = Optim.only_fg!(FG)

param_guess = zeros(61)

result = Optim.optimize(obj_fg,
param_guess,
Optim.LBFGS(),
Optim.Options(
iterations = 1)
)

"""
runs both the adjoint problem and a finite difference check in one go
"""
# S30, dS30, diffs30, enzyme_deriv30 = run_adjoint_plusfd(
#     output=false,
#     L_ratio=1,
#     g=9.81,
#     H=500,
#     wind_forcing_x="double_gyre",
#     Lx=3840e3,
#     seasonal_wind_x=false,
#     topography="flat",
#     bc="nonperiodic",
#     bottom_drag="quadratic",
#     α=2,
#     # νB=1000,
#     nx=128,
#     Ndays=30,
#     # initial_cond="ncfile",
#     # initpath="./run_0001/"
# )

# 109, 129 is the entry with the largest relative error after 30 days

# 72, 120 is the entry with the largest relative error after 10 days

"""
runs the finite difference check, requires prior computation of derivatives with enzyme
"""
# diffs, deriv = finite_difference_only(dS30, 20, 20, output=false,
#     L_ratio=1,
#     g=9.81,
#     H=500,
#     wind_forcing_x="double_gyre",
#     Lx=3840e3,
#     seasonal_wind_x=false,
#     topography="flat",
#     bc="nonperiodic",
#     bottom_drag="quadratic",
#     α=2,
#     # νB=1000,
#     nx=128,
#     Ndays=30,
#     # initial_cond="ncfile",
#     # initpath="./run_0001/"
# )

# @save "technicalpaper_75km_timeavgobj_onlyfinalmonth_everytimestep_startingfromspinup_finalprimal_struct_12months_120524.jld2" S
# @save "technicalpaper_75km_timeavgobj_onlyfinalmonth_everytimestep_startingfromspinup_finaladjoint_struct_12months_120524.jld2" dS
# @save "technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_fourtimesviscosity_startingfromspinup_fdcheck_vector_500mdepth_12months_float32start_120324.jld2" diffs

"""
For saving all derivatives computed with Enzyme (runs the backwards pass step by step)
"""
# S, dS, derivatives, states = enzyme_derivatives(
#     output=false,
#     L_ratio=1,
#     g=9.81,
#     H=500,
#     wind_forcing_x="double_gyre",
#     Lx=3840e3,
#     seasonal_wind_x=false,
#     topography="flat",
#     bc="nonperiodic",
#     bottom_drag="quadratic",
#     α=2,
#     # νB=1000,
#     nx=128,
#     Ndays=30,
#     # initial_cond="ncfile",
#     # initpath="./run_0001/"
# );

"""
Just runs the shallow water model
"""
# _, energy = ShallowWaters.run_model(
#     output=true,
#     L_ratio=1,
#     g=9.81,
#     H=500,
#     wind_forcing_x="double_gyre",
#     Lx=3840e3,
#     seasonal_wind_x=false,
#     topography="flat",
#     bc="nonperiodic",
#     bottom_drag="quadratic",
#     α=2,
#     nx=50,
#     Ndays=12*30*10
#     # initial_cond="ncfile",
#     # initpath="./data_files_gamma0.3/10yearspinup_128_noslipbc_fromrest_float32params"
# )
