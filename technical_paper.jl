"""
Three files for technical paper:
    1) technical_paper_integration.jl - contains the time stepping loop, checkpointed
        only include once
    2) technical_paper_plotting.jl - just a bunch of random plots
    3) technical_paper.jl - running experiments
"""


include("technical_paper_integration.jl")

function run_adjoint_plusfd(::Type{T}=Float32;     # number format
    kwargs...                               # all additional parameters
    ) where {T<:AbstractFloat}

    P = ShallowWaters.Parameter(T=T;kwargs...)
    S = ShallowWaters.model_setup(P)


    dS = Enzyme.Compiler.make_zero(Core.Typeof(S), IdDict(), S)
    snaps = Int(floor(sqrt(S.grid.nt)))
    revolve = Revolve{ShallowWaters.ModelSetup}(S.grid.nt,
        snaps;
        verbose=1,
        gc=true,
        write_checkpoints=false,
        write_checkpoints_filename = "technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_startingfromspinup_everytimestep_500m_4months_float32_112124",
        write_checkpoints_period = 224
    )

    autodiff(Enzyme.ReverseWithPrimal, checkpointed_integration, Duplicated(S, dS), Const(revolve))

    enzyme_deriv = dS.Prog.u[62,20]

    @show enzyme_deriv

    # steps = [50, 40, 30, 20, 10, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
    steps = [1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]

    S_outer = ShallowWaters.model_setup(P)

    snaps = Int(floor(sqrt(S_outer.grid.nt)))
    revolve = Revolve{ShallowWaters.ModelSetup}(S_outer.grid.nt, snaps;
        verbose=1,
        gc=true,
        write_checkpoints=false
    )

    J_outer = checkpointed_integration(S_outer, revolve)

    diffs = []

    for s in steps

        S_inner = ShallowWaters.model_setup(P)

        S_inner.Prog.u[62,20] += s

        J_inner = checkpointed_integration(S_inner, revolve)

        push!(diffs, (J_inner - J_outer) / s)

    end

    return diffs, enzyme_deriv, S, dS

end

diffs, enzyme_deriv, S, dS = run_adjoint_plusfd(
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
    α=2,
    # νB 
    nx=128,
    Ndays=30
    # initial_cond="ncfile",
    # initpath="./data_files_gamma0.3/10yearspinup_128_noslipbc_noforcing_float64prog"
)

# @save "technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_startingfromspinup_finalprimal_struct_500mdepth_4months_float64start_112224.jld2" S
# @save "technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_startingfromspinup_finaladjoint_struct_500mdepth_4months_float64start_112224.jld2" dS
# @save "technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_startingfromspinup_fdcheck_vector_500mdepth_4months_float64start_112224.jld2" diffs


# create_adjoint_gif()

# _, energy = ShallowWaters.run_model(T = Float64,
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
#     # νB 
#     nx=128,
#     Ndays=10*12*30
#     # initial_cond="ncfile",
#     # initpath="eddy-stresses/data_files_gamma0.3/128_10yearspinup_noforcing_noslipbc_float64params"
# )