# New structure with variables related to checkpointing,
# will also make it so that the parameters in S.Parameters
# are all constant, nothing changes in time
mutable struct exp1_Chkp{T1,T2}
    S::ShallowWaters.ModelSetup{T1,T2}      # model structure
    # data::Matrix{Float32}                 # computed data
    data_steps::StepRange{Int, Int}         # location of data points temporally
    J::Float64                              # objective function value
    j::Int                                  # for keeping track of location in data
    i::Int                                  # timestep iterator
    t::Int64                                # model time
end

function outer(S, dS, revolve)
    autodiff(Enzyme.ReverseWithPrimal, checkpointed_integration, Duplicated(S, dS), Const(revolve))
    nothing
end

function exp1_compute_objective()

end

function exp1_compute_gradient(G, data, data_steps, Ndays)

    # Type precision
    T = Float32

    P = ShallowWaters.Parameter(T=T;
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
        nn_forcing_dissipation=false,
        handwritten=false,
        N=1,
        α=2,
        nx=128,
        Ndays=Ndays,
        initial_cond="rest",
        # initpath="./data_files_forkf/128_spinup_noforcing/"
    )

    S = ShallowWaters.model_setup(P)

    data_steps = (S.grid.nt - (7*224)):1:S.grid.nt

    snaps = Int(floor(sqrt(S.grid.nt)))
    revolve = Revolve{exp1_Chkp{T, T}}(S.grid.nt,
        snaps;
        verbose=1,
        gc=true,
        write_checkpoints=false,
        write_checkpoints_filename = "",
        write_checkpoints_period = 224
    )

    chkp = exp1_Chkp{T, T}(S,
        data_steps,
        0.0,
        1,
        1,
        0.0
    )
    dchkp = Enzyme.make_zero(chkp)

    chkp_prim = deepcopy(chkp)
    J = exp1_checkpointed_integration(chkp_prim, revolve)
    println("Cost without AD: $J")

    @time J = autodiff(
        set_runtime_activity(Enzyme.ReverseWithPrimal),
        exp1_checkpointed_integration,
        Active,
        Duplicated(chkp, dchkp),
        Const(revolve)
    )[2]
    println("Cost with AD: $J")

    # Get gradient
    @unpack u, v, η = dchkp.S.Prog
    G .= [vec(u); vec(v); vec(η)]

    return nothing

end

# """
# runs the optim experiment
# """
# obj_fg = Optim.only_fg!(FG)

# param_guess = 0.001.*randn(61)

# result = Optim.optimize(obj_fg,
# param_guess,
# Optim.LBFGS(),
# Optim.Options(
# iterations = 1)
# )

# """
# runs both the adjoint problem and a finite difference check in one go
# """
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
#     Ndays=1,
#     # initial_cond="ncfile",
#     # initpath="./run_0001/"
# )

# # 109, 129 is the entry with the largest relative error after 30 days

# # 72, 120 is the entry with the largest relative error after 10 days

# """
# runs the finite difference check, requires prior computation of derivatives with enzyme
# """
# diffs, deriv = finite_difference_only(dS30, 20, 20,
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
#     nn_forcing_dissipation=true,
#     handwritten=false,
#     α=2,
#     # νB=1000,
#     nx=128,
#     Ndays=30,
#     # initial_cond="ncfile",
#     # initpath="./run_0001/"
# )

# # @save "technicalpaper_75km_timeavgobj_onlyfinalmonth_everytimestep_startingfromspinup_finalprimal_struct_12months_120524.jld2" S
# # @save "technicalpaper_75km_timeavgobj_onlyfinalmonth_everytimestep_startingfromspinup_finaladjoint_struct_12months_120524.jld2" dS
# # @save "technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_fourtimesviscosity_startingfromspinup_fdcheck_vector_500mdepth_12months_float32start_120324.jld2" diffs

# """
# For saving all derivatives computed with Enzyme (runs the backwards pass step by step)
# """
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
# )

# """
# Just runs the shallow water model
# """
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
