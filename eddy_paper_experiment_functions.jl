"""
Three files for technical paper:
    1) eddy_paper_integration.jl - contains the time stepping loop, checkpointed
        only include once
    2) eddy_paper_plotting.jl - just a bunch of random plots
    3) eddy_paper.jl - running experiments
"""

function outer(S, dS)
    autodiff(Enzyme.ReverseWithPrimal, integration, Duplicated(S, dS))
    nothing
end

function run_adjoint(::Type{T}=Float32;
    kwargs...
    ) where {T<:AbstractFloat}

    P = ShallowWaters.Parameter(T=T;kwargs...)
    @show typeof(P)
    S = ShallowWaters.model_setup(P)

    dS = Enzyme.Compiler.make_zero(S)
    snaps = Int(floor(sqrt(S.grid.nt)))
    revolve = Revolve{ShallowWaters.ModelSetup}(S.grid.nt,
        snaps;
        verbose=1,
        gc=true,
        write_checkpoints=false,
        write_checkpoints_filename = "technicalpaper_checkingderivatives_30dayrun_withcheckpointing_120924",
        write_checkpoints_period = 224
    )

    # autodiff(Enzyme.ReverseWithPrimal, checkpointed_integration, Duplicated(S, dS), Const(revolve))
    # autodiff(Enzyme.ReverseWithPrimal, integration, Duplicated(S, dS))
    # outer()

    S = Reactant.to_rarray(S)
    dS = Reactant.to_rarray(dS)
    revolve = Reactant.to_rarray(revolve)
    compiled_outer = @compile outer(S, dS)

    compiled_outer = outer
    compiled_outer(S, dS)
    return S, dS

end

function run_adjoint_plusfd(::Type{T}=Float32;
    kwargs...
    ) where {T<:AbstractFloat}

    P = ShallowWaters.Parameter(T=T;kwargs...)
    S = ShallowWaters.model_setup(P)

    S.Diag.NNVars.weights_corner=0.01 .* randn(2, 22)
    S.Diag.NNVars.weights_center=0.01 .* randn(1, 17)

    dS = Enzyme.Compiler.make_zero(Core.Typeof(S), IdDict(), S)
    snaps = Int(floor(sqrt(S.grid.nt)))
    revolve = Revolve{ShallowWaters.ModelSetup}(S.grid.nt,
        snaps;
        verbose=1,
        gc=true,
        write_checkpoints=false,
        write_checkpoints_filename = "technicalpaper_checkingderivatives_30dayrun_withcheckpointing_120924",
        write_checkpoints_period = 224
    )

    # autodiff(Enzyme.ReverseWithPrimal, checkpointed_integration, Duplicated(S, dS), Const(revolve))
    autodiff(Enzyme.ReverseWithPrimal, integration, Duplicated(S, dS))

    enzyme_deriv = dS.Prog.u[72,120]

    @show enzyme_deriv

    steps = [50, 40, 30, 20, 10, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
    # steps = [3, 2, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]

    S_outer = ShallowWaters.model_setup(P)

    snaps = Int(floor(sqrt(S_outer.grid.nt)))
    revolve = Revolve{ShallowWaters.ModelSetup}(S_outer.grid.nt, snaps;
        verbose=1,
        gc=true,
        write_checkpoints=false
    )

    # J_outer = checkpointed_integration(S_outer, revolve)
    J_outer = integration(S_outer)

    diffs = []

    for s in steps

        S_inner = ShallowWaters.model_setup(P)

        S_inner.Prog.u[72,120] += s

        # J_inner = checkpointed_integration(S_inner, revolve)
        J_inner = integration(S_inner)

        push!(diffs, (J_inner - J_outer) / s)

    end

    return S, dS, diffs, enzyme_deriv

end

function finite_difference_only(dS,x_coord,y_coord;kwargs...)

    enzyme_deriv = dS.Prog.u[x_coord,y_coord]

    steps = [50, 40, 30, 20, 10, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
    # steps = [3, 2, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]

    P = ShallowWaters.Parameter(T=dS.parameters.T;kwargs...)
    S_outer = ShallowWaters.model_setup(P)

    snaps = Int(floor(sqrt(S_outer.grid.nt)))
    revolve = Revolve{ShallowWaters.ModelSetup}(S_outer.grid.nt, snaps;
        verbose=1,
        gc=true,
        write_checkpoints=false
    )

    J_outer = integration(S_outer)

    diffs = []

    for s in steps

        S_inner = ShallowWaters.model_setup(output=false,
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
        # νB=1000,
        nx=128,
        Ndays=30,
        # initial_cond="ncfile",
        # initpath="./run_0001/"
        )
        S_inner.Prog.u[x_coord, y_coord] += s

        # J_inner = checkpointed_integration(S_inner, revolve)
        J_inner = integration(S_inner)

        push!(diffs, (J_inner - J_outer) / s)

    end

    return diffs, enzyme_deriv

end

"""

"""
function enzyme_derivatives(::Type{T}=Float32;     # number format
    kwargs...                               # all additional parameters
    ) where {T<:AbstractFloat}

    P = ShallowWaters.Parameter(T=T;kwargs...)
    S = ShallowWaters.model_setup(P)

    # run the forward model and save all the states

    # setup
    Diag = S.Diag
    Prog = S.Prog

    @unpack u,v,η,sst = Prog
    @unpack u0,v0,η0 = Diag.RungeKutta
    @unpack u1,v1,η1 = Diag.RungeKutta
    @unpack du,dv,dη = Diag.Tendencies
    @unpack du_sum,dv_sum,dη_sum = Diag.Tendencies
    @unpack du_comp,dv_comp,dη_comp = Diag.Tendencies

    @unpack um,vm = Diag.SemiLagrange

    @unpack dynamics,RKo,RKs,tracer_advection = S.parameters
    @unpack time_scheme,compensated = S.parameters
    @unpack RKaΔt,RKbΔt = S.constants
    @unpack Δt_Δ,Δt_Δs = S.constants

    @unpack nt,dtint = S.grid
    @unpack nstep_advcor,nstep_diff,nadvstep,nadvstep_half = S.grid

    S.parameters.data = zeros(S.grid.nt)

    # calculate layer thicknesses for initial conditions
    ShallowWaters.thickness!(Diag.VolumeFluxes.h,η,S.forcing.H)
    ShallowWaters.Ix!(Diag.VolumeFluxes.h_u,Diag.VolumeFluxes.h)
    ShallowWaters.Iy!(Diag.VolumeFluxes.h_v,Diag.VolumeFluxes.h)
    ShallowWaters.Ixy!(Diag.Vorticity.h_q,Diag.VolumeFluxes.h)

    # calculate PV terms for initial conditions
    urhs = convert(Diag.PrognosticVarsRHS.u,u)
    vrhs = convert(Diag.PrognosticVarsRHS.v,v)
    ηrhs = convert(Diag.PrognosticVarsRHS.η,η)

    ShallowWaters.advection_coriolis!(urhs,vrhs,ηrhs,Diag,S)
    ShallowWaters.PVadvection!(Diag,S)

    # propagate initial conditions
    copyto!(u0,u)
    copyto!(v0,v)
    copyto!(η0,η)

    # store initial conditions of sst for relaxation
    copyto!(Diag.SemiLagrange.sst_ref,sst)

    states_for_enzyme = []
    for j = 1:S.grid.nt
        push!(states_for_enzyme, deepcopy(S))
        push!(states, deepcopy(S))
        _ = onestep(S, j)
    end

    dS = Enzyme.make_zero(S)
    derivatives = []

    # run the backwards problem and save all the adjoint computed derivatives
    for k = Base.reverse(1:S.grid.nt)
        St = pop!(states_for_enzyme)

        autodiff(Enzyme.Reverse, onestep, Duplicated(St, dS), Const(k))

        if k in 1:224:S.grid.nt
            push!(derivatives, deepcopy(dS))
        end

    end

    return S, dS, derivatives, states

end

"""
Optim experiment functions
"""
function optim_run(::Type{T}=Float32;
    kwargs...
    ) where {T<:AbstractFloat}

end