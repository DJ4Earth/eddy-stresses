"""
Three files for technical paper:
    1) technical_paper_integration.jl - contains the time stepping loop, checkpointed
        only include once
    2) technical_paper_plotting.jl - just a bunch of random plots
    3) technical_paper.jl - running experiments
"""


include("../ShallowWaters.jl/src/ShallowWaters.jl")
using .ShallowWaters

using Enzyme
using Checkpointing

mutable struct MyPrognosticVars{T<:AbstractFloat}
    u::Array{T,2}           # u-velocity
    nu::Array{T,2}           # sea surface height / interface displacement
end

function checkpointed_integration(S, scheme)
    h = S.Diag.VolumeFluxes.h
    H = S.forcing.H
    nu = S.Prog.η
    @inbounds for i in eachindex(nu)
        h[i] = nu[i]
    end
    
    @checkpoint_struct scheme S for S.parameters.i = 1:S.grid.nt
        Diag = S.Diag
        Prog = S.Prog

        halo = S.grid.halo

        # undo scaling as well
        @views ucut = S.Prog.u[halo+1:end-halo,halo+1:end-halo]
        ηcut = S.Prog.η

        temp = MyPrognosticVars{Float64}(ucut, ηcut) 

        energy_lr = first(temp.u.^2)

        S.parameters.J += energy_lr
    end

    return S.parameters.J

end

function mymodel_setup(P::Parameter)
    T = P.T
    Tprog = P.Tprog

    G = ShallowWaters.Grid{T,Tprog}(P)
    C = ShallowWaters.Constants{T,Tprog}(P,G)
    F = ShallowWaters.Forcing{T}(P,G)

    Prog = ShallowWaters.initial_conditions(Tprog,G,P,C)
    Diag = ShallowWaters.preallocate(T,Tprog,G)

    S = ShallowWaters.ModelSetup{T,Tprog}(P,G,C,F,Prog,Diag,0)

    return S

end
function run_adjoint_plusfd(::Type{T}=Float32;     # number format
    kwargs...                               # all additional parameters
    ) where {T<:AbstractFloat}

    P = ShallowWaters.Parameter(T=T;kwargs...)
    S = mymodel_setup(P)


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
