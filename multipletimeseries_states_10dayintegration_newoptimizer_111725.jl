# New structure with variables related to checkpointing,
# will also make it so that the parameters in S.Parameters
# are all constant, nothing changes in time
mutable struct multistatenlp_Chkp{T, S} <: AbstractNLPModel{T,S}
    meta::NLPModelMeta{T,S}
    counters::Counters
    S::ShallowWaters.ModelSetup{T,T}        # model structure
    initial_cond::Array{Array{T,2}, 1}
    # data::Array{Array{Array{T,2}, 1}, 1}    # computed data
    data::Array{Array{T, 3}, 1}
    data_steps::StepRange{Int, Int}         # location of data points temporally
    J::Float64                              # objective function value
    j::Int                                  # for keeping track of location in data
    i::Int                                  # timestep iterator
    t::Int64                                # model time
    avg_eta::Array{T,2}                     # time-averaged eta from integration
    data_avg_eta::Array{T,2}                # time-averaged eta from data
end

# for running with checkpointing
function cpintegrate(chkp, scheme)::Float64

    # calculate layer thicknesses for initial conditions
    ShallowWaters.thickness!(chkp.S.Diag.VolumeFluxes.h, chkp.S.Prog.η, chkp.S.forcing.H)
    ShallowWaters.Ix!(chkp.S.Diag.VolumeFluxes.h_u, chkp.S.Diag.VolumeFluxes.h)
    ShallowWaters.Iy!(chkp.S.Diag.VolumeFluxes.h_v, chkp.S.Diag.VolumeFluxes.h)
    ShallowWaters.Ixy!(chkp.S.Diag.Vorticity.h_q, chkp.S.Diag.VolumeFluxes.h)

    # calculate PV terms for initial conditions
    urhs = chkp.S.Diag.PrognosticVarsRHS.u .= chkp.S.Prog.u
    vrhs = chkp.S.Diag.PrognosticVarsRHS.v .= chkp.S.Prog.v
    ηrhs = chkp.S.Diag.PrognosticVarsRHS.η .= chkp.S.Prog.η

    ShallowWaters.advection_coriolis!(urhs, vrhs, ηrhs, chkp.S.Diag, chkp.S)
    ShallowWaters.PVadvection!(chkp.S.Diag, chkp.S)

    # propagate initial conditions
    copyto!(chkp.S.Diag.RungeKutta.u0, chkp.S.Prog.u)
    copyto!(chkp.S.Diag.RungeKutta.v0, chkp.S.Prog.v)
    copyto!(chkp.S.Diag.RungeKutta.η0, chkp.S.Prog.η)

    # store initial conditions of sst for relaxation
    copyto!(chkp.S.Diag.SemiLagrange.sst_ref, chkp.S.Prog.sst)

    # run integration loop with checkpointing
    chkp.j = 1
    @ad_checkpoint scheme for chkp.i = 1:chkp.S.grid.nt

        t = chkp.t
        i = chkp.i

        # ghost point copy for boundary conditions
        ShallowWaters.ghost_points!(chkp.S.Prog.u, chkp.S.Prog.v, chkp.S.Prog.η, chkp.S)
        copyto!(chkp.S.Diag.RungeKutta.u1, chkp.S.Prog.u)
        copyto!(chkp.S.Diag.RungeKutta.v1, chkp.S.Prog.v)
        copyto!(chkp.S.Diag.RungeKutta.η1, chkp.S.Prog.η)

        if chkp.S.parameters.compensated
            fill!(chkp.S.Diag.Tendencies.du_sum, zero(chkp.S.parameters.Tprog))
            fill!(chkp.S.Diag.Tendencies.dv_sum, zero(chkp.S.parameters.Tprog))
            fill!(chkp.S.Diag.Tendencies.dη_sum, zero(chkp.S.parameters.Tprog))
        end

        for rki = 1:chkp.S.parameters.RKo
            if rki > 1
                ShallowWaters.ghost_points!(
                    chkp.S.Diag.RungeKutta.u1,
                    chkp.S.Diag.RungeKutta.v1,
                    chkp.S.Diag.RungeKutta.η1,
                    chkp.S
                )
            end

            # type conversion for mixed precision
            u1rhs = chkp.S.Diag.PrognosticVarsRHS.u .= chkp.S.Diag.RungeKutta.u1
            v1rhs = chkp.S.Diag.PrognosticVarsRHS.v .= chkp.S.Diag.RungeKutta.v1
            η1rhs = chkp.S.Diag.PrognosticVarsRHS.η .= chkp.S.Diag.RungeKutta.η1

            ShallowWaters.rhs!(u1rhs, v1rhs, η1rhs, chkp.S.Diag, chkp.S, t)          # momentum only
            ShallowWaters.continuity!(u1rhs, v1rhs, η1rhs, chkp.S.Diag, chkp.S, t)   # continuity equation

            if rki < chkp.S.parameters.RKo
                ShallowWaters.caxb!(
                    chkp.S.Diag.RungeKutta.u1,
                    chkp.S.Prog.u,
                    chkp.S.constants.RKbΔt[rki],
                    chkp.S.Diag.Tendencies.du
                )
                ShallowWaters.caxb!(
                    chkp.S.Diag.RungeKutta.v1,
                    chkp.S.Prog.v,
                    chkp.S.constants.RKbΔt[rki],
                    chkp.S.Diag.Tendencies.dv
                )
                ShallowWaters.caxb!(
                    chkp.S.Diag.RungeKutta.η1,
                    chkp.S.Prog.η,
                    chkp.S.constants.RKbΔt[rki],
                    chkp.S.Diag.Tendencies.dη
                )
            end

            if chkp.S.parameters.compensated
                ShallowWaters.axb!(chkp.S.Diag.Tendencies.du_sum, chkp.S.constants.RKaΔt[rki], chkp.S.Diag.Tendencies.du)
                ShallowWaters.axb!(chkp.S.Diag.Tendencies.dv_sum, chkp.S.constants.RKaΔt[rki], chkp.S.Diag.Tendencies.dv)
                ShallowWaters.axb!(chkp.S.Diag.Tendencies.dη_sum, chkp.S.constants.RKaΔt[rki], chkp.S.Diag.Tendencies.dη)
            else
                ShallowWaters.axb!(
                    chkp.S.Diag.RungeKutta.u0,
                    chkp.S.constants.RKaΔt[rki],
                    chkp.S.Diag.Tendencies.du
                )
                ShallowWaters.axb!(
                    chkp.S.Diag.RungeKutta.v0,
                    chkp.S.constants.RKaΔt[rki],
                    chkp.S.Diag.Tendencies.dv
                )
                ShallowWaters.axb!(
                    chkp.S.Diag.RungeKutta.η0,
                    chkp.S.constants.RKaΔt[rki],
                    chkp.S.Diag.Tendencies.dη
                )
            end
        end

        if chkp.S.parameters.compensated
            ShallowWaters.axb!(chkp.S.Diag.Tendencies.du_sum, -1, chkp.S.Diag.Tendencies.du_comp)
            ShallowWaters.axb!(chkp.S.Diag.Tendencies.dv_sum, -1, chkp.S.Diag.Tendencies.dv_comp)
            ShallowWaters.axb!(chkp.S.Diag.Tendencies.dη_sum, -1, chkp.S.Diag.Tendencies.dη_comp)

            ShallowWaters.axb!(chkp.S.Diag.RungeKutta.u0, 1, chkp.S.Diag.Tendencies.du_sum)
            ShallowWaters.axb!(chkp.S.Diag.RungeKutta.v0, 1, chkp.S.Diag.Tendencies.dv_sum)
            ShallowWaters.axb!(chkp.S.Diag.RungeKutta.η0, 1, chkp.S.Diag.Tendencies.dη_sum)

            ShallowWaters.dambmc!(
                chkp.S.Diag.Tendencies.du_comp,
                chkp.S.Diag.RungeKutta.u0,
                chkp.S.Prog.u,
                chkp.S.Diag.Tendencies.du_sum
            )
            ShallowWaters.dambmc!(
                chkp.S.Diag.Tendencies.dv_comp,
                chkp.S.Diag.RungeKutta.v0,
                chkp.S.Prog.v,
                chkp.S.Diag.Tendencies.dv_sum
            )
            ShallowWaters.dambmc!(
                chkp.S.Diag.Tendencies.dη_comp,
                chkp.S.Diag.RungeKutta.η0,
                chkp.S.Prog.η,
                chkp.S.Diag.Tendencies.dη_sum
            )
        end

        ShallowWaters.ghost_points!(
            chkp.S.Diag.RungeKutta.u0,
            chkp.S.Diag.RungeKutta.v0,
            chkp.S.Diag.RungeKutta.η0,
            chkp.S
        )

        u0rhs = chkp.S.Diag.PrognosticVarsRHS.u .= chkp.S.Diag.RungeKutta.u0
        v0rhs = chkp.S.Diag.PrognosticVarsRHS.v .= chkp.S.Diag.RungeKutta.v0
        η0rhs = chkp.S.Diag.PrognosticVarsRHS.η .= chkp.S.Diag.RungeKutta.η0

        if chkp.S.parameters.dynamics == "nonlinear" && chkp.S.grid.nstep_advcor > 0 && (i % chkp.S.grid.nstep_advcor) == 0
            ShallowWaters.UVfluxes!(u0rhs, v0rhs, η0rhs, chkp.S.Diag, chkp.S)
            ShallowWaters.advection_coriolis!(u0rhs, v0rhs, η0rhs, chkp.S.Diag, chkp.S)
        end

        if (chkp.i % chkp.S.grid.nstep_diff) == 0
            ShallowWaters.bottom_drag!(u0rhs, v0rhs, η0rhs, chkp.S.Diag, chkp.S)
            ShallowWaters.diffusion!(u0rhs, v0rhs, chkp.S.Diag, chkp.S)
            ShallowWaters.add_drag_diff_tendencies!(
                chkp.S.Diag.RungeKutta.u0,
                chkp.S.Diag.RungeKutta.v0,
                chkp.S.Diag,
                chkp.S
            )
            ShallowWaters.ghost_points_uv!(
                chkp.S.Diag.RungeKutta.u0,
                chkp.S.Diag.RungeKutta.v0,
                chkp.S
            )
        end

        t += chkp.S.grid.dtint

        u0rhs = chkp.S.Diag.PrognosticVarsRHS.u .= chkp.S.Diag.RungeKutta.u0
        v0rhs = chkp.S.Diag.PrognosticVarsRHS.v .= chkp.S.Diag.RungeKutta.v0
        ShallowWaters.tracer!(i, u0rhs, v0rhs, chkp.S.Prog, chkp.S.Diag, chkp.S)

        if chkp.i in chkp.data_steps

            temp = ShallowWaters.PrognosticVars{Float64}(ShallowWaters.remove_halo(
                chkp.S.Prog.u,
                chkp.S.Prog.v,
                chkp.S.Prog.η,
                chkp.S.Prog.sst,
                chkp.S
            )...)

            # time-average eta
            chkp.avg_eta += temp.η
            chkp.data_avg_eta += chkp.data[3][:,:,chkp.j]

            chkp.J += sum((temp.u .- chkp.data[1][:,:,chkp.j]).^2) + sum((temp.v .- chkp.data[2][:,:,chkp.j]).^2)

            chkp.j += 1

        end

        copyto!(chkp.S.Prog.u, chkp.S.Diag.RungeKutta.u0)
        copyto!(chkp.S.Prog.v, chkp.S.Diag.RungeKutta.v0)
        copyto!(chkp.S.Prog.η, chkp.S.Diag.RungeKutta.η0)

    end

    # add the time-averaged ssh to the loss function
    chkp.J += sum((chkp.avg_eta .- chkp.data_avg_eta).^2) / (chkp.j * 128^2)

    return chkp.J

end

# for running without checkpointing
function integrate(chkp)::Float64

    # calculate layer thicknesses for initial conditions
    ShallowWaters.thickness!(chkp.S.Diag.VolumeFluxes.h, chkp.S.Prog.η, chkp.S.forcing.H)
    ShallowWaters.Ix!(chkp.S.Diag.VolumeFluxes.h_u, chkp.S.Diag.VolumeFluxes.h)
    ShallowWaters.Iy!(chkp.S.Diag.VolumeFluxes.h_v, chkp.S.Diag.VolumeFluxes.h)
    ShallowWaters.Ixy!(chkp.S.Diag.Vorticity.h_q, chkp.S.Diag.VolumeFluxes.h)

    # calculate PV terms for initial conditions
    urhs = chkp.S.Diag.PrognosticVarsRHS.u .= chkp.S.Prog.u
    vrhs = chkp.S.Diag.PrognosticVarsRHS.v .= chkp.S.Prog.v
    ηrhs = chkp.S.Diag.PrognosticVarsRHS.η .= chkp.S.Prog.η

    ShallowWaters.advection_coriolis!(urhs, vrhs, ηrhs, chkp.S.Diag, chkp.S)
    ShallowWaters.PVadvection!(chkp.S.Diag, chkp.S)

    # propagate initial conditions
    copyto!(chkp.S.Diag.RungeKutta.u0, chkp.S.Prog.u)
    copyto!(chkp.S.Diag.RungeKutta.v0, chkp.S.Prog.v)
    copyto!(chkp.S.Diag.RungeKutta.η0, chkp.S.Prog.η)

    # store initial conditions of sst for relaxation
    copyto!(chkp.S.Diag.SemiLagrange.sst_ref, chkp.S.Prog.sst)

    # run integration loop with checkpointing
    chkp.j = 1
    for chkp.i = 1:chkp.S.grid.nt

        t = chkp.t
        i = chkp.i

        # ghost point copy for boundary conditions
        ShallowWaters.ghost_points!(chkp.S.Prog.u, chkp.S.Prog.v, chkp.S.Prog.η, chkp.S)
        copyto!(chkp.S.Diag.RungeKutta.u1, chkp.S.Prog.u)
        copyto!(chkp.S.Diag.RungeKutta.v1, chkp.S.Prog.v)
        copyto!(chkp.S.Diag.RungeKutta.η1, chkp.S.Prog.η)

        if chkp.S.parameters.compensated
            fill!(chkp.S.Diag.Tendencies.du_sum, zero(chkp.S.parameters.Tprog))
            fill!(chkp.S.Diag.Tendencies.dv_sum, zero(chkp.S.parameters.Tprog))
            fill!(chkp.S.Diag.Tendencies.dη_sum, zero(chkp.S.parameters.Tprog))
        end

        for rki = 1:chkp.S.parameters.RKo
            if rki > 1
                ShallowWaters.ghost_points!(
                    chkp.S.Diag.RungeKutta.u1,
                    chkp.S.Diag.RungeKutta.v1,
                    chkp.S.Diag.RungeKutta.η1,
                    chkp.S
                )
            end

            # type conversion for mixed precision
            u1rhs = chkp.S.Diag.PrognosticVarsRHS.u .= chkp.S.Diag.RungeKutta.u1
            v1rhs = chkp.S.Diag.PrognosticVarsRHS.v .= chkp.S.Diag.RungeKutta.v1
            η1rhs = chkp.S.Diag.PrognosticVarsRHS.η .= chkp.S.Diag.RungeKutta.η1

            ShallowWaters.rhs!(u1rhs, v1rhs, η1rhs, chkp.S.Diag, chkp.S, t)          # momentum only
            ShallowWaters.continuity!(u1rhs, v1rhs, η1rhs, chkp.S.Diag, chkp.S, t)   # continuity equation

            if rki < chkp.S.parameters.RKo
                ShallowWaters.caxb!(
                    chkp.S.Diag.RungeKutta.u1,
                    chkp.S.Prog.u,
                    chkp.S.constants.RKbΔt[rki],
                    chkp.S.Diag.Tendencies.du
                )
                ShallowWaters.caxb!(
                    chkp.S.Diag.RungeKutta.v1,
                    chkp.S.Prog.v,
                    chkp.S.constants.RKbΔt[rki],
                    chkp.S.Diag.Tendencies.dv
                )
                ShallowWaters.caxb!(
                    chkp.S.Diag.RungeKutta.η1,
                    chkp.S.Prog.η,
                    chkp.S.constants.RKbΔt[rki],
                    chkp.S.Diag.Tendencies.dη
                )
            end

            if chkp.S.parameters.compensated
                ShallowWaters.axb!(chkp.S.Diag.Tendencies.du_sum, chkp.S.constants.RKaΔt[rki], chkp.S.Diag.Tendencies.du)
                ShallowWaters.axb!(chkp.S.Diag.Tendencies.dv_sum, chkp.S.constants.RKaΔt[rki], chkp.S.Diag.Tendencies.dv)
                ShallowWaters.axb!(chkp.S.Diag.Tendencies.dη_sum, chkp.S.constants.RKaΔt[rki], chkp.S.Diag.Tendencies.dη)
            else
                ShallowWaters.axb!(
                    chkp.S.Diag.RungeKutta.u0,
                    chkp.S.constants.RKaΔt[rki],
                    chkp.S.Diag.Tendencies.du
                )
                ShallowWaters.axb!(
                    chkp.S.Diag.RungeKutta.v0,
                    chkp.S.constants.RKaΔt[rki],
                    chkp.S.Diag.Tendencies.dv
                )
                ShallowWaters.axb!(
                    chkp.S.Diag.RungeKutta.η0,
                    chkp.S.constants.RKaΔt[rki],
                    chkp.S.Diag.Tendencies.dη
                )
            end
        end

        if chkp.S.parameters.compensated
            ShallowWaters.axb!(chkp.S.Diag.Tendencies.du_sum, -1, chkp.S.Diag.Tendencies.du_comp)
            ShallowWaters.axb!(chkp.S.Diag.Tendencies.dv_sum, -1, chkp.S.Diag.Tendencies.dv_comp)
            ShallowWaters.axb!(chkp.S.Diag.Tendencies.dη_sum, -1, chkp.S.Diag.Tendencies.dη_comp)

            ShallowWaters.axb!(chkp.S.Diag.RungeKutta.u0, 1, chkp.S.Diag.Tendencies.du_sum)
            ShallowWaters.axb!(chkp.S.Diag.RungeKutta.v0, 1, chkp.S.Diag.Tendencies.dv_sum)
            ShallowWaters.axb!(chkp.S.Diag.RungeKutta.η0, 1, chkp.S.Diag.Tendencies.dη_sum)

            ShallowWaters.dambmc!(
                chkp.S.Diag.Tendencies.du_comp,
                chkp.S.Diag.RungeKutta.u0,
                chkp.S.Prog.u,
                chkp.S.Diag.Tendencies.du_sum
            )
            ShallowWaters.dambmc!(
                chkp.S.Diag.Tendencies.dv_comp,
                chkp.S.Diag.RungeKutta.v0,
                chkp.S.Prog.v,
                chkp.S.Diag.Tendencies.dv_sum
            )
            ShallowWaters.dambmc!(
                chkp.S.Diag.Tendencies.dη_comp,
                chkp.S.Diag.RungeKutta.η0,
                chkp.S.Prog.η,
                chkp.S.Diag.Tendencies.dη_sum
            )
        end

        ShallowWaters.ghost_points!(
            chkp.S.Diag.RungeKutta.u0,
            chkp.S.Diag.RungeKutta.v0,
            chkp.S.Diag.RungeKutta.η0,
            chkp.S
        )

        u0rhs = chkp.S.Diag.PrognosticVarsRHS.u .= chkp.S.Diag.RungeKutta.u0
        v0rhs = chkp.S.Diag.PrognosticVarsRHS.v .= chkp.S.Diag.RungeKutta.v0
        η0rhs = chkp.S.Diag.PrognosticVarsRHS.η .= chkp.S.Diag.RungeKutta.η0

        if chkp.S.parameters.dynamics == "nonlinear" && chkp.S.grid.nstep_advcor > 0 && (i % chkp.S.grid.nstep_advcor) == 0
            ShallowWaters.UVfluxes!(u0rhs, v0rhs, η0rhs, chkp.S.Diag, chkp.S)
            ShallowWaters.advection_coriolis!(u0rhs, v0rhs, η0rhs, chkp.S.Diag, chkp.S)
        end

        if (chkp.i % chkp.S.grid.nstep_diff) == 0
            ShallowWaters.bottom_drag!(u0rhs, v0rhs, η0rhs, chkp.S.Diag, chkp.S)
            ShallowWaters.diffusion!(u0rhs, v0rhs, chkp.S.Diag, chkp.S)
            ShallowWaters.add_drag_diff_tendencies!(
                chkp.S.Diag.RungeKutta.u0,
                chkp.S.Diag.RungeKutta.v0,
                chkp.S.Diag,
                chkp.S
            )
            ShallowWaters.ghost_points_uv!(
                chkp.S.Diag.RungeKutta.u0,
                chkp.S.Diag.RungeKutta.v0,
                chkp.S
            )
        end

        t += chkp.S.grid.dtint

        u0rhs = chkp.S.Diag.PrognosticVarsRHS.u .= chkp.S.Diag.RungeKutta.u0
        v0rhs = chkp.S.Diag.PrognosticVarsRHS.v .= chkp.S.Diag.RungeKutta.v0
        ShallowWaters.tracer!(i, u0rhs, v0rhs, chkp.S.Prog, chkp.S.Diag, chkp.S)

        if chkp.i in chkp.data_steps

            temp = ShallowWaters.PrognosticVars{Float64}(ShallowWaters.remove_halo(
                chkp.S.Prog.u,
                chkp.S.Prog.v,
                chkp.S.Prog.η,
                chkp.S.Prog.sst,
                chkp.S
            )...)

            # time-average eta
            chkp.avg_eta += temp.η
            chkp.data_avg_eta += chkp.data[3][:,:,chkp.j]

            chkp.J += sum((temp.u .- chkp.data[1][:,:,chkp.j]).^2) + sum((temp.v .- chkp.data[2][:,:,chkp.j]).^2)

            chkp.j += 1

        end

        copyto!(chkp.S.Prog.u, chkp.S.Diag.RungeKutta.u0)
        copyto!(chkp.S.Prog.v, chkp.S.Diag.RungeKutta.v0)
        copyto!(chkp.S.Prog.η, chkp.S.Diag.RungeKutta.η0)

    end

    # add the time-averaged ssh to the loss function
    chkp.J += sum((chkp.avg_eta .- chkp.data_avg_eta).^2) / (chkp.j * 128^2)

    return chkp.J

end

function NLPModels.obj(model, param_guess)

    # Type precision
    T = model.S.parameters.T

    # ensure that the model is reset
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
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=model.S.parameters.Ndays
    )

    model.S = ShallowWaters.model_setup(P)
    model.J = 0.
    model.j = 1
    model.i = 1
    model.t = 0

    data = model.data
    data_steps = model.data_steps
    initial_cond = model.initial_cond

    model.S.Prog.u .= initial_cond[1]
    model.S.Prog.v .= initial_cond[2]
    model.S.Prog.η .= initial_cond[3]

    current = 1
    for m in (model.S.Diag.CNNVars.model_Su, model.S.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                sz = prod(size(array))
                array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
                current += sz
            end
        end
    end

    model.J = integrate(model)
    println("Norm when running integrate: ", model.J)

    return model.J

end

function NLPModels.grad!(model, param_guess, G)

    # Type precision
    T = model.S.parameters.T

    # ensure that the model is reset
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
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=model.S.parameters.Ndays
    )
    model.S = ShallowWaters.model_setup(P)
    model.J = 0.
    model.j = 1
    model.i = 1
    model.t = 0

    S = model.S
    data = model.data
    data_steps = model.data_steps
    initial_cond = model.initial_cond

    S.Prog.u .= initial_cond[1]
    S.Prog.v .= initial_cond[2]
    S.Prog.η .= initial_cond[3]

    current = 1
    for m in (S.Diag.CNNVars.model_Su, S.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    dmodel = Enzyme.make_zero(model)

    snaps = Int(floor(sqrt(model.S.grid.nt)))
    revolve = Revolve(
        snaps;
        verbose=0,
        gc=true,
        write_checkpoints=false
    )

    J = autodiff(
        set_runtime_activity(Enzyme.ReverseWithPrimal),
        cpintegrate,
        Active,
        Duplicated(model, dmodel),
        Const(revolve)
    )[2]

    # Get gradient
    current = 1
    for m in (dmodel.S.Diag.CNNVars.model_Su, dmodel.S.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                sz = prod(size(array))
                G[current:(current + sz - 1)] .= vec(array)
                current += sz
            end
        end
    end

    return nothing

end

function multistatenlp_Chkp{T}(Ndays,param_guess,lower_bound,upper_bound) where {T<:AbstractFloat}

    Plr = ShallowWaters.Parameter(T=T,
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
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=Ndays
    )

    Slr = ShallowWaters.model_setup(Plr)

    Phr = ShallowWaters.Parameter(T=T,
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
        N=1,
        α=2,
        nx=1024,
        Ndays=Ndays
    )
    Shr = ShallowWaters.model_setup(Phr)

    # daily information
    # uhrcg = load_object("./spinup_files/coarsegrainedu_30days_dailysaves_071525.jld2")
    # vhrcg = load_object("./spinup_files/coarsegrainedv_30days_dailysaves_071525.jld2")
    # etahrcg = load_object("./spinup_files/coarsegrainedeta_30days_dailysaves_071525.jld2")
    # data_steps = 225:224:Slr.grid.nt
    # data = [uhrcg[2:11], vhrcg[2:11], etahrcg[2:11]]

    # hourly information
    # every 8 hours is when the timesteps matchup, so I'm doing that frequency for online data
    coarse_grained_hrstates = load_object("./offline_files/hrstates_filtered_downsized_hourly_tendays_uveta_beginsatonehour_102825.jld2")
    uhrcg = coarse_grained_hrstates[1]
    vhrcg = coarse_grained_hrstates[2]
    etahrcg = coarse_grained_hrstates[3]
    data_steps = 75:75:Slr.grid.nt
    data = [uhrcg[:,:,8:8:end], vhrcg[:,:,8:8:end], etahrcg[:,:,8:8:end]]

    u0 = load_object("./spinup_files/coarsegrained_1024_10yearstate_061925.jld2")[1]
    v0 = load_object("./spinup_files/coarsegrained_1024_10yearstate_061925.jld2")[2]
    eta0 = load_object("./spinup_files/coarsegrained_1024_10yearstate_061925.jld2")[3]

    initial_cond = [u0, v0, eta0]

    lvar = lower_bound .* ones(Float64, Lux.parameterlength(Slr.Diag.CNNVars.model_Su) + Lux.parameterlength(Slr.Diag.CNNVars.model_Sv))
    uvar = upper_bound .* ones(Float64, Lux.parameterlength(Slr.Diag.CNNVars.model_Su) + Lux.parameterlength(Slr.Diag.CNNVars.model_Sv))
    meta = NLPModelMeta(Lux.parameterlength(Slr.Diag.CNNVars.model_Su) + Lux.parameterlength(Slr.Diag.CNNVars.model_Sv);
        ncon=0,
        nnzh=0,
        x0=param_guess,
        lvar=lvar,
        uvar=uvar
    )
    counters = Counters()

    return multistatenlp_Chkp{T, typeof(param_guess)}(meta, Counters(), Slr, initial_cond, data, data_steps, 0.0, 1, 1, 0.0, zeros(128,128), zeros(128,128))

end

function run_multistate()

    # result = nothing
    # for ndays = [1, 2, 4, 6, 8, 10]

    #     if ndays === 1
    #         # the initial guess for weights will be the result from the offline problem
    #         param_guess = load_object("./offline_workingresults_5-25-25NN_nobias_1000iterations_1e-3obj_1e-2grad_102025.jld2").solution
    #     else
    #         param_guess = result.solution
    #     end
    #     nlp = multistatenlp_Chkp{Float64}(ndays,param_guess)
    #     qn_options = MadNLP.QuasiNewtonOptions(;max_history=100)
    #     result = madnlp(
    #         nlp;
    #         # linear_solver=LapackCPUSolver,
    #         hessian_approximation=MadNLP.CompactLBFGS,
    #         quasi_newton_options=qn_options,
    #         max_iter=100
    #     )

    # end

    T = Float64

    Plr = ShallowWaters.Parameter(T=T,
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
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=Ndays
    )

    Slr = ShallowWaters.model_setup(Plr)

    param_guess = zeros(Lux.parameterlength(Slr.Diag.CNNVars.model_Su) + Lux.parameterlength(Slr.Diag.CNNVars.model_Sv))
    current = 1
    for model in (Slr.Diag.CNNVars.model_Su, Slr.Diag.CNNVars.model_Sv)
        for layers in model[1]
            for array in layers
                    sz = prod(size(array))
                    param_guess[current:(current + sz - 1)] .= vec(array)
                    current += sz
            end
        end
    end

    # param_guess = load_object("./offline_files/offlineresult_3-25-25_1e-3objective_relu_activation.jld2").solution
    param_guess = load_object("./offline_files/results/offline_5snapshots_111025/result_offline_5snapshots_150iterations_geluactivation_111725.jld2").solution

    # lvar is by default -Inf * ones(Float64, nvar)
    # uvar is by default Inf * ones(Float64, nvar)
    ndays = 1
    lower_bound = -10000
    upper_bound = 10000
    nlp = multistatenlp_Chkp{Float64}(ndays,param_guess,lower_bound,upper_bound)

    # qn_options = MadNLP.QuasiNewtonOptions(;max_history=200)
    result = madnlp(
        nlp;
        # linear_solver=LapackCPUSolver,
        hessian_approximation=MadNLP.CompactLBFGS,
        # quasi_newton_options=qn_options,
        max_iter=100
    )

    # ipopt(nlp, hessian_approximation="limited-memory", limited_memory_max_history=50, max_iter=3)

    jldsave("online_fivedays_100iterations_8hourdata_result_200maxhistory.jld2", result=result)

    return nothing

end

function finite_difference_withnlp(Ndays, xcoord, ycoord)

    # Type precision
    T = Float64

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
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=Ndays
    )

    S0 = ShallowWaters.model_setup(P)

    snaps = Int(floor(sqrt(S0.grid.nt)))
    revolve = Revolve(
        snaps;
        verbose=0,
        gc=true,
        write_checkpoints=false,
        write_checkpoints_filename = "",
        write_checkpoints_period = 2274
    )

    coarse_grained_hrstates = load_object("./hrstates_filtered_downsized_hourly_tendays_uveta_beginsatonehour_102825")
    uhrcg = coarse_grained_hrstates[1]
    vhrcg = coarse_grained_hrstates[2]
    etahrcg = coarse_grained_hrstates[3]
    data_steps = 75:75:S0.grid.nt
    data = [uhrcg[:,:,8:8:end], vhrcg[:,:,8:8:end], etahrcg[:,:,8:8:end]]

    u0 = load_object("./spinup_files/coarsegrained_1024_10yearstate_061925.jld2")[1]
    v0 = load_object("./spinup_files/coarsegrained_1024_10yearstate_061925.jld2")[2]
    eta0 = load_object("./spinup_files/coarsegrained_1024_10yearstate_061925.jld2")[3]

    initial_cond = [u0, v0, eta0]

    param_guess = load_object("./offlineresult_3-25-25_1e-3objective.jld2").solution

    meta = NLPModelMeta(Lux.parameterlength(S0.Diag.CNNVars.model_Su) + Lux.parameterlength(S0.Diag.CNNVars.model_Sv);
        ncon=0,
        nnzh=0,
        x0=param_guess
    )

    S1 = deepcopy(S0)
    chkp1 = multistatenlp_Chkp{T, typeof(param_guess)}(meta,
        Counters(),
        S1,
        initial_cond,
        data,
        data_steps,
        0.0,
        1,
        1,
        0.0,
        zeros(128,128),
        zeros(128,128)
    )
    dchkp1 = Enzyme.make_zero(chkp1)

    # Enzyme deriv
    J = autodiff(
        set_runtime_activity(Enzyme.ReverseWithPrimal),
        multistate_checkpointed_integration,
        Active,
        Duplicated(chkp1, dchkp1),
        Const(revolve)
    )[2]

    enzyme_deriv = dchkp1.S.Diag.CNNVars.model_Su[1][1][1][3, 2, 2, 25] 
    println("Loss when using Enzyme + Checkpointing: $J")
    println("Enzyme derivative: $enzyme_deriv")

    S2 = deepcopy(S0)
    chkp2 = multistatenlp_Chkp{T, typeof(param_guess)}(meta,
        Counters(),
        S2,
        initial_cond,
        data,
        data_steps,
        0.0,
        1,
        1,
        0.0,
        zeros(128,128),
        zeros(128,128)
    )

    @time unperturbed_loss = multistate_checkpointed_integration(chkp2, revolve)
    println("Loss when not using Enzyme: $unperturbed_loss")

    steps = [100, 50, 30, 20, 10, 1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
    diffs = []
    for s in steps

        S3 = deepcopy(S0)
        chkp3 = multistatenlp_Chkp{T, typeof(param_guess)}(meta,
            Counters(),
            S3,
            initial_cond,
            data,
            data_steps,
            0.0,
            1,
            1,
            0.0,
            zeros(128,128),
            zeros(128,128)
        )

        chkp3.S.Diag.CNNVars.model_Su[1][1][1][3, 2, 2, 25] += s

        J = multistate_checkpointed_integration(chkp3, revolve)
        push!(diffs, (J - unperturbed_loss) / s)

    end

    println("Finite difference result: $diffs")

end

# how to save with jld2

# jldsave("exp3_minimizer_initcond_forcing_adjoint_042925.jld2",
#     u = reshape(result.minimizer[1:17292], 131, 132),
#     v = reshape(result.minimizer[17293:34584], 132, 131),
#     eta = reshape(result.minimizer[34585:end-1], 130, 130),
#     Fx0 = result.minimizer[end]
# )
