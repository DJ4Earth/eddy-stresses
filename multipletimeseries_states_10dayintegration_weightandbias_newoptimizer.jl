# New structure with variables related to checkpointing,
# will also make it so that the parameters in S.Parameters
# are all constant, nothing changes in time
mutable struct multistatenlp_Chkp{T, S} <: AbstractNLPModel{T,S}
    meta::NLPModelMeta{T,S}
    counters::Counters
    S::ShallowWaters.ModelSetup{T,T}        # model structure
    initial_cond::Array{Array{T,2}, 1}
    data::Array{Array{Array{T,2}, 1}, 1}    # computed data
    data_steps::StepRange{Int, Int}         # location of data points temporally
    J::Float64                              # objective function value
    j::Int                                  # for keeping track of location in data
    i::Int                                  # timestep iterator
    t::Int64                                # model time
end

# for running with checkpointing
function multistate_checkpointed_integration(chkp, scheme)

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

         temp = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(
            chkp.S.Prog.u,
            chkp.S.Prog.v,
            chkp.S.Prog.η,
            chkp.S.Prog.sst,
            chkp.S
        )...)

        chkp.J += sum((temp.u .- chkp.data[1][chkp.j]).^2 + (temp.v .- chkp.data[2][chkp.j]).^2)

        # storing the objective function over time
        # S.parameters.data[S.parameters.i] = S.parameters.J / length((S.grid.nt - 30*224):1:S.parameters.i)

        chkp.j += 1

    end

    ##### time-averaging the objective function #######
    # chkp.J = chkp.J / length((chkp.S.grid.nt - 7*224):1:chkp.S.grid.nt) # time-averaging
    ##########################################################

    copyto!(chkp.S.Prog.u, chkp.S.Diag.RungeKutta.u0)
    copyto!(chkp.S.Prog.v, chkp.S.Diag.RungeKutta.v0)
    copyto!(chkp.S.Prog.η, chkp.S.Diag.RungeKutta.η0)

    end

    return chkp.J

end

# for running without checkpointing
function multistate_integration(chkp)

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

         temp = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(
            chkp.S.Prog.u,
            chkp.S.Prog.v,
            chkp.S.Prog.η,
            chkp.S.Prog.sst,
            chkp.S
        )...)

        chkp.J += sum((temp.u - chkp.data[1][chkp.j]).^2) + sum((temp.v - chkp.data[2][chkp.j]).^2)

        # storing the objective function over time
        # S.parameters.data[S.parameters.i] = S.parameters.J / length((S.grid.nt - 30*224):1:S.parameters.i)

        chkp.j += 1

    end

    ##### time-averaging the objective function #######
    # chkp.J = chkp.J / length((chkp.S.grid.nt - 7*224):1:chkp.S.grid.nt) # time-averaging
    ##########################################################

    copyto!(chkp.S.Prog.u, chkp.S.Diag.RungeKutta.u0)
    copyto!(chkp.S.Prog.v, chkp.S.Diag.RungeKutta.v0)
    copyto!(chkp.S.Prog.η, chkp.S.Diag.RungeKutta.η0)

    end

    return chkp.J

end

function NLPModels.obj(model, param_guess)

    # Type precision
    T = Float64

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

    model.J = multistate_integration(model)

    return model.J

end

function NLPModels.grad!(model, param_guess, G)

    # Type precision
    T = Float64

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

    J = @time autodiff(
        set_runtime_activity(Enzyme.ReverseWithPrimal),
        multistate_integration,
        Active,
        Duplicated(model, dmodel)
    )[2]

    # Get gradient
    G = zeros(Lux.parameterlength(S.Diag.CNNVars.model_Su) + Lux.parameterlength(S.Diag.CNNVars.model_Sv))
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

function multistatenlp_Chkp{T}(Ndays,param_guess) where {T<:AbstractFloat}

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
        Ndays=Ndays
    )

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
        N=1,
        α=2,
        nx=1024,
        Ndays=Ndays
    )

    # daily information
    # hrstates = load_object("./spinup_files/1024_coarsegrained_tendays_dailysaves_062425.jld2")
    # data = hrstates[2:end]

    uhrcg = load_object("./spinup_files/coarsegrainedu_30days_dailysaves_071525.jld2")
    vhrcg = load_object("./spinup_files/coarsegrainedv_30days_dailysaves_071525.jld2")
    etahrcg = load_object("./spinup_files/coarsegrainedeta_30days_dailysaves_071525.jld2")

    data_steps = 225:225:Slr.grid.nt

    data = [uhrcg[2:11], vhrcg[2:11], etahrcg[2:11]]
    u0 = load_object("./spinup_files/coarsegrained_1024_10yearstate_061925.jld2")[1]
    v0 = load_object("./spinup_files/coarsegrained_1024_10yearstate_061925.jld2")[2]
    eta0 = load_object("./spinup_files/coarsegrained_1024_10yearstate_061925.jld2")[3]

    initial_cond = [u0, v0, eta0]

    meta = NLPModelMeta(Lux.parameterlength(Slr.Diag.CNNVars.model_Su) + Lux.parameterlength(Slr.Diag.CNNVars.model_Sv);
        ncon=0,
        nnzh=0,
        x0=param_guess
    )
    counters = Counters()

    return multistatenlp_Chkp{T, typeof(param_guess)}(meta, Counters(), Slr, initial_cond, data, data_steps, 0.0, 1, 1, 0.0)

end

function run_multistate()

    result = nothing
    for ndays = [1, 2, 4, 6, 8, 10]

        if ndays === 1
            param_guess = load_object("./result_offline_newoptimizer_onesnapshot_zhangetalmodel_082825.jld2").solution
        else
            param_guess = result.solution
        end
        nlp = multistatenlp_Chkp{Float64}(ndays,param_guess)
        qn_options = MadNLP.QuasiNewtonOptions(;max_history=20)
        result = madnlp(nlp;
            hessian_approximation=MadNLP.CompactLBFGS,
            quasi_newton_options=qn_options,
            tol=1e-3
        )

    end

    return result

end

# how to save with jld2

# jldsave("exp3_minimizer_initcond_forcing_adjoint_042925.jld2",
#     u = reshape(result.minimizer[1:17292], 131, 132),
#     v = reshape(result.minimizer[17293:34584], 132, 131),
#     eta = reshape(result.minimizer[34585:end-1], 130, 130),
#     Fx0 = result.minimizer[end]
# )
