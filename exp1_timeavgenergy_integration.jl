# don't quite remember what this is from, but presumably from computing 
# individual derivatives using Enzyme and saving them when I don't use 
# checkpointing
function set_value(S)
    S.parameters.νB *= 1.0
end

# for running with checkpointing
function exp1_checkpointed_integration(chkp, scheme)::Float64

    # calculate layer thicknesses for initial conditions
    ShallowWaters.thickness!(chkp.S.Diag.VolumeFluxes.h, chkp.S.Prog.η, chkp.S.forcing.H)
    ShallowWaters.Ix!(chkp.S.Diag.VolumeFluxes.h_u, chkp.S.Diag.VolumeFluxes.h)
    ShallowWaters.Iy!(chkp.S.Diag.VolumeFluxes.h_v, chkp.S.Diag.VolumeFluxes.h)
    ShallowWaters.Ixy!(chkp.S.Diag.Vorticity.h_q, chkp.S.Diag.VolumeFluxes.h)

    # calculate PV terms for initial conditions
    urhs = convert(chkp.S.Diag.PrognosticVarsRHS.u, chkp.S.Prog.u)
    vrhs = convert(chkp.S.Diag.PrognosticVarsRHS.v, chkp.S.Prog.v)
    ηrhs = convert(chkp.S.Diag.PrognosticVarsRHS.η, chkp.S.Prog.η)

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
    @checkpoint_struct scheme chkp for chkp.i = 1:chkp.S.grid.nt

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
            u1rhs = convert(chkp.S.Diag.PrognosticVarsRHS.u, chkp.S.Diag.RungeKutta.u1)
            v1rhs = convert(chkp.S.Diag.PrognosticVarsRHS.v, chkp.S.Diag.RungeKutta.v1)
            η1rhs = convert(chkp.S.Diag.PrognosticVarsRHS.η, chkp.S.Diag.RungeKutta.η1)

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

        u0rhs = convert(chkp.S.Diag.PrognosticVarsRHS.u, chkp.S.Diag.RungeKutta.u0)
        v0rhs = convert(chkp.S.Diag.PrognosticVarsRHS.v, chkp.S.Diag.RungeKutta.v0)
        η0rhs = convert(chkp.S.Diag.PrognosticVarsRHS.η, chkp.S.Diag.RungeKutta.η0)

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

    u0rhs = convert(chkp.S.Diag.PrognosticVarsRHS.u, chkp.S.Diag.RungeKutta.u0)
    v0rhs = convert(chkp.S.Diag.PrognosticVarsRHS.v, chkp.S.Diag.RungeKutta.v0)
    ShallowWaters.tracer!(i, u0rhs, v0rhs, chkp.S.Prog, chkp.S.Diag, chkp.S)

    #### Energy objective function, time averaged
    if chkp.i in chkp.data_steps

        temp = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(
            chkp.S.Prog.u,
            chkp.S.Prog.v,
            chkp.S.Prog.η,
            chkp.S.Prog.sst,
            chkp.S
        )...)

        energy_lr = (sum(temp.u.^2) + sum(temp.v.^2)) / (chkp.S.grid.nx * chkp.S.grid.ny)

        chkp.J = chkp.J + energy_lr

        # storing the objective function over time
        # S.parameters.data[S.parameters.i] = S.parameters.J / length((S.grid.nt - 30*224):1:S.parameters.i)
    
    end

    ##### time-averaging the objective function #######
    chkp.J = chkp.J / length((chkp.S.grid.nt - 7*224):1:chkp.S.grid.nt) # time-averaging
    ##########################################################

    copyto!(chkp.S.Prog.u, chkp.S.Diag.RungeKutta.u0)
    copyto!(chkp.S.Prog.v, chkp.S.Diag.RungeKutta.v0)
    copyto!(chkp.S.Prog.η, chkp.S.Diag.RungeKutta.η0)
    end

    return chkp.J

end

# for running without checkpointing
function exp1_integration(S)

    # calculate layer thicknesses for initial conditions
    ShallowWaters.thickness!(chkp.S.Diag.VolumeFluxes.h, chkp.S.Prog.η, chkp.S.forcing.H)
    ShallowWaters.Ix!(chkp.S.Diag.VolumeFluxes.h_u, chkp.S.Diag.VolumeFluxes.h)
    ShallowWaters.Iy!(chkp.S.Diag.VolumeFluxes.h_v, chkp.S.Diag.VolumeFluxes.h)
    ShallowWaters.Ixy!(chkp.S.Diag.Vorticity.h_q, chkp.S.Diag.VolumeFluxes.h)

    # calculate PV terms for initial conditions
    urhs = convert(chkp.S.Diag.PrognosticVarsRHS.u, chkp.S.Prog.u)
    vrhs = convert(chkp.S.Diag.PrognosticVarsRHS.v, chkp.S.Prog.v)
    ηrhs = convert(chkp.S.Diag.PrognosticVarsRHS.η, chkp.S.Prog.η)

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

        t = chkp.S.t
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
            u1rhs = convert(chkp.S.Diag.PrognosticVarsRHS.u, chkp.S.Diag.RungeKutta.u1)
            v1rhs = convert(chkp.S.Diag.PrognosticVarsRHS.v, chkp.S.Diag.RungeKutta.v1)
            η1rhs = convert(chkp.S.Diag.PrognosticVarsRHS.η, chkp.S.Diag.RungeKutta.η1)

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

        u0rhs = convert(chkp.S.Diag.PrognosticVarsRHS.u, chkp.S.Diag.RungeKutta.u0)
        v0rhs = convert(chkp.S.Diag.PrognosticVarsRHS.v, chkp.S.Diag.RungeKutta.v0)
        η0rhs = convert(chkp.S.Diag.PrognosticVarsRHS.η, chkp.S.Diag.RungeKutta.η0)

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

    u0rhs = convert(chkp.S.Diag.PrognosticVarsRHS.u, chkp.S.Diag.RungeKutta.u0)
    v0rhs = convert(chkp.S.Diag.PrognosticVarsRHS.v, chkp.S.Diag.RungeKutta.v0)
    ShallowWaters.tracer!(i, u0rhs, v0rhs, chkp.S.Prog, chkp.S.Diag, chkp.S)

    #### Energy objective function, time averaged
    temp = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(chkp.S.Prog.u,
    chkp.S.Prog.v,
    chkp.Prog.η,
    chkp.S.Prog.sst,S)...)

    if chkp.i in (chkp.S.grid.nt - 7*224):1:chkp.S.grid.nt

        energy_lr = (sum(temp.u.^2) + sum(temp.v.^2)) / (chkp.S.grid.nx * chkp.S.grid.ny)

        chkp.J = chkp.J + energy_lr

        # storing the objective function over time
        # S.parameters.data[S.parameters.i] = S.parameters.J / length((S.grid.nt - 30*224):1:S.parameters.i)

    end

    ##### time-averaging the objective function #######
    S.parameters.J = S.parameters.J / length((S.grid.nt - 7*224):1:S.grid.nt) # time-averaging
    ##########################################################

    copyto!(chkp.S.Prog.u, chkp.S.Diag.RungeKutta.u0)
    copyto!(chkp.S.Prog.v, chkp.S.Diag.RungeKutta.v0)
    copyto!(chkp.S.Prog.η, chkp.S.Diag.RungeKutta.η0)
    end

    return chkp.J

end

# for running a single step of integration, returns the value of the objective function
function onestep(S, step)

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
    t = S.t
    i = S.parameters.i

    # ghost point copy for boundary conditions
    ShallowWaters.ghost_points!(u,v,η,S)
    copyto!(u1,u)
    copyto!(v1,v)
    copyto!(η1,η)


    if compensated
        fill!(du_sum,zero(Tprog))
        fill!(dv_sum,zero(Tprog))
        fill!(dη_sum,zero(Tprog))
    end

    for rki = 1:RKo
        if rki > 1
            ShallowWaters.ghost_points!(u1,v1,η1,S)
        end

        # type conversion for mixed precision
        u1rhs = convert(Diag.PrognosticVarsRHS.u,u1)
        v1rhs = convert(Diag.PrognosticVarsRHS.v,v1)
        η1rhs = convert(Diag.PrognosticVarsRHS.η,η1)

        ShallowWaters.rhs!(u1rhs,v1rhs,η1rhs,Diag,S,t)          # momentum only
        ShallowWaters.continuity!(u1rhs,v1rhs,η1rhs,Diag,S,t)   # continuity equation

        if rki < RKo
            ShallowWaters.caxb!(u1,u,RKbΔt[rki],du)   #u1 .= u .+ RKb[rki]*Δt*du
            ShallowWaters.caxb!(v1,v,RKbΔt[rki],dv)   #v1 .= v .+ RKb[rki]*Δt*dv
            ShallowWaters.caxb!(η1,η,RKbΔt[rki],dη)   #η1 .= η .+ RKb[rki]*Δt*dη
        end

        if compensated      # accumulate tendencies
            ShallowWaters.axb!(du_sum,RKaΔt[rki],du)
            ShallowWaters.axb!(dv_sum,RKaΔt[rki],dv)
            ShallowWaters.axb!(dη_sum,RKaΔt[rki],dη)
        else    # sum RK-substeps on the go
            ShallowWaters.axb!(u0,RKaΔt[rki],du)          #u0 .+= RKa[rki]*Δt*du
            ShallowWaters.axb!(v0,RKaΔt[rki],dv)          #v0 .+= RKa[rki]*Δt*dv
            ShallowWaters.axb!(η0,RKaΔt[rki],dη)          #η0 .+= RKa[rki]*Δt*dη
        end
    end

    if compensated
        # add compensation term to total tendency
        ShallowWaters.axb!(du_sum,-1,du_comp)
        ShallowWaters.axb!(dv_sum,-1,dv_comp)
        ShallowWaters.axb!(dη_sum,-1,dη_comp)

        ShallowWaters.axb!(u0,1,du_sum)   # update prognostic variable with total tendency
        ShallowWaters.axb!(v0,1,dv_sum)
        ShallowWaters.axb!(η0,1,dη_sum)

        ShallowWaters.dambmc!(du_comp,u0,u,du_sum)    # compute new compensation
        ShallowWaters.dambmc!(dv_comp,v0,v,dv_sum)
        ShallowWaters.dambmc!(dη_comp,η0,η,dη_sum)
    end

    ShallowWaters.ghost_points!(u0,v0,η0,S)

    # type conversion for mixed precision
    u0rhs = convert(Diag.PrognosticVarsRHS.u,u0)
    v0rhs = convert(Diag.PrognosticVarsRHS.v,v0)
    η0rhs = convert(Diag.PrognosticVarsRHS.η,η0)

    # ADVECTION and CORIOLIS TERMS
    # although included in the tendency of every RK substep,
    # only update every nstep_advcor steps if nstep_advcor > 0
    if dynamics == "nonlinear" && nstep_advcor > 0 && (i % nstep_advcor) == 0
        ShallowWaters.UVfluxes!(u0rhs,v0rhs,η0rhs,Diag,S)
        ShallowWaters.advection_coriolis!(u0rhs,v0rhs,η0rhs,Diag,S)
    end

    # DIFFUSIVE TERMS - SEMI-IMPLICIT EULER
    # use u0 = u^(n+1) to evaluate tendencies, add to u0 = u^n + rhs
    # evaluate only every nstep_diff time steps
    if (S.parameters.i % nstep_diff) == 0
        ShallowWaters.bottom_drag!(u0rhs,v0rhs,η0rhs,Diag,S)
        ShallowWaters.diffusion!(u0rhs,v0rhs,Diag,S)
        ShallowWaters.add_drag_diff_tendencies!(u0,v0,Diag,S)
        ShallowWaters.ghost_points_uv!(u0,v0,S)
    end

    t += dtint

    # TRACER ADVECTION
    u0rhs = convert(Diag.PrognosticVarsRHS.u,u0) 
    v0rhs = convert(Diag.PrognosticVarsRHS.v,v0)
    ShallowWaters.tracer!(i,u0rhs,v0rhs,Prog,Diag,S)

    #### Energy objective function, time averaged ###################################
    temp = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(S.Prog.u,
    S.Prog.v,
    S.Prog.η,
    S.Prog.sst,S)...)

    if step in (S.grid.nt - 7*224):1:S.grid.nt # currently accumulating over one week of the one month integration

        energy_lr = (sum(temp.u.^2) + sum(temp.v.^2)) / (S.grid.nx * S.grid.ny)

        S.parameters.J = S.parameters.J + energy_lr

        # storing the objective function over time
        S.parameters.data[step] = S.parameters.J / length((S.grid.nt - 30*224):1:step)
    
    end
    ##################################################################################

    # Copy back from substeps
    copyto!(u,u0)
    copyto!(v,v0)
    copyto!(η,η0)

    ##### use if time-averaging the objective function #######
    S.parameters.J = S.parameters.J / length((S.grid.nt - 7*224):1:S.grid.nt) # time-averaging
    ##########################################################

    # Copy back from substeps
    copyto!(u,u0)
    copyto!(v,v0)
    copyto!(η,η0)

    return S.parameters.J

end