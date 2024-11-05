# Assuming the information we have is solely a spacially averaged kinetic energy

include("../ShallowWaters.jl/src/ShallowWaters.jl")
using .ShallowWaters

using Enzyme#main
using Checkpointing, HDF5, Serialization
using Plots, NetCDF, JLD2

Enzyme.API.looseTypeAnalysis!(true)

using Parameters
using Optim
using LaTeXStrings

function checkpointed_integration(S, scheme)

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

    # run integration loop with checkpointing
    loop(S, scheme)

    return S.parameters.J

end

function loop(S,scheme)

    @checkpoint_struct scheme S for S.parameters.i = 1:S.grid.nt

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

        # Copy back from substeps
        copyto!(u,u0)
        copyto!(v,v0)
        copyto!(η,η0)

    end

    ##### Energy objective function ###########
    temp = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(S.Prog.u,
    S.Prog.v,
    S.Prog.η,
    S.Prog.sst,S)...)

    energy_lr = (sum(temp.u.^2) + sum(temp.v.^2)) #/ (S.grid.nx * S.grid.ny)
    S.parameters.J = energy_lr
    ###########################################

    return nothing

end

function check_derivative(Ndays, H)

    S = ShallowWaters.model_setup(output=false,
        L_ratio=1,
        g=9.81,
        H=H,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        α=2,
        nx=128,
        Ndays = Ndays,
        initial_cond="ncfile",
        initpath="./data_files_gamma0.3/128_spinup_noforcing_noslipbc/"
        )

    dS = Enzyme.Compiler.make_zero(Core.Typeof(S), IdDict(), S)
    snaps = Int(floor(sqrt(S.grid.nt)))
    revolve = Revolve{ShallowWaters.ModelSetup}(S.grid.nt,
        snaps;
        verbose=1,
        gc=true,
        write_checkpoints=false
        # write_checkpoints_filename = "technicalpaper_5000m_period286_1yearintegration_110124.h5",
        # write_checkpoints_period = 286
    )

    autodiff(Enzyme.ReverseWithPrimal, checkpointed_integration, Duplicated(S, dS), Const(revolve))

    @save "technicalpaper_finalprimal_struct_500mdepth_1year_nocDfield_110524.jld2" S
    @save "technicalpaper_finaladjoint_struct_500mdepth_1year_nocDfield_110524.jld2" dS

    enzyme_deriv = dS.constants.cD

    @show enzyme_deriv

    steps = [50, 40, 30, 20, 10, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]

    S_outer = ShallowWaters.model_setup(output=false,
        L_ratio=1,
        g=9.81,
        H=H,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        α=2,
        nx=128,
        Ndays = Ndays,
        initial_cond="ncfile",
        initpath="./data_files_gamma0.3/128_spinup_noforcing_noslipbc/"
        )

    snaps = Int(floor(sqrt(S_outer.grid.nt)))
    revolve = Revolve{ShallowWaters.ModelSetup}(S_outer.grid.nt, snaps; 
        verbose=1, 
        gc=true, 
        write_checkpoints=false
    )

    J_outer = checkpointed_integration(S_outer, revolve)

    diffs = []

    for s in steps

        S_inner = ShallowWaters.model_setup(output=false,
            L_ratio=1,
            g=9.81,
            H=H,
            wind_forcing_x="double_gyre",
            Lx=3840e3,
            seasonal_wind_x=false,
            topography="flat",
            bc="nonperiodic",
            bottom_drag="quadratic",
            α=2,
            nx=128,
            Ndays = Ndays,
            initial_cond="ncfile",
            initpath="./data_files_gamma0.3/128_spinup_noforcing_noslipbc/"
            )

        S_inner.constants.cD += s

        J_inner = checkpointed_integration(S_inner, revolve)

        push!(diffs, (J_inner - J_outer) / s)

    end

    return diffs, enzyme_deriv, S, dS

end

"""
Mostly figure generation, I just wanted to be able to run include("technical_paper.jl")
without all of this also running. The function deserialize opens any saved checkpoints
"""
function deserialize(x)
    s = IOBuffer(x)
    Serialization.deserialize(s)
end

function stuff()

    # loss function is final spatially averaged energy
    # initial condition sensitivity

    adjoint500 = h5open("adjoint_technicalpaper_5000m_period286_1yearintegration_110124.h5", "r")
    blob = read(adjoint500["1"])
    adj_500_1 = deserialize(blob)

    state_derivs = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(dS2.Prog.u,
    dS2.Prog.v,
    dS2.Prog.η,
    dS2.Prog.sst,dS2)...)

    one = heatmap(LinRange(0, 3840, 127),
        LinRange(0, 3840, 128),
        state_derivs.u[:, :]',
        c=:balance,
        xlabel="x (km)",
        xguidefontsize=13,
        ylabel="y (km)",
        yguidefontsize=13,
        title=L"\partial J / \partial u(t_0)",
        plot_titlefontsize=13,
        colorbar_title=L"m",
        colorbar_titlefontsize=13,
        clim=(-3e17,3e17),
        colorbar=:false,
        dpi=300
    )

    two = heatmap(LinRange(0, 3840, 128),
        LinRange(0, 3840, 127),
        state_derivs.v[:, :]',
        c=:balance,
        xlabel="x (km)",
        xguidefontsize=13,
        ylabel="y (km)",
        yguidefontsize=13,
        title=L"\partial J / \partial v(t_0)",
        plot_titlefontsize=13,
        colorbar_title=L"      ",
        colorbar_titlefontsize=13,
        clim=(-3e17,3e17),
        dpi=300
    )

    plot(one, two, layout=grid(1,2,
        widths=(4/8,4/8)),
        size=(950,400),
        margin=5mm
    )

    # wind stress sensitivity

    wind_stress_derivative = heatmap(LinRange(0, 3840, 127),
    LinRange(0, 3840, 128),
    dS.forcing.Fx',
    c=:balance,
    xlabel="x (km)",
    xguidefontsize=13,
    ylabel="y (km)",
    yguidefontsize=13,
    title=L"\partial J / F_x",
    plot_titlefontsize=13,
    colorbar_title=L"m",
    colorbar_titlefontsize=13,
    colorbar=:false,
    clim=(-50000,50000),
    dpi=300,
    size=(500,500)
    )

    # bottom drag coefficient sensitivity 
    heatmap(LinRange(0, 3840, 127),
        LinRange(0, 3840, 128),
        dS.constants.cDu[2:end-1,2:end-1]',
        c=:balance,
        xlabel="x (km)",
        xguidefontsize=13,
        ylabel="y (km)",
        yguidefontsize=13,
        title=L"\partial J / c_D^u(x,y)",
        plot_titlefontsize=13,
        colorbar_title=L"m",
        colorbar_titlefontsize=13,
        dpi=300,
        size=(500,500)
    )

    heatmap(LinRange(0, 3840, 128),
        LinRange(0, 3840, 127),
        dS.constants.cDv[2:end-1,2:end-1]',
        c=:balance,
        xlabel="x (km)",
        xguidefontsize=13,
        ylabel="y (km)",
        yguidefontsize=13,
        title=L"\partial J / c_D^v(x,y)",
        plot_titlefontsize=13,
        colorbar_title=L"m",
        colorbar_titlefontsize=13,
        dpi=300,
        size=(500,500),
        clim=(-1e6,1e6),
        colorbar=:false
    )

end

function create_adjoint_gif()

    primal_fid = h5open("primal_technicalpaper_5000m_period286_1yearintegration_110124.h5")
    adj_fid = h5open("adjoint_technicalpaper_5000m_period286_1yearintegration_110124.h5", "r")
    # states = ncread("../data_files_gamma0.3/1024_spinup/eta.nc", "eta")

    u_anim = Animation()
    v_anim = Animation()
    eta_anim = Animation()

    for j = 81797:-286:1
    # for j = 1:3651


        blob = read(adj_fid[string(j)])
        adj_chkp = deserialize(blob)

        temp = ShallowWaters.PrognosticVars{Float32}(
            ShallowWaters.remove_halo(adj_chkp.Prog.u,
            adj_chkp.Prog.v,
            adj_chkp.Prog.η,
            adj_chkp.Prog.sst,adj_chkp)...
        )

        frame(eta_anim, heatmap(temp.η',
            # title=L"\partial \mathcal{E}(t_f)/\partial u(%$j)",
            title=L"\partial \mathcal{E} / \partial \eta(%$j)",
            xlabel=L"x",
            ylabel=L"y",
            c=:balance,
            dpi=300,
            colorbar_title="         ",
            xlabelfontsize=14,
            ylabelfontsize=14,
            clim=(-15000,15000),
            colorbar_titlefontsize=14,
            colorbar_tickfontsize=8,
            xtickfontsize=8,
            ytickfontsize=8)
        )

        # frame(v_anim, heatmap(temp.v',
        #     title=L"\partial \mathcal{E}(t_f)/\partial v(%$j)",
        #     clim=(-1e9, 1e9),
        #     legend=:none,
        #     xlabel=L"x",
        #     ylabel=L"y",
        #     c=:balance,
        #     dpi=300,
        #     colorbar_title="         ",
        #     xlabelfontsize=14,
        #     ylabelfontsize=14,
        #     colorbar_titlefontsize=14,
        #     colorbar_tickfontsize=8,
        #     xtickfontsize=8,
        #     ytickfontsize=8)
        # )

        # frame(eta_anim, heatmap(temp.η',
        #     title=L"\partial \mathcal{E}(t_f)/\partial \eta(%$j)",
        #     clim=(-50, 50),
        #     legend=:none,
        #     xlabel=L"x",
        #     ylabel=L"y",
        #     c=:balance,
        #     dpi=300,
        #     colorbar_title="         ",
        #     xlabelfontsize=14,
        #     ylabelfontsize=14,
        #     colorbar_titlefontsize=14,
        #     colorbar_tickfontsize=8,
        #     xtickfontsize=8,
        #     ytickfontsize=8)
        # )

    end

    gif(eta_anim, "eta_integration2_110424.gif", fps = 8)
    # gif(u_anim, "du_integration_365_energy_withclosure_fps7_031424.png", fps = 7)
    # gif(v_anim, "dv_integration_365_energy_withclosure_fps7_031424.png", fps = 7)


end

"""
The last few lines are about running the above functions
"""

diffs, enzyme_deriv, S, dS = check_derivative(1*365, 500)

# @save "technicalpaper_finalprimal_struct_500mdepth_1year_nocDfield_110424.jld2" S
# @save "technicalpaper_finaladjoint_struct_500mdepth_1year_nocDfield_110424.jld2" dS
@save "technicalpaper_fdcheck_vector_500mdepth_1year_nocDfield_110424.jld2" diffs

# create_adjoint_gif()