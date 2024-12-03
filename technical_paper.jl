# Assuming the information we have is solely a spacially averaged kinetic energy

include("../ShallowWaters.jl/src/ShallowWaters.jl")
using .ShallowWaters

using Enzyme
using Checkpointing, HDF5, Serialization
using NetCDF, JLD2, CairoMakie

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

    S.parameters.data = zeros(S.grid.nt + 1)

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

        #### Energy objective function, time averaged

        if S.parameters.i in (S.grid.nt - 30*224):1:S.grid.nt

            temp = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(S.Prog.u,
            S.Prog.v,
            S.Prog.η,
            S.Prog.sst,S)...)

            energy_lr = (sum(temp.u.^2) + sum(temp.v.^2)) / (S.grid.nx * S.grid.ny)

            S.parameters.J += energy_lr

        end
        #############################################

        # Storing the energy over time
        temp = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(S.Prog.u,
        S.Prog.v,
        S.Prog.η,
        S.Prog.sst,S)...)

        S.parameters.data[S.parameters.i] = (sum(temp.u.^2) + sum(temp.v.^2)) / (S.grid.nx * S.grid.ny)

        # Copy back from substeps
        copyto!(u,u0)
        copyto!(v,v0)
        copyto!(η,η0)

    end

    ##### use if time-averaging the objective function #######
    S.parameters.J = S.parameters.J / length((S.grid.nt - 30*224):1:S.grid.nt)
    ##########################################################

    ##### Energy objective function, not time averaged ###########
    # temp = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(S.Prog.u,
    # S.Prog.v,
    # S.Prog.η,
    # S.Prog.sst,S)...)

    # energy_lr = (sum(temp.u.^2) + sum(temp.v.^2)) / (S.grid.nx * S.grid.ny)
    # S.parameters.J = energy_lr
    ###########################################

    return nothing

end

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
        write_checkpoints=true,
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

"""
Mostly figure generation, I just wanted to be able to run include("technical_paper.jl")
without all of this also running. The function deserialize opens any saved checkpoints
"""
function deserialize(x)
    s = IOBuffer(x)
    Serialization.deserialize(s)
end

function stuff()

    f = Figure()
    ax = Axis(f[1,1])
    h = heatmap!(ax, 0:128:3840,
    0:127:3840,
    temp.u,
    colormap=:balance,
    colorrange=(-maximum(temp.u), maximum(temp.u))
    )
    Colorbar(f[1,2], h)
    
    # investigating derivatives 

    # norm of the prognostic variables derivatives
    detanorm = load_object("./technicalpaper_datafiles/finaldatafiles/tech")
    dunorm = load_object("./technicalpaper_timeaveragedobjective_starting3monthsin_dunorm_dividedbynxny_500m_1year_111424.jld2")
    dvnorm = load_object("./technicalpaper_timeaveragedobjective_starting3monthsin_dvnorm_dividedbynxny_500m_1year_111424.jld2")

    timestep1 = 60:-0.99726:1
    f = Figure(size = (1000, 500));
    ax1 = Axis(f[1, 1],
        title = "Norm of adjoint derivative, two-month integration from spinup",
        xlabel = "t (days)",
        ylabel = L"||\partial J / \partial x||",
    )
    timestep2 = 120:-0.99726:1
    ax2 = Axis(f[2,1],
    title = "Norm of adjoint derivative, four-month integration from spinup",
    xlabel = "t (days)",
    # yscale=log10,
    ylabel = L"||\partial J / \partial x||",
    )
    timestep3 = 360:-0.99726:1
    ax3 = Axis(f[3, 1],
    title = "Norm of adjoint derivative, twelve-month integration from spinup",
    xlabel = "t (days)",
    # yscale = log,
    ylabel = L"||\partial J / \partial x||")


    lines!(ax1, timestep1, dunorm2, label = L"\partial J / \partial u(t)")
    lines!(ax1, timestep1, dvnorm2, label = L"\partial J / \partial v(t)")
    axislegend(ax1, position = :rt)
    lines!(ax2, timestep2, dunorm4, label = L"\partial J / \partial u(t)")
    lines!(ax2, timestep2, dvnorm4, label = L"\partial J / \partial v(t)")
    axislegend(ax2, position = :rt)
    lines!(ax3, timestep3, dunorm12, label = L"\partial J / \partial u(t)")
    lines!(ax3, timestep3, dvnorm12, label = L"\partial J / \partial v(t)")
    axislegend(ax3, position = :rt)

    save("technicalpaper_normprog_divnxny_1000m_114224.png", f)

    # norm of the forcing and parameter derivatives

    dcDnorm = load_object("./technicalpaper_timeaveragedobjective_dcDnorm_dividedbynxny_1000m_1year_111424.jld2");
    dFxnorm = load_object("./technicalpaper_timeaveragedobjective_dFxnorm_dividedbynxny_1000m_1year_1114224.jld2");

    timestep = 365:-0.99726:1
    f = Figure();
    ax1 = Axis(f[1, 1],
        title = "Norm of adjoint derivative w.r.t. wind-stress, time-averaged objective, H = 1000m",
        xlabel = "t (days)",
        ylabel = L"||\partial J / \partial F_x||"
    );
    lines!(ax1, timestep, dFxnorm, label = L"\partial J / \partial F_x")
    axislegend(ax1, position = :rt);

    ax2 = Axis(f[2, 1],
    title = "Norm of adjoint derivative w.r.t. bottom drag, time-averaged objective, H = 1000m",
    xlabel = "t (days)",
    ylabel = L"||\partial J / \partial c_D||");
    lines!(ax2, timestep, dcDnorm, label = L"\partial J / \partial c_D");
    axislegend(ax2, position = :rt);


    save("technicalpaper_parametersensitivities_normcDFx_divnxny_1000mdepth_111424.png", f)

    # computing time averaged velocity/variance/mean flow

    u = ncread("./30km_1year_output/u.nc", "u")
    v = ncread("./30km_1year_output/v.nc", "v")

    avg_u = zeros(127,128)
    avg_v = zeros(128,127)
    for t = 1:366
        avg_u += u[:, :, t]
        avg_v += v[:, :, t]
    end

    f = Figure()
    ax1 = Axis(f[1, 1],
    title = "Average x-velocity",
    xlabel = "x",
    ylabel = "y"
    )
    heatmap!(ax1, 0:127:3840,
        0:128:3840,
        avg_u,
        colorrange = (-150, 150),
        colormap=:balance
    )

    f = Figure()
    ax2 = Axis(f[1, 1],
    title = "Time averaged y-velocity",
    xlabel = "x",
    ylabel = "y"
    )
    h = heatmap!(ax2, 0:128:3840,
    0:127:3840,
    avg_v,
    colorrange = (-600, 600),
    colormap=:balance
    )
    Colorbar(f[1,2], h)
    save("technicalpaper_avgv_111324.png")

    # loss function is final spatially averaged energy
    # initial condition sensitivity

    primal500 = h5open("primal_technicalpaper_check_110524.h5.h5", "r")
    blob = read(primal500["1"])
    prim = deserialize(blob)

    adj50012 = h5open("./technicalpaper_datafiles/finaldatafiles/adjoint_technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_500m_12months_float32_112024.h5", "r")
    day = read(adj50012["225"])
    month = read(adj50012["6721"])
    sixmonth = read(adj50012["40321"])
    twelvemonth = read(adj50012["80641"])

    adj1 = deserialize(day)
    adj30 = deserialize(month)
    adj180 = deserialize(sixmonth)
    adj360 = deserialize(twelvemonth)

    fig = Figure(size=(900, 900));

    dS360 = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(adj360.Prog.u,
    adj360.Prog.v,
    adj360.Prog.η, 
    adj360.Prog.sst,adj360)...)
    ax1 = Axis(fig[1,1], 
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Objective function sensitivity to u -0 month",
    );
    hm1 = CairoMakie.heatmap!(ax1, 0:127:3840,
        0:128:3840,
        dS360.u,
        colormap=:balance,
        colorrange = (-maximum(dS360.u), maximum(dS360.u))
    );
    Colorbar(fig[2, 1], 
        hm1,
        vertical=false
    )

    dS180 = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(adj180.Prog.u,
    adj180.Prog.v,
    adj180.Prog.η, 
    adj180.Prog.sst,adj180)...)
    ax2 = Axis(fig[1,2],
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Objective function sensitivity to u -6 month",
    );
    hm2 = CairoMakie.heatmap!(ax2, 0:128:3840,
        0:128:3840,
        dS180.u,
        colormap=:balance,
        colorrange = (-maximum(dS180.u),maximum(dS180.u))
    );
    Colorbar(fig[2, 2], hm2, vertical=false)

    dS30 = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(adj30.Prog.u,
    adj30.Prog.v,
    adj30.Prog.η, 
    adj30.Prog.sst,adj30)...)
    ax3 = Axis(fig[3,1], 
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Objective function sensitivity to u -11 month",
    );
    hm1 = CairoMakie.heatmap!(ax3, 0:127:3840,
        0:128:3840,
        dS30.u,
        colormap=:balance,
        colorrange = (-maximum(dS30.u),maximum(dS30.u))
    );
    Colorbar(fig[4, 1], 
        hm1,
        vertical=false
    )

    dS1 = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(adj1.Prog.u,
    adj1.Prog.v,
    adj1.Prog.η, 
    adj1.Prog.sst,adj1)...)
    ax4 = Axis(fig[3,2],
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Objective function sensitivity to u -12 month",
    );
    hm2 = CairoMakie.heatmap!(ax4, 0:128:3840,
        0:128:3840,
        dS1.u,
        colormap=:balance,
        colorrange = (-maximum(dS1.u), maximum(dS1.u))
    );
    Colorbar(fig[4, 2], hm2, vertical=false)

    fig

    adj50012 = load_object("./technicalpaper_datafiles/finaldatafiles/technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_startingfromspinup_finaladjoint_struct_500mdepth_12months_float32start_112024.jld2")
    dS50012 = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(adj50012.Prog.u,
    adj50012.Prog.v,
    adj50012.Prog.η, 
    adj50012.Prog.sst,adj50012)...)

    # Initial condition

    fig = Figure(size=(900, 450));

    ax1 = Axis(fig[1,1], 
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = L"\partial J / \partial u(t_0, x, y)",
    );
    hm1 = CairoMakie.heatmap!(ax1, 0:127:3840,
        0:128:3840,
        dS50012.u,
        colormap=:balance,
        colorrange = (-maximum(dS50012.u), maximum(dS50012.u))
    );
    Colorbar(fig[2, 1], 
        hm1,
        vertical=false
    )

    ax2 = Axis(fig[1,2],
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = L"\partial J / \partial v(t_0, x, y)",
    );
    hm2 = CairoMakie.heatmap!(ax2, 0:128:3840,
        0:128:3840,
        dS50012.v,
        colormap=:balance,
        colorrange = (-maximum(dS50012.v),maximum(dS50012.v))
    );
    Colorbar(fig[2, 2], hm2, vertical=false)
    
    ax3 = Axis(fig[1,1], 
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = L"\partial J / \partial F_x",
    );
    hm1 = CairoMakie.heatmap!(ax3, 0:127:3840,
        0:128:3840,
        adj50012.forcing.Fx,
        colormap=:balance,
        colorrange = (-maximum(adj50012.forcing.Fx), maximum(adj50012.forcing.Fx))
    );
    Colorbar(fig[1, 2], 
        hm1,
        vertical=true
    )

    ax4 = Axis(fig[1,1],
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = L"\partial J / \partial c_D",
    );
    hm2 = CairoMakie.heatmap!(ax4, 0:128:3840,
        0:128:3840,
        adj50012.constants.cDfield[2:end-1,2:end-1],
        colormap=:balance,
        colorrange = (-maximum(adj50012.constants.cDfield), maximum(adj50012.constants.cDfield))
    );
    Colorbar(fig[1, 2], hm2, vertical=true)

    fig

    # wind stress sensitivity
    # Makie plot
    fig = Figure(size = (1000,500));
    ax1 = Axis(fig[1,1], 
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Wind-stress sensitivity, H = 500m, time-averaged objective"
    );
    hm1 = CairoMakie.heatmap!(ax1, 0:128:3840,
        0:128:3840,
        adj3.forcing.Fx,
        colormap=:balance,
        colorrange=(-20,20)
    );
    Colorbar(fig[:, end+1], hm1);

    ax2 = Axis(fig[1,3], 
    xlabel = "x (km)",
    ylabel = "y (km)",
    title = "Wind stress sensitivity, H = 500m, non time-averaged objective"
    );
    hm2 = CairoMakie.heatmap!(ax2, 0:128:3840,
    0:128:3840,
    adj2.forcing.Fx,
    colorrange = (-3e8, 3e8),
    colormap=:balance,
    );
    Colorbar(fig[1, 4], hm2);
    fig


    # bottom drag coefficient sensitivity
    fig = Figure(size=(1000,500));
    ax1 = Axis(fig[1,1], 
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Bottom drag sensitivity, H = 500m, time-averaged objective"
    );
    hm1 = CairoMakie.heatmap!(ax1, 0:128:3840,
        0:128:3840,
        adj3.constants.cDfield[2:end-1,2:end-1],
        colorrange = (-600, 600),
        colormap=:balance,
    );
    Colorbar(fig[1, 2], hm1);

    ax2 = Axis(fig[1,3], 
    xlabel = "x (km)",
    ylabel = "y (km)",
    title = "Bottom drag sensitivity, H = 500m, non time-averaged objective"
    );
    hm2 = CairoMakie.heatmap!(ax2, 0:128:3840,
    0:128:3840,
    adj2.constants.cDfield[2:end-1,2:end-1],
    colorrange = (-4e9, 4e9),
    colormap=:balance,
    );
    Colorbar(fig[1, 4], hm2);
    fig

    #### Energy plot 

    t = LinRange(0, 11, 888686)
    t1 = LinRange(0, 1, 80640)
    t2 = LinRange(10, 11, 80641)
    fig = Figure(size=(1000,700));
    ax1 = Axis(fig[1,1], title="Energy during 10-year spinup + 1 year",
    ylabel="Energy")
    ax2 = Axis(fig[2,1], title="Energy during first year of spinup", ylabel="Energy")
    ax3 = Axis(fig[3,1], title="Energy 1-year beyond spinup", ylabel="Energy", xlabel="Years")
    lines!(ax1, t, energy)
    lines!(ax1, t1, energy[1:360*224], label="First year")
    lines!(ax1, t2, energy[808046:end], label="One year beyond spinup")
    axislegend(ax1, position=:rb)
    lines!(ax2, t1, energy[1:360*224])
    lines!(ax3,t2, energy[808046:end])


    fig = Figure();
    ax = Axis(fig[1,1])
    lines!(ax, energy)
    lines!(ax, energy2)
    lines!(ax, S3.parameters.data[738328:end-1])

end

function create_adjoint_gif()

    # primal_fid = h5open("technicalpaper.h5")
    adj_fid = h5open("./technicalpaper_datafiles/finaldatafiles/adjoint_technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_500m_12months_float32_112024.h5")
    # primal_fid = h5open("./primal_technicalpaper_timeavgobj_500m_2months_111524.h5")
    # states = ncread("../data_files_gamma0.3/1024_spinup/eta.nc", "eta")

    # unorm_anim = Animation()
    # vnorm_anim = Animation()
    # etanorm_anim = Animation()

    dunorm4 = []
    dvnorm4 = []
    detanorm4 = []

    dcDnorm4 = []
    dFxnorm4 = []

    J = []

    final = 1:224:224*30*12
    for j = final[end]:-224:1
    # for j = 1:3651
    # for j = 1:224:81761


        blob = read(adj_fid[string(j)])
        adj_chkp = deserialize(blob)

        # push!(J, primal_chkp.parameters.J)


        temp = ShallowWaters.PrognosticVars{Float32}(
            ShallowWaters.remove_halo(adj_chkp.Prog.u,
            adj_chkp.Prog.v,
            adj_chkp.Prog.η,
            adj_chkp.Prog.sst,adj_chkp)...
        )

        push!(dcDnorm4, sum(adj_chkp.constants.cDfield.^2) / 128^2)
        push!(dFxnorm4, sum(adj_chkp.forcing.Fx.^2) / (127 * 128))
        push!(dunorm4, sum(temp.u.^2) / (127*128))
        push!(dvnorm4, sum(temp.v.^2) / (127*128))
        push!(detanorm4, sum(temp.η.^2) / (128 * 128))

        # frame(eta_anim, heatmap(temp.η',
        #     # title=L"\partial \mathcal{E}(t_f)/\partial u(%$j)",
        #     title=L"\partial \mathcal{E} / \partial \eta(%$j)",
        #     xlabel=L"x",
        #     ylabel=L"y",
        #     c=:balance,
        #     dpi=300,
        #     colorbar_title="         ",
        #     xlabelfontsize=14,
        #     ylabelfontsize=14,
        #     clim=(-15000,15000),
        #     colorbar_titlefontsize=14,
        #     colorbar_tickfontsize=8,
        #     xtickfontsize=8,
        #     ytickfontsize=8)
        # )

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

    # gif(eta_anim, "eta_integration2_110424.gif", fps = 8)
    # gif(u_anim, "du_integration_365_energy_withclosure_fps7_031424.png", fps = 7)
    # gif(v_anim, "dv_integration_365_energy_withclosure_fps7_031424.png", fps = 7)

    @save "technicalpaper_timeaveragedobjective_finalmonthobj_dcDnormrest_500m_12months_111824.jld2" dcDnorm12rest
    @save "technicalpaper_timeaveragedobjective_finalmonthobj_dFxnormrest_500m_12months_111824.jld2" dFxnorm12rest
    @save "technicalpaper_timeaveragedobjective_finalmonthobj_dunormrest_500m_12months_111824.jld2" dunorm12rest
    @save "technicalpaper_timeaveragedobjective_finalmonthobj_dvnormrest_500m_12months_111824.jld2" dvnorm12rest
    @save "technicalpaper_timeaveragedobjective_finalmonthobj_detanormrest_500m_12months_111824.jld2" detanorm12rest

end

function energy_plot()

    f = Figure();
    ax = Axis(f[1,1], xlabel="Timestep", ylabel="Energy")
    for j = 1:5

        S = ShallowWaters.model_setup(output=false,
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
        nx=128,
        Ndays=12*30,
        initial_cond="ncfile",
        initpath="./data_files_gamma0.3/10yearspinup_128_noslipbc_fromrest_float32params"
        )

        snaps = Int(floor(sqrt(S.grid.nt)))
        revolve = Revolve{ShallowWaters.ModelSetup}(S.grid.nt, snaps;
            verbose=1,
            gc=true,
            write_checkpoints=false
        )

        S.Prog.u = S.Prog.u + 0.001 .* randn(Float32, 131, 132)
        S.Prog.v = S.Prog.v + 0.001 .* randn(Float32, 132, 131)
        S.Prog.η = S.Prog.η + 0.001 .* randn(Float32, 130, 130)

        _ = checkpointed_integration(S, revolve)

        lines!(ax, S.parameters.data[1:end-1])
        hlines!(ax, sum(S.parameters.data) / S.grid.nt)

    end

end

"""
The last few lines are about running the above functions
"""

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
    Ndays=10,
    initial_cond="ncfile",
    initpath="./data_files_gamma0.3/10yearspinup_128_noslipbc_fromrest_float32params"
)

# @save "technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_startingfromspinup_finalprimal_struct_500mdepth_4months_float32start_112124.jld2" S
# @save "technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_startingfromspinup_finaladjoint_struct_500mdepth_4months_float32start_112124.jld2" dS
# @save "technicalpaper_timeavgobj_onlyfinalmonth_everytimestep_startingfromspinup_fdcheck_vector_500mdepth_4months_float32start_112124.jld2" diffs


# create_adjoint_gif()