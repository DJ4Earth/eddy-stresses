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

        if S.parameters.i in 1:224:S.grid.nt

            temp = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(S.Prog.u,
            S.Prog.v,
            S.Prog.η,
            S.Prog.sst,S)...)

            energy_lr = (sum(temp.u.^2) + sum(temp.v.^2)) / (S.grid.nx * S.grid.ny)

            S.parameters.J += energy_lr

        end
        #############################################

        # Copy back from substeps
        copyto!(u,u0)
        copyto!(v,v0)
        copyto!(η,η0)

    end

    ##### use if time-averaging the objective function ########
    S.parameters.J = S.parameters.J / length(1:224:S.grid.nt)
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
        write_checkpoints_filename = "technicalpaper_timeavgobj_1000m_1year_111424.h5",
        write_checkpoints_period = 224
    )

    autodiff(Enzyme.ReverseWithPrimal, checkpointed_integration, Duplicated(S, dS), Const(revolve))

    enzyme_deriv = dS.constants.cDfield[50,50]

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

        S_inner.constants.cDfield[50,50] += s

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
    
    # investigating derivatives 

    # norm of the prognostic variables derivatives
    detanorm = load_object("./technicalpaper_timeaveragedobjective_detanorm_dividedbynxny_1000m_1year_111424.jld2");
    dunorm = load_object("./technicalpaper_timeaveragedobjective_dunorm_dividedbynxny_1000m_1year_111424.jld2");
    dvnorm = load_object("./technicalpaper_timeaveragedobjective_dvnorm_dividedbynxny_1000m_1year_111424.jld2");

    timestep = 365:-0.99726:1
    f = Figure();
    ax1 = Axis(f[1, 1],
        title = "Norm of adjoint derivative w.r.t. u, v, H = 1000m, time-averaged objective",
        xlabel = "t (days)",
        ylabel = L"||\partial J / \partial x||",
    )
    ax2 = Axis(f[2, 1],
    title = "Norm of adjoint derivative w.r.t. eta, H = 1000m, time-averaged objective",
    xlabel = "t (days)",
    ylabel = L"||\partial J / \partial x||",

    )
    lines!(ax2, timestep, detanorm, label = L"\partial J / \partial \eta(t)")
    axislegend(ax2, position = :rt)
    lines!(ax1, timestep, dunorm, label = L"\partial J / \partial u(t)")
    lines!(ax1, timestep, dvnorm, label = L"\partial J / \partial v(t)")
    axislegend(ax1, position = :rt)

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

    adj500 = h5open("adjoint_technicalpaper_check_500m_1year_110624.h5.h5", "r")
    blob = read(adj500["1"])
    adj = deserialize(blob)

    dS = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(adj.Prog.u,
    adj.Prog.v,
    adj.Prog.η,
    adj.Prog.sst,adj)...)

    adj5002 = load_object("./technicalpaper_datafiles/timeaveragedobjective/technicalpaper_timeavgobj_finaladjoint_struct_500mdepth_1year_111324.jld2")
    dS500m1 = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(adj5002.Prog.u,
    adj5002.Prog.v,
    adj5002.Prog.η, 
    adj5002.Prog.sst,adj5002)...)

    adj2 = load_object("./technicalpaper_datafiles/nontimeaveragedobjective/technicalpaper_finaladjoint_struct_500mdepth_1year_correctedcDfield_110624.jld2")
    dS2 = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(adj2.Prog.u,
    adj2.Prog.v,
    adj2.Prog.η, 
    adj2.Prog.sst,adj2)...)


    # Initial condition

    derivs500 = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(dS500m1.Prog.u,
    dS500m1.Prog.v,
    dS500m1.Prog.η, 
    dS500m1.Prog.sst,dS500m1)...)

    derivs5km = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(dS5km1.Prog.u,
    dS5km1.Prog.v,
    dS5km1.Prog.η, 
    dS5km1.Prog.sst,dS5km1)...)

    fig = Figure(size=(2200, 1500));

    ax1 = Axis(fig[1,1], 
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Initial x-velocity sensitivity, H = 500m, time-averaged objective function",
    );
    hm1 = CairoMakie.heatmap!(ax1, 0:127:3840,
        0:128:3840,
        dS3.u,
        colormap=:balance,
        colorrange = (-5e5, 5e5)
    );
    Colorbar(fig[2, 1], 
        hm1,
        vertical=false
    )

    ax2 = Axis(fig[1,2],
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Initial y-velocity sensitivity, H = 500m, time-averaged objective function",
    );
    hm2 = CairoMakie.heatmap!(ax2, 0:128:3840,
        0:128:3840,
        dS3.v,
        colormap=:balance,
        colorrange = (-5e5, 5e5)
    );
    Colorbar(fig[2, 2], hm2, vertical=false)

    ax3 = Axis(fig[3,1], 
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Initial x-velocity sensitivity, H = 500m, non time-averaged objective function",
    );
    hm1 = CairoMakie.heatmap!(ax3, 0:127:3840,
        0:128:3840,
        dS2.u,
        colormap=:balance,
        colorrange = (-5e19, 5e19)
    );
    Colorbar(fig[4, 1], 
        hm1,
        vertical=false
    )

    ax4 = Axis(fig[3,2],
        xlabel = "x (km)",
        ylabel = "y (km)",
        title = "Initial y-velocity sensitivity, H = 500m, non time-averaged objective function",
    );
    hm2 = CairoMakie.heatmap!(ax4, 0:128:3840,
        0:128:3840,
        dS2.v,
        colormap=:balance,
        colorrange = (-5e19, 5e19)
    );
    Colorbar(fig[4, 2], hm2, vertical=false)

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


end

function create_adjoint_gif()

    # primal_fid = h5open("technicalpaper.h5")
    adj_fid = h5open("./adjoint_technicalpaper_timeavgobj_1000m_1year_111424.h5", "r")
    # states = ncread("../data_files_gamma0.3/1024_spinup/eta.nc", "eta")

    # unorm_anim = Animation()
    # vnorm_anim = Animation()
    # etanorm_anim = Animation()

    dunorm = []
    dvnorm = []
    detanorm = []

    dcDnorm = []
    dFxnorm = []


    # for j = 81797:-286:1
    # for j = 1:3651
    for j = 81761:-224:1


        blob = read(adj_fid[string(j)])
        adj_chkp = deserialize(blob)

        temp = ShallowWaters.PrognosticVars{Float32}(
            ShallowWaters.remove_halo(adj_chkp.Prog.u,
            adj_chkp.Prog.v,
            adj_chkp.Prog.η,
            adj_chkp.Prog.sst,adj_chkp)...
        )

        push!(dcDnorm, sum(adj_chkp.constants.cDfield.^2) / 128^2)
        push!(dFxnorm, sum(adj_chkp.forcing.Fx.^2) / (127 * 128))
        push!(dunorm, sum(temp.u.^2) / (127*128))
        push!(dvnorm, sum(temp.v.^2) / (127*128))
        push!(detanorm, sum(temp.η.^2) / (128 * 128))

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

    @save "technicalpaper_timeaveragedobjective_dcDnorm_dividedbynxny_1000m_1year_111424.jld2" dcDnorm
    @save "technicalpaper_timeaveragedobjective_dFxnorm_dividedbynxny_1000m_1year_1114224.jld2" dFxnorm
    @save "technicalpaper_timeaveragedobjective_dunorm_dividedbynxny_1000m_1year_111424.jld2" dunorm
    @save "technicalpaper_timeaveragedobjective_dvnorm_dividedbynxny_1000m_1year_111424.jld2" dvnorm
    @save "technicalpaper_timeaveragedobjective_detanorm_dividedbynxny_1000m_1year_111424.jld2" detanorm

end

"""
The last few lines are about running the above functions
"""

diffs, enzyme_deriv, S, dS = run_adjoint_plusfd(
    output=true,
    L_ratio=1,
    g=9.81,
    H=1000,
    wind_forcing_x="double_gyre",
    Lx=3840e3,
    seasonal_wind_x=false,
    topography="flat",
    bc="nonperiodic",
    bottom_drag="quadratic",
    α=2,
    nx=128,
    Ndays=365,
    initial_cond="ncfile",
    initpath="./data_files_gamma0.3/128_spinup_noforcing_noslipbc/"
)
# @save "technicalpaper_timeavgobj_finalprimal_struct_1000mdepth_1year_111424.jld2" S
# @save "technicalpaper_timeavgobj_finaladjoint_struct_1000mdepth_1year_111424.jld2" dS
# @save "technicalpaper_timeavgobj_fdcheck_vector_1000mdepth_1year_11424.jld2" diffs


# create_adjoint_gif()