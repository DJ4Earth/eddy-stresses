# this is just to save 10 days worth of coarse-grained high resolution data for the 
# online problems

include("eddy_paper.jl")

function hourly_save_run(S_true)

    ufiltered = zeros(1023, 1024, 240)
    vfiltered = zeros(1024, 1023, 240)
    etafiltered = zeros(1024, 1024, 240)

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    # setup that happens prior to the actual integration
    Diag = S_true.Diag
    Prog = S_true.Prog

    @unpack u,v,η,sst = Prog
    @unpack u0,v0,η0 = Diag.RungeKutta
    @unpack u1,v1,η1 = Diag.RungeKutta
    @unpack du,dv,dη = Diag.Tendencies
    @unpack du_sum,dv_sum,dη_sum = Diag.Tendencies
    @unpack du_comp,dv_comp,dη_comp = Diag.Tendencies

    @unpack um,vm = Diag.SemiLagrange

    @unpack dynamics,RKo,RKs,tracer_advection = S_true.parameters
    @unpack time_scheme,compensated = S_true.parameters
    @unpack RKaΔt,RKbΔt = S_true.constants
    @unpack Δt_Δ,Δt_Δs = S_true.constants

    @unpack nt,dtint = S_true.grid
    @unpack nstep_advcor,nstep_diff,nadvstep,nadvstep_half = S_true.grid

    # calculate layer thicknesses for initial conditions
    ShallowWaters.thickness!(Diag.VolumeFluxes.h,η,S_true.forcing.H)
    ShallowWaters.Ix!(Diag.VolumeFluxes.h_u,Diag.VolumeFluxes.h)
    ShallowWaters.Iy!(Diag.VolumeFluxes.h_v,Diag.VolumeFluxes.h)
    ShallowWaters.Ixy!(Diag.Vorticity.h_q,Diag.VolumeFluxes.h)

    # calculate PV terms for initial conditions
    urhs = Diag.PrognosticVarsRHS.u .= u
    vrhs = Diag.PrognosticVarsRHS.v .= v
    ηrhs = Diag.PrognosticVarsRHS.η .= η

    ShallowWaters.advection_coriolis!(urhs,vrhs,ηrhs,Diag,S_true)
    ShallowWaters.PVadvection!(Diag,S_true)

    # propagate initial conditions
    copyto!(u0,u)
    copyto!(v0,v)
    copyto!(η0,η)

    # store initial conditions of sst for relaxation
    copyto!(Diag.SemiLagrange.sst_ref,sst)
    j = 1
    t = 0           # model time
    for i = 1:S_true.grid.nt

        Diag = S_true.Diag
        Prog = S_true.Prog
    
        @unpack u,v,η,sst = Prog
        @unpack u0,v0,η0 = Diag.RungeKutta
        @unpack u1,v1,η1 = Diag.RungeKutta
        @unpack du,dv,dη = Diag.Tendencies
        @unpack du_sum,dv_sum,dη_sum = Diag.Tendencies
        @unpack du_comp,dv_comp,dη_comp = Diag.Tendencies
    
        @unpack um,vm = Diag.SemiLagrange
    
        @unpack dynamics,RKo,RKs,tracer_advection = S_true.parameters
        @unpack time_scheme,compensated = S_true.parameters
        @unpack RKaΔt,RKbΔt = S_true.constants
        @unpack Δt_Δ,Δt_Δs = S_true.constants
    
        @unpack nt,dtint = S_true.grid
        @unpack nstep_advcor,nstep_diff,nadvstep,nadvstep_half = S_true.grid

        # ghost point copy for boundary conditions
        ShallowWaters.ghost_points!(u,v,η,S_true)
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
                ShallowWaters.ghost_points!(u1,v1,η1,S_true)
            end

            # type conversion for mixed precision
            u1rhs = Diag.PrognosticVarsRHS.u .= u1
            v1rhs = Diag.PrognosticVarsRHS.v .= v1
            η1rhs = Diag.PrognosticVarsRHS.η .= η1

            ShallowWaters.rhs!(u1rhs,v1rhs,η1rhs,Diag,S_true,t)          # momentum only
            ShallowWaters.continuity!(u1rhs,v1rhs,η1rhs,Diag,S_true,t)   # continuity equation

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

        ShallowWaters.ghost_points!(u0,v0,η0,S_true)

        # type conversion for mixed precision
        u0rhs = Diag.PrognosticVarsRHS.u .= u0
        v0rhs = Diag.PrognosticVarsRHS.v .= v0
        η0rhs = Diag.PrognosticVarsRHS.η .= η0

        # ADVECTION and CORIOLIS TERMS
        # although included in the tendency of every RK substep,
        # only update every nstep_advcor steps if nstep_advcor > 0
        if dynamics == "nonlinear" && nstep_advcor > 0 && (i % nstep_advcor) == 0
            ShallowWaters.UVfluxes!(u0rhs,v0rhs,η0rhs,Diag,S_true)
            ShallowWaters.advection_coriolis!(u0rhs,v0rhs,η0rhs,Diag,S_true)
        end

        # DIFFUSIVE TERMS - SEMI-IMPLICIT EULER
        # use u0 = u^(n+1) to evaluate tendencies, add to u0 = u^n + rhs
        # evaluate only every nstep_diff time steps
        if (i % nstep_diff) == 0
            ShallowWaters.bottom_drag!(u0rhs,v0rhs,η0rhs,Diag,S_true)
            ShallowWaters.diffusion!(u0rhs,v0rhs,Diag,S_true)
            ShallowWaters.add_drag_diff_tendencies!(u0,v0,Diag,S_true)
            ShallowWaters.ghost_points_uv!(u0,v0,S_true)
        end

        # TRACER ADVECTION
        u0rhs = Diag.PrognosticVarsRHS.u .= u0
        v0rhs = Diag.PrognosticVarsRHS.v .= v0
        ShallowWaters.tracer!(i,u0rhs,v0rhs,Prog,Diag,S_true)

        if i ∈ 75:75:S_true.grid.nt
            temp = ShallowWaters.PrognosticVars{S_true.parameters.Tprog}(
                ShallowWaters.remove_halo(u,v,η,sst,S_true)...)
            ufiltered[:,:,j] .= imfilter(temp.u, reflect(ker))
            vfiltered[:,:,j] .= imfilter(temp.v, reflect(ker))
            etafiltered[:,:,j] .= imfilter(temp.η, reflect(ker))
            j+=1
        end

        # Copy back from substeps
        copyto!(u,u0)
        copyto!(v,v0)
        copyto!(η,η0)

        t += dtint

    end

    return [ufiltered, vfiltered, etafiltered]

end


function run()

    T = Float64
    Ndays = 10

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
        N=1,
        Ndays=10,
        α=2,
        nx=1024,
        initial_cond="ncfile",
        initpath="./spinup_files/1024_spinup_noslip/"
    )

    S = ShallowWaters.model_setup(P)

    filtered_states = hourly_save_run(S)

    jldsave("coarsegrained_hrstates_uveta_10days_imfilter_102825.jld2", filtered_states=filtered_states)

    return nothing

end

function downsize()

    cgstates = load_object("./coarsegrained_hrstates_uveta_10days_imfilter_102825.jld2")

    ucg = cgstates[1]
    vcg = cgstates[2]
    etacg = cgstates[3]

    ucgdownsized = (ucg[8:8:end, 4:8:end, :] .+ ucg[8:8:end, 5:8:end, :]) ./ 2
    vcgdownsize = (vcg[4:8:end, 8:8:end, :] .+ vcg[5:8:end, 8:8:end, :]) ./ 2
    etacgdownsized = (etacg[4:8:end,4:8:end,:] .+ etacg[5:8:end,5:8:end,:] .+ etacg[4:8:end,5:8:end,:] .+ etacg[5:8:end,4:8:end,:]) ./ 4

    # the following is to create the coarse grained uv term for computing the off-diagonal entries in T

    uhrh = cat(zeros(T,1023+2*halo,halo),cat(zeros(T,halo,1024),hru[:,:,1],zeros(T,halo,1024),dims=1),zeros(T,1023+2*halo,halo),dims=2)
    vhrh = cat(zeros(T,1024+2*halo,halo),cat(zeros(T,halo,1023),hrv[:,:,1],zeros(T,halo,1023),dims=1),zeros(T,1024+2*halo,halo),dims=2)

    # moving to hr corner grid and cut off the halo

    uhrq = ShallowWaters.Iy(uhrh)[2:end-1,2:end-1]
    vhrq = ShallowWaters.Ix(vhrh)[2:end-1,2:end-1]

    uhrT = zeros(1024,1024)
    vhrT = zeros(1024,1024)

    ShallowWaters.Ixy!(uhrT,uhrq)
    ShallowWaters.Ixy!(vhrT,vhrq)

    uvbar = ShallowWaters.coarse_grain_eta(uhrT .* vhrT, 1024, SZB)

end

# computing and saving true tensors for hourly data
function hourly_Ts()

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    T11filtered = zeros(1024, 1024, 240)
    T22filtered = zeros(1024, 1024, 240)
    T12filtered = zeros(1025, 1025, 240)

    cgstates = load_object("./offline_files/coarsegrained_hrstates_uveta_10days_imfilter_102825.jld2")
    

    ufiltered[:,:,j] .= imfilter(temp.u, reflect(ker))
    vfiltered[:,:,j] .= imfilter(temp.v, reflect(ker))
    etafiltered[:,:,j] .= imfilter(temp.η, reflect(ker))


end