# this is just to save 10 days worth of coarse-grained high resolution data for the 
# online problems

include("eddy_paper.jl")

function hourly_save_run(S_true)

    hrstates = []
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

        # if i ∈ 8*75:8*75:S_true.grid.nt
        if i ∈ 225:224:S_true.grid.nt
            temp = ShallowWaters.PrognosticVars{S_true.parameters.Tprog}(
                ShallowWaters.remove_halo(u,v,η,sst,S_true)...)
            push!(hrstates, temp)
            # ufiltered[:,:,j] .= imfilter(temp.u, reflect(ker))
            # vfiltered[:,:,j] .= imfilter(temp.v, reflect(ker))
            # etafiltered[:,:,j] .= imfilter(temp.η, reflect(ker))
            j+=1
        end

        # Copy back from substeps
        copyto!(u,u0)
        copyto!(v,v0)
        copyto!(η,η0)

        t += dtint

    end

    return hrstates
    # return [ufiltered, vfiltered, etafiltered]

end

function filter()

    # u1 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day1-766saves/u.nc", "u");
    # v1 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day1-766saves/v.nc", "v");
    # eta1 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day1-766saves/eta.nc", "eta");

    # u2 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day766-end/u.nc", "u");
    # v2 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day766-end/v.nc", "v");
    # eta2 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day766-end/eta.nc", "eta");

    # u = ncread("./dissipation_constant/spinup_files/1024_7years_startfrom3yearpostspinup_weeklysaves/u.nc", "u");
    # v = ncread("./dissipation_constant/spinup_files/1024_7years_startfrom3yearpostspinup_weeklysaves/v.nc", "v");
    # eta = ncread("./dissipation_constant/spinup_files/1024_7years_startfrom3yearpostspinup_weeklysaves/eta.nc", "eta");

    u1 = ncread("./dissipation_constant/spinup_files/1024_spinup_newwindamp/1024_spinup_newwindamp_days1-374/u.nc", "u");
    v1 = ncread("./dissipation_constant/spinup_files/1024_spinup_newwindamp/1024_spinup_newwindamp_days1-374/v.nc", "v");
    eta1 = ncread("./dissipation_constant/spinup_files/1024_spinup_newwindamp/1024_spinup_newwindamp_days1-374/eta.nc", "eta");

    u2 = ncread("./dissipation_constant/spinup_files/1024_spinup_newlatandwindamp/u.nc", "u");
    v2 = ncread("./dissipation_constant/spinup_files/1024_spinup_newlatandwindamp/v.nc", "v");
    eta2 = ncread("./dissipation_constant/spinup_files/1024_spinup_newlatandwindamp/eta.nc", "eta");

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    ufiltered = zeros(1023, 1024, 366)
    vfiltered = zeros(1024, 1023, 366)
    etafiltered = zeros(1024, 1024, 366)

    for j = 1:381
        ufiltered[:,:,j] .= imfilter(u[:,:,j], reflect(ker))
        vfiltered[:,:,j] .= imfilter(v[:,:,j], reflect(ker))
        etafiltered[:,:,j] .= imfilter(eta[:,:,j], reflect(ker))
    end

    for j = 2:end
        ufiltered[:,:,j+766] .= imfilter(u2[:,:,j+1], reflect(ker))
        vfiltered[:,:,j+766] .= imfilter(v2[:,:,j+1], reflect(ker))
        etafiltered[:,:,j+766] .= imfilter(eta2[:,:,j+1], reflect(ker))
    end

    jldsave("1024_filtered_uveta_imfilter_3years_postspinup_dailysaves_correctedsetup.jld2", uveta = [ufiltered, vfiltered, etafiltered])

end

function run()

    T = Float64
    Ndays = 10

    P = ShallowWaters.Parameter(T=T;
        output=true,
        output_dt=8,
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
        Ndays=Ndays,
        α=2,
        nx=1024,
        initial_cond="ncfile",
        initpath="./spinup_files/1024_spinup_noslip/"
    )

    S = ShallowWaters.model_setup(P)

    filtered_states = hourly_save_run(S)

    jldsave("coarsegrained_hrstates_uveta_10days_hourlysaves_imfilter_111025.jld2", filtered_states=filtered_states)

    return nothing

    ∈
end

function downsize()

    cgstates = load_object("./dissipation_constant/spinup_files/1024_filtered_uveta_imfilter_7years_startfrom3yearpostspinup_weeklysaves.jld2");

    ucg = cgstates[1];
    vcg = cgstates[2];
    etacg = cgstates[3];

    ucgdownsized = (ucg[8:8:end, 4:8:end] .+ ucg[8:8:end, 5:8:end]) ./ 2;
    vcgdownsized = (vcg[4:8:end, 8:8:end] .+ vcg[5:8:end, 8:8:end]) ./ 2;
    etacgdownsized = (etacg[4:8:end,4:8:end] .+ etacg[5:8:end,5:8:end] .+ etacg[4:8:end,5:8:end] .+ etacg[5:8:end,4:8:end]) ./ 4;

    jldsave("1024_filtered_downsized_uveta_imfilter_7years_startfrom3yearpostspinup_weeklysaves.jld2", uveta = [ucgdownsized, vcgdownsized, etacgdownsized])

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

    T = Float64
    S_true = ShallowWaters.model_setup(T=T; output=false,
        output_dt = 8,
        L_ratio=1,
        g=9.81,
        H=500,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        diffusion="Smagorinsky",        # this is the only new parameter to be adjusted in the new spinups
        tracer_advection=false,
        tracer_relaxation=false,
        N=1,
        α=2,
        nx=1024,
        Ndays=30,
        initial_cond="rest",
        initpath="./spinup_files/1024_spinup_noslip"
    );
    halo = S_true.grid.halo

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    T11true = zeros(1024, 1024, 273)
    T22true = zeros(1024, 1024, 273)
    T12true = zeros(1025, 1025, 273)

    u = ncread("./dissipation_smagorinsky/spinup_files_newdissipation/1024_1yearpostspinup_smag_noslipbc_8hoursaves/u.nc", "u")
    v = ncread("./dissipation_smagorinsky/spinup_files_newdissipation/1024_1yearpostspinup_smag_noslipbc_8hoursaves/v.nc", "v")
    eta = ncread("./dissipation_smagorinsky/spinup_files_newdissipation/1024_1yearpostspinup_smag_noslipbc_8hoursaves/eta.nc", "eta")

    # cgstates = load_object("./offline_files/hrstates_filtered_downsized_hourly_tendays_uveta_beginsatonehour_111925.jld2")

    for j = 1:273

        # the following is to create the coarse grained uv term for computing the off-diagonal entries in T

        uhrh = cat(zeros(T,1023+2*halo,halo),cat(zeros(T,halo,1024),u[:,:,j],zeros(T,halo,1024),dims=1),zeros(T,1023+2*halo,halo),dims=2)
        vhrh = cat(zeros(T,1024+2*halo,halo),cat(zeros(T,halo,1023),v[:,:,j],zeros(T,halo,1023),dims=1),zeros(T,1024+2*halo,halo),dims=2)

        # moving to hr corner grid and cut off the halo

        uhrq = ShallowWaters.Iy(uhrh)[2:end-1,2:end-1]
        vhrq = ShallowWaters.Ix(vhrh)[2:end-1,2:end-1]

        uhrT = zeros(1024,1024)
        vhrT = zeros(1024,1024)

        ShallowWaters.Ixy!(uhrT,uhrq)
        ShallowWaters.Ixy!(vhrT,vhrq)

        # ubar = uhrT
        # vbar = vhrT

        ubar = imfilter(uhrT, reflect(ker))
        vbar = imfilter(vhrT, reflect(ker))

        # usqbar = ubar.^2
        # vsqbar = vbar.^2

        usqbar = imfilter(uhrT.^2, reflect(ker))
        vsqbar = imfilter(vhrT.^2, reflect(ker))

        # uvbar = uhrq .* vhrq
        uvbar = imfilter(uhrq .* vhrq, reflect(ker))
        ubarvbar = imfilter(uhrq, reflect(ker)) .* imfilter(vhrq, reflect(ker))

        T11true[:,:,j] .= ubar .* ubar - usqbar
        T22true[:,:,j] .= vbar .* vbar - vsqbar
        T12true[:,:,j] .= ubarvbar - uvbar

    end

    jldsave("trueTs_90days_8hoursaves_T11T22T12.jld2", Ts=[T11true, T22true, T12true])

    return T11true, T22true, T12true

end

function coarsen_Ts()

    Ts = load_object("./dissipation_smagorinsky/spinup_files_newdissipation/trueTs_90days_8hoursaves_T11T22T12.jld2")

    T11 = Ts[1]
    T22 = Ts[2]
    T12 = Ts[3]

    T11downsized = zeros(128, 128, 273)
    T22downsized = zeros(128, 128, 273)
    T12downsized = zeros(129, 129, 273)
    for j = 1:273
        T11downsized[:,:,j] = (T11[4:8:end,4:8:end,j] .+ T11[5:8:end,5:8:end,j] .+ T11[4:8:end,5:8:end,j] .+ T11[5:8:end,4:8:end,j]) ./ 4
        T22downsized[:,:,j] = (T22[4:8:end,4:8:end,j] .+ T22[5:8:end,5:8:end,j] .+ T22[4:8:end,5:8:end,j] .+ T22[5:8:end,4:8:end,j]) ./ 4
        T12downsized[:,:,j] = T12[1:8:end,1:8:end,j]
    end

    jldsave("trueTs_downsized_90days_8hoursaves_T11T22T12.jld2", Ts=[T11downsized, T22downsized, T12downsized])

end

# prior way I computed the true Ts, hopefully the sameish as above

# ucg = load_object("./offline_files/coarsegrained_hr_ubar_ubarsq_t1_foroffline_101525.jld2")
# vcg = load_object("./offline_files/coarsegrained_hr_vbar_vbarsq_t1_foroffline_101525.jld2")
# uvbar = load_object("./offline_files/coarsegrained_hr_uvbar_t1_foroffline_101525.jld2")

# ubar = ucg[1]
# ubarsq = ucg[2]

# vbar = vcg[1]
# vbarsq = vcg[2]

# ubarh = cat(zeros(T,nux+2*halo,halo),cat(zeros(T,halo,nuy),ubar,zeros(T,halo,nuy),dims=1),zeros(T,nux+2*halo,halo),dims=2)
# ubarhsq = cat(zeros(T,nux+2*halo,halo),cat(zeros(T,halo,nuy),ubarsq,zeros(T,halo,nuy),dims=1),zeros(T,nux+2*halo,halo),dims=2)

# vbarh = cat(zeros(T,nvx+2*halo,halo),cat(zeros(T,halo,nvy),vbar,zeros(T,halo,nvy),dims=1),zeros(T,nvx+2*halo,halo),dims=2)
# vbarhsq = cat(zeros(T,nvx+2*halo,halo),cat(zeros(T,halo,nvy),vbarsq,zeros(T,halo,nvy),dims=1),zeros(T,nvx+2*halo,halo),dims=2)

# uvbarh = cat(zeros(T,nx+2*haloη,haloη),cat(zeros(T,haloη,ny),uvbar,zeros(T,haloη,ny),dims=1),zeros(T,nx+2*haloη,haloη),dims=2)

# T12_true = ShallowWaters.Iy(ubarh)[2:end-1,2:end-1] .* ShallowWaters.Ix(vbarh)[2:end-1,2:end-1] - ShallowWaters.Ixy(uvbarh)

# T11_true = (ShallowWaters.Ixy(ShallowWaters.Iy(ubarh)[2:end-1,2:end-1])).^2 - ShallowWaters.Ixy(ShallowWaters.Iy(ubarhsq)[2:end-1,2:end-1])
# T22_true = ShallowWaters.Ixy((ShallowWaters.Ix(vbarh)[2:end-1,2:end-1])).^2 - ShallowWaters.Ixy(ShallowWaters.Ix(vbarhsq)[2:end-1,2:end-1])
