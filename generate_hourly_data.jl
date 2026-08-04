"""
this script is outdated and no longer really used
"""
include("eddy_paper.jl")

function filter()

    ker = ImageFiltering.Kernel.gaussian((30e3/3750, 30e3/3750), (65, 65))

    ufiltered = zeros(1023, 1024, 1096);
    vfiltered = zeros(1024, 1023, 1096);
    # etafiltered = zeros(1024, 1024, 1096);

    for j = 1:1096
        ufiltered[:,:,j] .= imfilter(uhrall[:,:,j], reflect(ker))
        vfiltered[:,:,j] .= imfilter(vhrall[:,:,j], reflect(ker))
        # etafiltered[:,:,j] .= imfilter(etahrall[:,:,j], reflect(ker),border=NA())
    end

    jldsave("1024_filtered_uv_imfilter_3years_postspinup_dailysaves_largerkernel.jld2", uv = [ufiltered, vfiltered])

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

end

function downsize()

    # cgstates = load_object("./dissipation_constant/spinup_files/1024_filtered_uveta_imfilter_7years_startfrom3yearpostspinup_weeklysaves.jld2");

    # ucg = cgstates[1];
    # vcg = cgstates[2];
    # etacg = cgstates[3];

    ucg = ufiltered;
    vcg = vfiltered;
    etacg = etafiltered;

    ucgdownsized = (ucg[8:8:end, 4:8:end,:] .+ ucg[8:8:end, 5:8:end,:]) ./ 2;
    vcgdownsized = (vcg[4:8:end, 8:8:end,:] .+ vcg[5:8:end, 8:8:end,:]) ./ 2;
    etacgdownsized = (etacg[4:8:end,4:8:end,:] .+ etacg[5:8:end,5:8:end,:] .+ etacg[4:8:end,5:8:end,:] .+ etacg[5:8:end,4:8:end,:]) ./ 4;

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
