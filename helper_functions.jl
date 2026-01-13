
function compute_hrS(u, v)

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
        tracer_advection=false,
        tracer_relaxation=false,
        N=1,
        α=2,
        nx=1024,
        Ndays=30,
        initial_cond="ncfile",
        initpath="./dissipation_constant/spinup_files/1024_spinup_noslip"
    );
    halo = S_true.grid.halo

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    T11true = zeros(1024, 1024)
    T22true = zeros(1024, 1024)
    T12true = zeros(1025, 1025)

    # cgstates = load_object("./offline_files/hrstates_filtered_downsized_hourly_tendays_uveta_beginsatonehour_111925.jld2")
    uhrh = cat(zeros(T,1023+2*halo,halo),cat(zeros(T,halo,1024),u,zeros(T,halo,1024),dims=1),zeros(T,1023+2*halo,halo),dims=2)
    vhrh = cat(zeros(T,1024+2*halo,halo),cat(zeros(T,halo,1023),v,zeros(T,halo,1023),dims=1),zeros(T,1024+2*halo,halo),dims=2)

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

    T11true[:,:] .= ubar .* ubar - usqbar
    T22true[:,:] .= vbar .* vbar - vsqbar
    T12true[:,:] .= ubarvbar - uvbar

    T11downsized = zeros(128, 128)
    T22downsized = zeros(128, 128)
    T12downsized = zeros(129, 129)

    T11downsized[:,:] .= (T11true[4:8:end,4:8:end] .+ T11true[5:8:end,5:8:end] .+ T11true[4:8:end,5:8:end] .+ T11true[5:8:end,4:8:end]) ./ 4
    T22downsized[:,:] .= (T22true[4:8:end,4:8:end] .+ T22true[5:8:end,5:8:end] .+ T22true[4:8:end,5:8:end] .+ T22true[5:8:end,4:8:end]) ./ 4
    T12downsized[:,:] .= T12true[1:8:end,1:8:end]
    
    Slr = ShallowWaters.model_setup(T=T; output=false,
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
        tracer_advection=false,
        tracer_relaxation=false,
        N=1,
        α=2,
        nx=128,
        Ndays=30
    );
    halo = Slr.grid.halo
    dT11dx = zeros(T,Slr.Diag.CNNVars.nux,Slr.Diag.CNNVars.nuy)         # derivative of T11 in the x-direction, u-grid
    dT12dy = zeros(T,Slr.Diag.CNNVars.nux+halo,Slr.Diag.CNNVars.nuy)    # derivative of T12 in the y-direction, u-grid
    dT12dx = zeros(T,Slr.Diag.CNNVars.nvx,Slr.Diag.CNNVars.nvy+halo)    # derivative of T12 in the x-direction, v-grid
    dT22dy = zeros(T,Slr.Diag.CNNVars.nvx,Slr.Diag.CNNVars.nvy)         # derivative of T22 in the y-direction, v-grid

    ShallowWaters.∂x!(dT11dx, T11downsized)
    ShallowWaters.∂y!(dT12dy, T12downsized)

    ShallowWaters.∂x!(dT12dx, T12downsized)
    ShallowWaters.∂y!(dT22dy, T22downsized)

    S_u = zeros(T,Slr.Diag.CNNVars.nux,Slr.Diag.CNNVars.nuy)             # total forcing in x-direction
    S_v = zeros(T,Slr.Diag.CNNVars.nvx,Slr.Diag.CNNVars.nvy)             # total forcing in y-direction

    @inbounds for j in 1:Slr.Diag.CNNVars.nuy
        for k in 1:Slr.Diag.CNNVars.nux
            S_u[k,j] = Slr.grid.scale * (dT11dx[k,j] + dT12dy[k+1,j])
        end
    end

    @inbounds for j in 1:Slr.Diag.CNNVars.nvy
        for k in 1:Slr.Diag.CNNVars.nvx
            S_v[k,j] = Slr.grid.scale * (dT22dy[k,j] + dT12dx[k,j+1])
        end
    end

    return S_u, S_v

end

function compute_true_hrS()

    uhr = ncread("./dissipation_constant/spinup_files/1024_postspinup_noslip_5years_061824/u.nc", "u");
    vhr = ncread("./dissipation_constant/spinup_files/1024_postspinup_noslip_5years_061824/v.nc", "v");
    etahr = ncread("./dissipation_constant/spinup_files/1024_postspinup_noslip_5years_061824/eta.nc", "eta");

    coarse_grained_hrstates = load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_3years_postspinup_dailysaves.jld2");
    uhrcg = coarse_grained_hrstates[1];
    vhrcg = coarse_grained_hrstates[2];
    etahrcg = coarse_grained_hrstates[3];

    Phr = ShallowWaters.Parameter(T=Float64,
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
        α=2,
        nx=1024,
        Ndays=1
    )
    Shr = ShallowWaters.model_setup(Phr);

    Plr = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        cfl=.898,
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
        α=2,
        nx=128,
        Ndays=3*365
    );
    Slr = ShallowWaters.model_setup(Plr);

    # this is where the total time derivative gets stored
    # u0,v0,η0 = Diag.RungeKutta

    S_u = zeros(127,128,1097)
    S_v = zeros(128,127,1097)

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))
    for n = 1:1097

        uhr_, vhr_, etahr_ = ShallowWaters.add_halo(uhr[:,:,n+1], vhr[:,:,n+1], etahr[:,:,n+1], Shr)
        uhrcg_, vhrcg_, etahrcg_ = ShallowWaters.add_halo(uhrcg[:,:,n+1], vhrcg[:,:,n+1], etahrcg[:,:,n+1], Slr)

        Shr_ = deepcopy(Shr)
        Slr_ = deepcopy(Slr)

        Shr_.Prog.u = uhr_
        Shr_.Prog.v = vhr_
        Shr_.Prog.η = etahr_

        Slr_.Prog.u = uhrcg_
        Slr_.Prog.v = vhrcg_
        Slr_.Prog.η = etahrcg_

        duhr, dvhr, detahr = single_step(Shr_, n*1800*Shr_.grid.dtint)
        duhrcg, dvhrcg, detahrcg = single_step(Slr_, n*225*Slr_.grid.dtint)

        ufiltered = imfilter(duhr, reflect(ker))
        vfiltered = imfilter(dvhr, reflect(ker))

        udownsized = (ufiltered[8:8:end, 4:8:end, :] .+ ufiltered[8:8:end, 5:8:end, :]) ./ 2;
        vdownsized = (vfiltered[4:8:end, 8:8:end, :] .+ vfiltered[5:8:end, 8:8:end, :]) ./ 2;

        S_u[:,:,n] .= duhrcg - udownsized
        S_v[:,:,n] .= dvhrcg - vdownsized

    end

end

function paddingu(x, Su, n1)

    T = Real
    nfft = nextfastfft(size(x))
    out = zeros(DSP.Util.fftabs2type(T), n1>>1 + 1)

    if prod(nfft) == length(x) && isa(x, StridedArray)
        input1 = x[:,:,t] # no need to pad
        input2 = Su
    else
        input1 = zeros(fftintype(T), nfft)
        input1[1:size(x,1), 1:size(x,2)] = x
        input2 = zeros(fftintype(T), nfft)
        input2[1:size(Su,1), 1:size(Su,2)] = Su
    end

    return out, input1, input2

end

function paddingv(x, Sv, n1)

    T = Real
    nfft = nextfastfft(size(x))
    out = zeros(DSP.Util.fftabs2type(T), n1>>1 + 1)

    if prod(nfft) == length(x) && isa(x, StridedArray)
        input1 = x # no need to pad
        input2 = Sv
    else
        input1 = zeros(fftintype(T), nfft)
        input1[1:size(x,1), 1:size(x,2)] = x
        input2 = zeros(fftintype(T), nfft)
        input2[1:size(Sv,1), 1:size(Sv,2)] = Sv
    end

    return out, input1, input2, nfft

end

function fft2pow2radial!(out::Array{T}, s_fft::Matrix{Complex{T}}, u_fft::Matrix{Complex{T}}, n1, n2; ptype=2) where T
    r = length(s_fft)
    nmin = min(n1, n2)

    n1max = n1 >> 1 + 1  # since rfft is used
    n1max != size(s_fft, 1) && throw(ArgumentError("fft size incorrect"))
    m1 = convert(T, 1/r)
    m2 = convert(T, 2/r)
    wavenum = 0          # wavenumber index
    kmax = length(out)   # the highest wavenumber
    wc = zeros(Int, kmax) # wave count for radial average
    if n1 == nmin        # scale the wavevector for non-square s_fft
        c2 = n1/n2
        c1 = one(c2)
    else
        c1 = n2/n1
        c2 = one(c1)
    end

    sqrt_muladd(a, k) = sqrt(muladd(a, a, k))

    @inbounds begin
        for j = 1:n2
            kj1 = ifelse(j <= n2>>1 + 1, j-1, -n2+j-1)
            kj2 = (kj1 * c2)^2

            wavenum = round(Int, sqrt_muladd(c1 * (1 - 1), kj2)) + 1
            if wavenum<=kmax
                out[wavenum] = muladd(real(conj(s_fft[1, j]) * u_fft[1, j]), m1, out[wavenum])
                wc[wavenum] += 1
            end
            for i = 2:n1max-1
                wavenum = round(Int, sqrt_muladd(c1 * (i - 1), kj2)) + 1
                if wavenum<=kmax
                    out[wavenum] = muladd(real(conj(s_fft[i, j]) * u_fft[i, j]), m2, out[wavenum])
                    wc[wavenum] += 2
                end
            end
            wavenum = round(Int, sqrt_muladd(c1 * (n1max - 1), kj2)) + 1
            if wavenum<=kmax
                out[wavenum] = muladd(real(conj(s_fft[n1max, j]) * u_fft[n1max, j]), ifelse(iseven(n1), m1, m2), out[wavenum])
                wc[wavenum] += ifelse(iseven(n1), 1, 2)
            end
        end
    end
    if ptype == 2
        for i = 1:kmax
            @inbounds out[i] /= wc[i]
        end
    end
    out
end

function single_step(S, t)

    # calculate layer thicknesses for initial conditions
    ShallowWaters.thickness!(S.Diag.VolumeFluxes.h, S.Prog.η, S.forcing.H)
    ShallowWaters.Ix!(S.Diag.VolumeFluxes.h_u, S.Diag.VolumeFluxes.h)
    ShallowWaters.Iy!(S.Diag.VolumeFluxes.h_v, S.Diag.VolumeFluxes.h)
    ShallowWaters.Ixy!(S.Diag.Vorticity.h_q, S.Diag.VolumeFluxes.h)

    # calculate PV terms for initial conditions
    urhs = S.Diag.PrognosticVarsRHS.u .= S.Prog.u
    vrhs = S.Diag.PrognosticVarsRHS.v .= S.Prog.v
    ηrhs = S.Diag.PrognosticVarsRHS.η .= S.Prog.η

    ShallowWaters.advection_coriolis!(urhs, vrhs, ηrhs, S.Diag, S)
    ShallowWaters.PVadvection!(S.Diag, S)

    # propagate initial conditions
    copyto!(S.Diag.RungeKutta.u0, S.Prog.u)
    copyto!(S.Diag.RungeKutta.v0, S.Prog.v)
    copyto!(S.Diag.RungeKutta.η0, S.Prog.η)

    # store initial conditions of sst for relaxation
    copyto!(S.Diag.SemiLagrange.sst_ref, S.Prog.sst)

    # run a single step of integration loop

    # ghost point copy for boundary conditions
    ShallowWaters.ghost_points!(S.Prog.u, S.Prog.v, S.Prog.η, S)
    copyto!(S.Diag.RungeKutta.u1, S.Prog.u)
    copyto!(S.Diag.RungeKutta.v1, S.Prog.v)
    copyto!(S.Diag.RungeKutta.η1, S.Prog.η)

    if S.parameters.compensated
        fill!(S.Diag.Tendencies.du_sum, zero(S.parameters.Tprog))
        fill!(S.Diag.Tendencies.dv_sum, zero(S.parameters.Tprog))
        fill!(S.Diag.Tendencies.dη_sum, zero(S.parameters.Tprog))
    end

    for rki = 1:S.parameters.RKo
        if rki > 1
            ShallowWaters.ghost_points!(
                S.Diag.RungeKutta.u1,
                S.Diag.RungeKutta.v1,
                S.Diag.RungeKutta.η1,
                S
            )
        end

        # type conversion for mixed precision
        u1rhs = S.Diag.PrognosticVarsRHS.u .= S.Diag.RungeKutta.u1
        v1rhs = S.Diag.PrognosticVarsRHS.v .= S.Diag.RungeKutta.v1
        η1rhs = S.Diag.PrognosticVarsRHS.η .= S.Diag.RungeKutta.η1

        ShallowWaters.rhs!(u1rhs, v1rhs, η1rhs, S.Diag, S, t)          # momentum only
        ShallowWaters.continuity!(u1rhs, v1rhs, η1rhs, S.Diag, S, t)   # continuity equation

        if rki < S.parameters.RKo
            ShallowWaters.caxb!(
                S.Diag.RungeKutta.u1,
                S.Prog.u,
                S.constants.RKbΔt[rki],
                S.Diag.Tendencies.du
            )
            ShallowWaters.caxb!(
                S.Diag.RungeKutta.v1,
                S.Prog.v,
                S.constants.RKbΔt[rki],
                S.Diag.Tendencies.dv
            )
            ShallowWaters.caxb!(
                S.Diag.RungeKutta.η1,
                S.Prog.η,
                S.constants.RKbΔt[rki],
                S.Diag.Tendencies.dη
            )
        end

        if S.parameters.compensated
            ShallowWaters.axb!(S.Diag.Tendencies.du_sum, S.constants.RKaΔt[rki], S.Diag.Tendencies.du)
            ShallowWaters.axb!(S.Diag.Tendencies.dv_sum, S.constants.RKaΔt[rki], S.Diag.Tendencies.dv)
            ShallowWaters.axb!(S.Diag.Tendencies.dη_sum, S.constants.RKaΔt[rki], S.Diag.Tendencies.dη)
        else
            ShallowWaters.axb!(
                S.Diag.RungeKutta.u0,
                S.constants.RKaΔt[rki],
                S.Diag.Tendencies.du
            )
            ShallowWaters.axb!(
                S.Diag.RungeKutta.v0,
                S.constants.RKaΔt[rki],
                S.Diag.Tendencies.dv
            )
            ShallowWaters.axb!(
                S.Diag.RungeKutta.η0,
                S.constants.RKaΔt[rki],
                S.Diag.Tendencies.dη
            )
        end
    end

    if S.parameters.compensated
        ShallowWaters.axb!(S.Diag.Tendencies.du_sum, -1, S.Diag.Tendencies.du_comp)
        ShallowWaters.axb!(S.Diag.Tendencies.dv_sum, -1, S.Diag.Tendencies.dv_comp)
        ShallowWaters.axb!(S.Diag.Tendencies.dη_sum, -1, S.Diag.Tendencies.dη_comp)

        ShallowWaters.axb!(S.Diag.RungeKutta.u0, 1, S.Diag.Tendencies.du_sum)
        ShallowWaters.axb!(S.Diag.RungeKutta.v0, 1, S.Diag.Tendencies.dv_sum)
        ShallowWaters.axb!(S.Diag.RungeKutta.η0, 1, S.Diag.Tendencies.dη_sum)

        ShallowWaters.dambmc!(
            S.Diag.Tendencies.du_comp,
            S.Diag.RungeKutta.u0,
            S.Prog.u,
            S.Diag.Tendencies.du_sum
        )
        ShallowWaters.dambmc!(
            S.Diag.Tendencies.dv_comp,
            S.Diag.RungeKutta.v0,
            S.Prog.v,
            S.Diag.Tendencies.dv_sum
        )
        ShallowWaters.dambmc!(
            S.Diag.Tendencies.dη_comp,
            S.Diag.RungeKutta.η0,
            S.Prog.η,
            S.Diag.Tendencies.dη_sum
        )
    end

    ShallowWaters.ghost_points!(
        S.Diag.RungeKutta.u0,
        S.Diag.RungeKutta.v0,
        S.Diag.RungeKutta.η0,
        S
    )

    u0rhs = S.Diag.PrognosticVarsRHS.u .= S.Diag.RungeKutta.u0
    v0rhs = S.Diag.PrognosticVarsRHS.v .= S.Diag.RungeKutta.v0
    η0rhs = S.Diag.PrognosticVarsRHS.η .= S.Diag.RungeKutta.η0

    # if S.parameters.dynamics == "nonlinear" && S.grid.nstep_advcor > 0 && (i % S.grid.nstep_advcor) == 0
        ShallowWaters.UVfluxes!(u0rhs, v0rhs, η0rhs, S.Diag, S)
        ShallowWaters.advection_coriolis!(u0rhs, v0rhs, η0rhs, S.Diag, S)
    # end

    # if (chkp.i % S.grid.nstep_diff) == 0
        ShallowWaters.bottom_drag!(u0rhs, v0rhs, η0rhs, S.Diag, S)
        ShallowWaters.diffusion!(u0rhs, v0rhs, S.Diag, S)
        ShallowWaters.add_drag_diff_tendencies!(
            S.Diag.RungeKutta.u0,
            S.Diag.RungeKutta.v0,
            S.Diag,
            S
        )
        ShallowWaters.ghost_points_uv!(
            S.Diag.RungeKutta.u0,
            S.Diag.RungeKutta.v0,
            S
        )
    # end

    t += S.grid.dtint

    P = ShallowWaters.PrognosticVars{S.parameters.Tprog}(ShallowWaters.remove_halo(u0rhs,
        v0rhs,
        η0rhs,
        S.Prog.sst,
        S)...
    )

    # u0rhs = S.Diag.PrognosticVarsRHS.u .= S.Diag.RungeKutta.u0
    # v0rhs = S.Diag.PrognosticVarsRHS.v .= S.Diag.RungeKutta.v0

    return P.u, P.v, P.η

end