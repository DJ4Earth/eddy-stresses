# computes the nonlinear advection term *approximation* using high-resolution snapshots
function compute_approxhrS()

    T = Float64
    Shr = ShallowWaters.model_setup(T=T; output=false,
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
        Ndays=1
    );

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
    halo = Shr.grid.halo

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    halo = Slr.grid.halo
    dT11dx = zeros(T,Slr.Diag.CNNVars.nux,Slr.Diag.CNNVars.nuy)         # derivative of T11 in the x-direction, u-grid
    dT12dy = zeros(T,Slr.Diag.CNNVars.nux+halo,Slr.Diag.CNNVars.nuy)    # derivative of T12 in the y-direction, u-grid
    dT12dx = zeros(T,Slr.Diag.CNNVars.nvx,Slr.Diag.CNNVars.nvy+halo)    # derivative of T12 in the x-direction, v-grid
    dT22dy = zeros(T,Slr.Diag.CNNVars.nvx,Slr.Diag.CNNVars.nvy)         # derivative of T22 in the y-direction, v-grid

    S_u = zeros(T,Slr.Diag.CNNVars.nux,Slr.Diag.CNNVars.nuy, 365)             # total forcing in x-direction
    S_v = zeros(T,Slr.Diag.CNNVars.nvx,Slr.Diag.CNNVars.nvy, 365)             # total forcing in y-direction

    uhrT = zeros(1024,1024)
    vhrT = zeros(1024,1024)

    s = Slr.grid.Δ^2
    for t = 1:365

        uhrh = cat(zeros(T,1023+2*halo,halo),cat(zeros(T,halo,1024),uhrall[:,:,t+1096],zeros(T,halo,1024),dims=1),zeros(T,1023+2*halo,halo),dims=2)
        vhrh = cat(zeros(T,1024+2*halo,halo),cat(zeros(T,halo,1023),vhrall[:,:,t+1096],zeros(T,halo,1023),dims=1),zeros(T,1024+2*halo,halo),dims=2)

        # moving to hr corner grid and cut off the halo

        uhrq = ShallowWaters.Iy(uhrh)[2:end-1,2:end-1]
        vhrq = ShallowWaters.Ix(vhrh)[2:end-1,2:end-1]

        ShallowWaters.Ixy!(uhrT,uhrq)
        ShallowWaters.Ixy!(vhrT,vhrq)

        ubar = imfilter(uhrT, reflect(ker))
        vbar = imfilter(vhrT, reflect(ker))

        usqbar = imfilter(uhrT.^2, reflect(ker))
        vsqbar = imfilter(vhrT.^2, reflect(ker))

        # uvbar = uhrq .* vhrq
        uvbar = imfilter(uhrq .* vhrq, reflect(ker))
        ubarvbar = imfilter(uhrq, reflect(ker)) .* imfilter(vhrq, reflect(ker))

        T11true = ubar .* ubar - usqbar
        T22true = vbar .* vbar - vsqbar
        T12true = ubarvbar - uvbar

        T11downsized = (T11true[4:8:end,4:8:end] .+ T11true[5:8:end,5:8:end] .+ T11true[4:8:end,5:8:end] .+ T11true[5:8:end,4:8:end]) ./ 4
        T22downsized = (T22true[4:8:end,4:8:end] .+ T22true[5:8:end,5:8:end] .+ T22true[4:8:end,5:8:end] .+ T22true[5:8:end,4:8:end]) ./ 4
        T12downsized = T12true[1:8:end,1:8:end]

        ShallowWaters.∂x!(dT11dx, T11downsized)
        ShallowWaters.∂y!(dT12dy, T12downsized)

        ShallowWaters.∂x!(dT12dx, T12downsized)
        ShallowWaters.∂y!(dT22dy, T22downsized)

        @inbounds for j in 1:Slr.Diag.CNNVars.nuy
            for k in 1:Slr.Diag.CNNVars.nux
                S_u[k,j,t] = (dT11dx[k,j] + dT12dy[k+1,j]) / s
            end
        end

        @inbounds for j in 1:Slr.Diag.CNNVars.nvy
            for k in 1:Slr.Diag.CNNVars.nvx
                S_v[k,j,t] = (dT22dy[k,j] + dT12dx[k,j+1]) / s
            end
        end

    end

    return S_u, S_v

end

# used to compute the SGS term S from a difference of time-derivatives (total tendencies)
function compute_true_hrS()

    uhr1 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day1-766saves/u.nc", "u");
    vhr1 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day1-766saves/v.nc", "v");
    etahr1 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day1-766saves/eta.nc", "eta");

    uhr2 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day766-end/u.nc", "u");
    vhr2 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day766-end/v.nc", "v");
    etahr2 = ncread("./dissipation_constant/spinup_files/1024_postspinup_3years_dailysaves_correctedsetup/1024_postspinup_day766-end/eta.nc", "eta");

    coarse_grained_hrstates = load_object("./dissipation_constant/spinup_files/1024_filtered_downsized_uveta_imfilter_3years_postspinup_dailysaves_correctedsetup.jld2");
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
    );
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

    for S in [Slr, Shr]
        # calculate layer thicknesses for initial conditions
        ShallowWaters.thickness!(S.Diag.VolumeFluxes.h, S.Prog.η, S.forcing.H)
        ShallowWaters.Ix!(S.Diag.VolumeFluxes.h_u, S.Diag.VolumeFluxes.h)
        ShallowWaters.Iy!(S.Diag.VolumeFluxes.h_v, S.Diag.VolumeFluxes.h)
        ShallowWaters.Ixy!(S.Diag.Vorticity.h_q, S.Diag.VolumeFluxes.h)
    end

    # this is where the total time derivative gets stored
    # u0,v0,η0 = Diag.RungeKutta

    S_u = zeros(127,128,365)
    S_v = zeros(128,127,365)

    hr_du = zeros(1023,1024,365)
    hr_dv = zeros(1024, 1023,365)
    hr_deta = zeros(1024,1024,365)

    hrcg_du = zeros(127,128,365)
    hrcg_dv = zeros(128,127,365)
    hrcg_deta = zeros(128,128,365)

    Shr_ = deepcopy(Shr)
    Slr_ = deepcopy(Slr)

    thr = 1800 * Shr.grid.dtint
    tlr = 225  * Slr.grid.dtint

    duhr = zeros(Shr.grid.nux, Shr.grid.nuy)
    dvhr = zeros(Shr.grid.nvx, Shr.grid.nvy)
    detahr = zeros(Shr.grid.nx, Shr.grid.ny)

    dulr = zeros(Slr.grid.nux, Slr.grid.nuy)
    dvlr = zeros(Slr.grid.nvx, Slr.grid.nvy)
    detalr = zeros(Slr.grid.nx, Slr.grid.ny)

    Mu = zeros(Slr.grid.nux, Slr.grid.nuy)
    Mv = zeros(Slr.grid.nvx, Slr.grid.nvy)
    ker = ImageFiltering.Kernel.gaussian((30e3/3750))
    for n = 1:365

        # if n ≤ 766
        #     uhr_, vhr_, etahr_ = ShallowWaters.add_halo(uhr1[:,:,n], vhr1[:,:,n], etahr1[:,:,n], Shr)
        # else
        #     uhr_, vhr_, etahr_ = ShallowWaters.add_halo(uhr2[:,:,n-766+1], vhr2[:,:,n-766+1], etahr2[:,:,n-766+1], Shr)
        # end

        uhr_, vhr_, etahr_ = ShallowWaters.add_halo(uhrall[:,:,n+1096], vhrall[:,:,n+1096], etahrall[:,:,n+1096], Shr)
        uhrcg_, vhrcg_, etahrcg_ = ShallowWaters.add_halo(uhrcgall[:,:,n+1096], vhrcgall[:,:,n], etahrcgall[:,:,n+1096], Slr)

        Shr_.Prog.u = uhr_
        Shr_.Prog.v = vhr_
        Shr_.Prog.η = etahr_

        Slr_.Prog.u = uhrcg_
        Slr_.Prog.v = vhrcg_
        Slr_.Prog.η = etahrcg_

        single_step!(duhr, dvhr, detahr, Shr_, (n + 1096 + 6)*thr)
        single_step!(dulr, dvlr, detalr, Slr_, (n + 1096 + 6)*tlr)

        @views hr_du[:,:,n] .= duhr
        @views hr_dv[:,:,n] .= dvhr
        @views hr_deta[:,:,n] .= detahr

        @views hrcg_du[:,:,n] .= dulr
        @views hrcg_dv[:,:,n] .= dvlr
        @views hrcg_deta[:,:,n] .= detalr

        # imfilter!(dufiltered, hr_du[:,:,n], reflect(ker))
        # imfilter!(dvfiltered, hr_dv[:,:,n], reflect(ker))

        # @views dudownsized = (dufiltered[8:8:end, 4:8:end, :] .+ dufiltered[8:8:end, 5:8:end, :]) ./ 2;
        # @views dvdownsized = (dvfiltered[4:8:end, 8:8:end, :] .+ dvfiltered[5:8:end, 8:8:end, :]) ./ 2;

        # @views S_u[:,:,n] .= hrcg_du[:,:,n] - dudownsized
        # @views S_v[:,:,n] .= hrcg_dv[:,:,n] - dvdownsized

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

function single_step!(du, dv, deta, S, t)

    # uold = copy(S.Prog.u)
    # vold = copy(S.Prog.v)
    # etaold = copy(S.Prog.η)

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

    u0rhs = S.Diag.PrognosticVarsRHS.u .= S.Diag.RungeKutta.u0
    v0rhs = S.Diag.PrognosticVarsRHS.v .= S.Diag.RungeKutta.v0
    # ShallowWaters.tracer!(i, u0rhs, v0rhs, chkp.S.Prog, chkp.S.Diag, chkp.S)

    # copyto!(S.Prog.u, S.Diag.RungeKutta.u0)
    # copyto!(S.Prog.v, S.Diag.RungeKutta.v0)
    # copyto!(S.Prog.η, S.Diag.RungeKutta.η0)

    halo = S.grid.halo
    haloη = S.grid.haloη

    du .= S.constants.scale_inv .* S.Diag.RungeKutta.u0[halo+1:end-halo, halo+1:end-halo] -
                S.constants.scale_inv .* S.Prog.u[halo+1:end-halo, halo+1:end-halo]

    dv .= S.constants.scale_inv .* S.Diag.RungeKutta.v0[halo+1:end-halo, halo+1:end-halo] -
                S.constants.scale_inv .* S.Prog.v[halo+1:end-halo, halo+1:end-halo]

    deta .= S.Diag.RungeKutta.η0[haloη+1:end-haloη, haloη+1:end-haloη] -
                S.Prog.η[haloη+1:end-haloη, haloη+1:end-haloη]

    return nothing

end

function single_step_diff!(Mu, Mv, S, t)

    # uold = copy(S.Prog.u)
    # vold = copy(S.Prog.v)
    # etaold = copy(S.Prog.η)

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

    u0rhs = S.Diag.PrognosticVarsRHS.u .= S.Diag.RungeKutta.u0
    v0rhs = S.Diag.PrognosticVarsRHS.v .= S.Diag.RungeKutta.v0
    # ShallowWaters.tracer!(i, u0rhs, v0rhs, chkp.S.Prog, chkp.S.Diag, chkp.S)

    # copyto!(S.Prog.u, S.Diag.RungeKutta.u0)
    # copyto!(S.Prog.v, S.Diag.RungeKutta.v0)
    # copyto!(S.Prog.η, S.Diag.RungeKutta.η0)

    Mu .= S.Diag.Smagorinsky.LLu1[:,2:end-1] + S.Diag.Smagorinsky.LLu2[2:end-1,:]
    Mv .= S.Diag.Smagorinsky.LLv1[:,2:end-1] + S.Diag.Smagorinsky.LLv2[2:end-1,:]

    return nothing

end

function save_viscosityterm()

    # u = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/u.nc", "u");
    # v = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/v.nc", "v");
    # eta = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/eta.nc", "eta");

    # u = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/u.nc", "u");
    # v = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/v.nc", "v");
    # eta = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/eta.nc", "eta");

    u = ncread("./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves/u.nc", "u");
    v = ncread("./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves/v.nc", "v");
    eta = ncread("./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves/eta.nc", "eta");

    # u = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_3years_dailysaves/u.nc", "u");
    # v = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_3years_dailysaves/v.nc", "v");
    # eta = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_3years_dailysaves/eta.nc", "eta");

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
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128,
        Ndays=3*365
    );
    Slr = ShallowWaters.model_setup(Plr);

    for S in [Slr]
        # calculate layer thicknesses for initial conditions
        ShallowWaters.thickness!(S.Diag.VolumeFluxes.h, S.Prog.η, S.forcing.H)
        ShallowWaters.Ix!(S.Diag.VolumeFluxes.h_u, S.Diag.VolumeFluxes.h)
        ShallowWaters.Iy!(S.Diag.VolumeFluxes.h_v, S.Diag.VolumeFluxes.h)
        ShallowWaters.Ixy!(S.Diag.Vorticity.h_q, S.Diag.VolumeFluxes.h)
    end

    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86daystart_3dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-6-8-10-13-15-18-23-28-33-38-41-44-48-51-53-58-63-65-68-73-78-83-86-88daystart_2dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution;
    # current = 1
    # for m in (S.Diag.CNNVars.model_Su, S.Diag.CNNVars.model_Sv)
    #     for layers in m[1]
    #         for array in layers
    #                 sz = prod(size(array))
    #                 array .= reshape(onlineweights[current:(current + sz - 1)], size(array)...)
    #                 current += sz
    #         end
    #     end
    # end

    Slr_ = deepcopy(Slr)

    tlr = 225 * Slr.grid.dtint
    Mu = zeros(Slr.grid.nux, Slr.grid.nuy)
    Mv = zeros(Slr.grid.nvx, Slr.grid.nvy)

    Mu_all = zeros(Slr.grid.nux, Slr.grid.nuy,1096)
    Mv_all = zeros(Slr.grid.nvx, Slr.grid.nvy,1096)

    for n = 1:1096

        u_, v_, eta_ = ShallowWaters.add_halo(u[:,:,n], v[:,:,n], eta[:,:,n], Slr)

        Slr_.Prog.u = u_
        Slr_.Prog.v = v_
        Slr_.Prog.η = eta_

        single_step_diff!(Mu, Mv, Slr_, n*tlr)

        @views Mu_all[:,:,n] .= Mu
        @views Mv_all[:,:,n] .= Mv

    end

end

