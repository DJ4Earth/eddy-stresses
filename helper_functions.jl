"""
computes the approximate nonlinear advection sgs forcing
"""
function compute_approxhrS()

    @views u = uhrall;
    @views v = vhrall;
    @views eta = etahrall;

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
        cfl=.898,
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

    S_u = zeros(T,Slr.Diag.CNNVars.nux,Slr.Diag.CNNVars.nuy, 1096)             # total forcing in x-direction
    S_v = zeros(T,Slr.Diag.CNNVars.nvx,Slr.Diag.CNNVars.nvy, 1096)             # total forcing in y-direction

    uhrT = zeros(1024,1024)
    vhrT = zeros(1024,1024)

    s = Slr.grid.Δ^2
    alpha = 0.1
    winu = tukey((Shr.grid.nux, Shr.grid.nuy), alpha)
    winv = tukey((Shr.grid.nvx,Shr.grid.nvy), alpha)
    totalstates = 1096
    for t = 1:totalstates

        uhrh = cat(zeros(T,1023+2*halo,halo),cat(zeros(T,halo,1024),winu.*u[:,:,t],zeros(T,halo,1024),dims=1),zeros(T,1023+2*halo,halo),dims=2)
        vhrh = cat(zeros(T,1024+2*halo,halo),cat(zeros(T,halo,1023),winv.*v[:,:,t],zeros(T,halo,1023),dims=1),zeros(T,1024+2*halo,halo),dims=2)

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

"""
these are functions related to computing the sum over shells of a Fourier transform, needed when computing the KE transfer
These are pulled from the DSP.jl Julia package
"""
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

"""
The rest of the functions here are related to computing sgs forcing terms
"""

"""
There are effectively two ways to compute S_tot, both are technically correct and its a question of which is the 
"best" way to go about it. Given the equation
    u_t = f(u),
we have 
    overline{u_t} = overline{f(u)}
the coarse-grained high-resolution equation and we also have,
    overline{u}_t = f(overline{u}) + S,
the equation for a coarse-grained u. Then, as we want them equal this gives
    S = S_{tot} = overline{f(u)} - f(overline{u})
the first way of defining the total SGS forcing. But this doesn't take into account the discretization process.
ShallowWaters uses RK4 by default for the time-stepping, and it's how I computed my high-resolution states. If we 
discretize the first equation, we have
    u_t = u_{t - Delta t} + Delta t f(u)
If we only did a first order time-stepping scheme, then this would be the same as the first way of computing 
S_{tot}. However, that's not the case, we have higher order terms that appear in RK4
    u_t = u_{t - Delta t} + four terms
which is *not* equal to the first order method, and thus not equal to the first way of computing S_{tot}. This leads to the 
second way of getting at S_{tot}, which takes into account these higher order terms in the RK4 step. This second way is what I 
used to compute S_{tot}, and what gave the figure currently in the manuscript. It's a question of which is better, I'd argue it's the 
second, because again we used RK4 for our high-resolution variables, and it makes sense that the NN learned something about this
"""
function compute_tendencies_withrk!(du, dv, deta, S, t)

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

function compute_tendencies_witheuler!(du, dv, deta, S, t)

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

    if S.parameters.RKo > 2
        error("you need to set it to use rk2")
    end

    for rki = 1:1
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

    # println("Norm of q:", norm(S.Diag.Vorticity.q))
    # println("Norm of U:", norm(S.Diag.VolumeFluxes.U))
    # println("Norm of V:", norm(S.Diag.VolumeFluxes.V))
    # println("Norm of qalpha:", norm(S.Diag.ArakawaHsu.qα))
    # println("Norm of qbeta:", norm(S.Diag.ArakawaHsu.qβ))
    # println("Norm of qgamma:", norm(S.Diag.ArakawaHsu.qγ))
    # println("Norm of qdelta:", norm(S.Diag.ArakawaHsu.qδ))
    # println("Norm of qhv:", norm(S.Diag.Vorticity.qhv))

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

function compute_momentum_new!(mom_u, mom_v, S, t)

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

    if S.parameters.RKo > 2
        error("you need to set it to use rk2")
    end

    for rki = 1:1
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
        # ShallowWaters.bottom_drag!(u0rhs, v0rhs, η0rhs, S.Diag, S)
        # ShallowWaters.diffusion!(u0rhs, v0rhs, S.Diag, S)
        # ShallowWaters.add_drag_diff_tendencies!(
        #     S.Diag.RungeKutta.u0,
        #     S.Diag.RungeKutta.v0,
        #     S.Diag,
        #     S
        # )
        # ShallowWaters.ghost_points_uv!(
        #     S.Diag.RungeKutta.u0,
        #     S.Diag.RungeKutta.v0,
        #     S
        # )
    # end

    t += S.grid.dtint

    u0rhs = S.Diag.PrognosticVarsRHS.u .= S.Diag.RungeKutta.u0
    v0rhs = S.Diag.PrognosticVarsRHS.v .= S.Diag.RungeKutta.v0
    # ShallowWaters.tracer!(i, u0rhs, v0rhs, chkp.S.Prog, chkp.S.Diag, chkp.S)

    # copyto!(S.Prog.u, S.Diag.RungeKutta.u0)
    # copyto!(S.Prog.v, S.Diag.RungeKutta.v0)
    # copyto!(S.Prog.η, S.Diag.RungeKutta.η0)

    # @unpack q,dvdx,dudy,h_q = S.Diag.Vorticity
    # @unpack f_q,ep = S.grid

    # m,n = size(S.Diag.Vorticity.q)
    # @inbounds for j ∈ 1:n
    #     for i ∈ 1:m
    #         S.Diag.Vorticity.q[i,j] = (S.grid.f_q[i,j] + S.Diag.Vorticity.dvdx[i+1,j+1] - S.Diag.Vorticity.dudy[i+1+S.grid.ep,j+1]) / S.Diag.Vorticity.h_q[i,j]
    #     end
    # end

    # # @unpack qα,qβ,qγ,qδ = S.Diag.ArakawaHsu.
    # ShallowWaters.AHα!(S.Diag.ArakawaHsu.qα,S.Diag.Vorticity.q)
    # ShallowWaters.AHβ!(S.Diag.ArakawaHsu.qβ,S.Diag.Vorticity.q)
    # ShallowWaters.AHγ!(S.Diag.ArakawaHsu.qγ,S.Diag.Vorticity.q)
    # ShallowWaters.AHδ!(S.Diag.ArakawaHsu.qδ,S.Diag.Vorticity.q)

    # @unpack U,V = S.Diag.VolumeFluxes.
    # @unpack qhv,qhu,ep = S.Diag.Vorticity
    # @unpack qα,qβ,qγ,qδ = S.Diag.ArakawaHsu

    # Linear combinations of q and V=hv to yield qhv
    # ep = S.grid.ep
    # m,n = size(S.Diag.Vorticity.qhv)
    # @inbounds for j ∈ 1:n
    #     for i ∈ 1:m
    #         S.Diag.Vorticity.qhv[i,j] = S.Diag.ArakawaHsu.qα[1-ep+i,j]*S.Diag.VolumeFluxes.V[2-ep+i,j+1] +
    #         S.Diag.ArakawaHsu.qβ[1-ep+i,j]*S.Diag.VolumeFluxes.V[1-ep+i,j+1] +
    #         S.Diag.ArakawaHsu.qγ[1-ep+i,j]*S.Diag.VolumeFluxes.V[1-ep+i,j] +
    #         S.Diag.ArakawaHsu.qδ[1-ep+i,j]*S.Diag.VolumeFluxes.V[2-ep+i,j]
    #     end
    # end

    # Linear combinations of q and U=hu to yield qhu
    # m,n = size(qhu)
    # @inbounds for j ∈ 1:n
    #     for i ∈ 1:m
    #         qhu[i,j] = qα[i,j]*U[i,j+1] + qβ[i+1,j]*U[i+1,j+1] + qγ[i+1,j+1]*U[i+1,j+2] + qδ[i,j+1]*U[i,j+2]
    #     end
    # end

    # m,n = size(p)
    # @boundscheck (m+ep,n+2) == size(KEu) || throw(BoundsError())
    # @boundscheck (m+2,n) == size(KEv) || throw(BoundsError())
    # @boundscheck (m,n) == size(η) || throw(BoundsError())

    # one_half_scale_inv = convert(T,0.5)*scale_inv

    # @inbounds for j ∈ 1:n
    #     for i ∈ 1:m
    #         p[i,j] = one_half_scale_inv*(KEu[i+ep,j+1] + KEv[i+1,j]) + g*η[i,j]
    #     end
    # end

    # ShallowWaters.∂x!(dpdx,p)
    # ShallowWaters.∂y!(dpdy,p)

    ep = S.grid.ep
    m,n = size(S.Diag.Tendencies.du) .- (2*S.grid.halo,2*S.grid.halo)
    for j = 1:n
        for i = 1:m
        mom_u[i,j] = S.Diag.Vorticity.qhv[i,j] - S.Diag.Bernoulli.dpdx[i+1-ep,j+1]
        end
    end

    m,n = size(S.Diag.Tendencies.dv) .- (2*S.grid.halo,2*S.grid.halo)
    for j = 1:n
        for i = 1:m
        mom_v[i,j] = -S.Diag.Vorticity.qhu[i,j] - S.Diag.Bernoulli.dpdy[i+1,j+1]
        end
    end

    mom_u .= (S.constants.scale_inv .* mom_u) ./ S.grid.Δ
    mom_v .= (S.constants.scale_inv .* mom_v) ./ S.grid.Δ

    return nothing

end

# function compute_momentum!(mom_u, mom_v, S, t)

#     # calculate layer thicknesses for initial conditions
#     ShallowWaters.thickness!(S.Diag.VolumeFluxes.h, S.Prog.η, S.forcing.H)
#     ShallowWaters.Ix!(S.Diag.VolumeFluxes.h_u, S.Diag.VolumeFluxes.h)
#     ShallowWaters.Iy!(S.Diag.VolumeFluxes.h_v, S.Diag.VolumeFluxes.h)
#     ShallowWaters.Ixy!(S.Diag.Vorticity.h_q, S.Diag.VolumeFluxes.h)

#     # propagate initial conditions
#     copyto!(S.Diag.RungeKutta.u0, S.Prog.u)
#     copyto!(S.Diag.RungeKutta.v0, S.Prog.v)
#     copyto!(S.Diag.RungeKutta.η0, S.Prog.η)

#     ShallowWaters.ghost_points!(S.Prog.u, S.Prog.v, S.Prog.η, S)
#     copyto!(S.Diag.RungeKutta.u1, S.Prog.u)
#     copyto!(S.Diag.RungeKutta.v1, S.Prog.v)
#     copyto!(S.Diag.RungeKutta.η1, S.Prog.η)

#     # type conversion for mixed precision
#     u1rhs = S.Diag.PrognosticVarsRHS.u .= S.Diag.RungeKutta.u1
#     v1rhs = S.Diag.PrognosticVarsRHS.v .= S.Diag.RungeKutta.v1
#     η1rhs = S.Diag.PrognosticVarsRHS.η .= S.Diag.RungeKutta.η1

#     u = S.Prog.u
#     v = S.Prog.v
#     halo = S.grid.halo

#     # this chunk will compute ShallowWaters.UVFluxes(u,v,\eta, Diag,S) #########
#     @unpack h,h_u,h_v,U,V = S.Diag.VolumeFluxes
#     @unpack H = S.forcing
#     @unpack ep = S.grid
#     @unpack scale_inv = S.constants

#     ShallowWaters.thickness!(h,S.Prog.η,H)
#     ShallowWaters.Ix!(h_u,h)
#     ShallowWaters.Iy!(h_v,h)

#     # mass or volume flux U,V = uh,vh
#     ShallowWaters.Uflux!(U,u,h_u,ep,scale_inv)
#     ShallowWaters.Vflux!(V,v,h_v,scale_inv)

#     # Next we want to compute ShallowWaters.advection_coriolis!(S.Prog.u, S.Prog.v, S.Prog.\eta, S.Diag, S) #########
#     @unpack h = S.Diag.VolumeFluxes
#     @unpack q,h_q,dvdx,dudy = S.Diag.Vorticity
#     @unpack u²,v²,KEu,KEv = S.Diag.Bernoulli
#     @unpack ep,f_q = S.grid

#     ShallowWaters.Ixy!(h_q,h)

#     # off-diagonals of stress tensor ∇(u,v)
#     ShallowWaters.∂x!(dvdx,v)
#     ShallowWaters.∂y!(dudy,u)

#     # non-linear part of the Bernoulli potential
#     ShallowWaters.speed!(u²,v²,u,v)
#     ShallowWaters.Ix!(KEu,u²)
#     ShallowWaters.Iy!(KEv,v²)

#     # the potential vorticity computation with the Coriolis force
#     m,n = size(q)
#     @inbounds for j ∈ 1:n
#         for i ∈ 1:m
#             q[i,j] = (f_q[i,j] + dvdx[i+1,j+1] - dudy[i+1+ep,j+1]) / h_q[i,j]
#         end
#     end

#     # Now we want to compute ShallowWaters.PVadvection!(S.Diag, S) ######
#     @unpack qα,qβ,qγ,qδ = S.Diag.ArakawaHsu
#     ShallowWaters.AHα!(qα,q)
#     ShallowWaters.AHβ!(qβ,q)
#     ShallowWaters.AHγ!(qγ,q)
#     ShallowWaters.AHδ!(qδ,q)

#     m,n = size(S.Diag.Vorticity.qhv)
#     qhv = S.Diag.Vorticity.qhv
#     @inbounds for j ∈ 1:n
#         for i ∈ 1:m
#             qhv[i,j] = qα[1-ep+i,j]*V[2-ep+i,j+1] + qβ[1-ep+i,j]*V[1-ep+i,j+1] + qγ[1-ep+i,j]*V[1-ep+i,j] + qδ[1-ep+i,j]*V[2-ep+i,j]
#         end
#     end

#     m,n = size(S.Diag.Vorticity.qhu)
#     qhu = S.Diag.Vorticity.qhu
#     @inbounds for j ∈ 1:n
#         for i ∈ 1:m
#             qhu[i,j] = qα[i,j]*U[i,j+1] + qβ[i+1,j]*U[i+1,j+1] + qγ[i+1,j+1]*U[i+1,j+2] + qδ[i,j+1]*U[i,j+2]
#         end
#     end

#     # lastly the piece from Bernoulli needed for the advection term
#     p = S.Diag.Bernoulli.p
#     scale_inv = S.constants.scale_inv
#     g = S.constants.g

#     m,n = size(p)
#     @boundscheck (m+ep,n+2) == size(KEu) || throw(BoundsError())
#     @boundscheck (m+2,n) == size(KEv) || throw(BoundsError())
#     @boundscheck (m,n) == size(S.Prog.η) || throw(BoundsError())

#     one_half_scale_inv = convert(Float64,0.5)*scale_inv

#     # computing the Bernoulli potential
#     @inbounds for j ∈ 1:n
#         for i ∈ 1:m
#             p[i,j] = one_half_scale_inv*(KEu[i+ep,j+1] + KEv[i+1,j]) + g*S.Prog.η[i,j]
#         end
#     end

#     ShallowWaters.∂x!(S.Diag.Bernoulli.dpdx, p)
#     ShallowWaters.∂y!(S.Diag.Bernoulli.dpdy, p)

#     m,n = size(S.Diag.Tendencies.du) .- (2*halo,2*halo)
#     for j = 1:n
#         for i = 1:m
#         mom_u[i,j] = S.Diag.Vorticity.qhv[i,j] - S.Diag.Bernoulli.dpdx[i+1-ep,j+1]
#         end
#     end

#     m,n = size(S.Diag.Tendencies.dv) .- (2*halo,2*halo)
#     for j = 1:n
#         for i = 1:m
#         mom_v[i,j] = -S.Diag.Vorticity.qhu[i,j] - S.Diag.Bernoulli.dpdy[i+1,j+1]
#         end
#     end

#     mom_u .= (S.constants.scale_inv .* mom_u) ./ S.grid.Δ
#     mom_v .= (S.constants.scale_inv .* mom_v) ./ S.grid.Δ

#     return nothing

# end

"""
Now we want to compute the rest of the momentum terms (including linear terms so that I can
compute the momentum budget)
"""

"""
Computes the viscosity term given a u, v in S
"""
function compute_viscosity_new!(Mu, Mv, S)

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

    @unpack dudx,dvdy,dvdx,dudy = S.Diag.Vorticity

    u = S.Prog.u
    v = S.Prog.v
    Diag = S.Diag

    ShallowWaters.∂x!(dudx,u)
    ShallowWaters.∂y!(dvdy,v)
    ShallowWaters.∂x!(dvdx,v)
    ShallowWaters.∂y!(dudy,u)

    # biharmonic diffusion
    ShallowWaters.stress_tensor!(u,v,Diag)
    ShallowWaters.smagorinsky_coeff!(Diag,S)
    ShallowWaters.viscous_tensor_smagorinsky!(Diag)

    @unpack LLu1,LLu2,LLv1,LLv2 = Diag.Smagorinsky
    @unpack S11,S12,S21,S22 = Diag.Smagorinsky

    ShallowWaters.∂x!(LLu1,S11)
    ShallowWaters.∂y!(LLu2,S12)
    ShallowWaters.∂x!(LLv1,S21)
    ShallowWaters.∂y!(LLv2,S22)

    Mu .= S.Diag.Smagorinsky.LLu1[:,2:end-1] + S.Diag.Smagorinsky.LLu2[2:end-1,:]
    Mv .= S.Diag.Smagorinsky.LLv1[:,2:end-1] + S.Diag.Smagorinsky.LLv2[2:end-1,:]

    return nothing

end

function compute_viscosity_old!(Mu, Mv, S, t)

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

"""
Computes the bottom drag, again need to provide a u, v in S
"""
function compute_bottomdrag!(B_u, B_v, S)

    Diag = S.Diag

    @unpack h,h_u,h_v = Diag.VolumeFluxes
    @unpack u²,v²,KEu,KEv = Diag.Bernoulli
    @unpack Bu,Bv,sqrtKE,sqrtKE_u,sqrtKE_v = Diag.Bottomdrag
    @unpack ep, halo = S.grid
    @unpack cD = S.constants
    @unpack H = S.forcing

    u = S.Prog.u
    v = S.Prog.v
    η = S.Prog.η

    ShallowWaters.thickness!(h,η,H)
    ShallowWaters.Ix!(h_u,h)
    ShallowWaters.Iy!(h_v,h)

    ShallowWaters.speed!(u²,v²,u,v)
    ShallowWaters.Ix!(KEu,u²)
    ShallowWaters.Iy!(KEv,v²)

    # sqrt of KE, which is actually the kinetic energy without the 0.5 factor
    m,n = size(sqrtKE)
    @boundscheck (m+ep,n+2) == size(KEu) || throw(BoundsError())
    @boundscheck (m+2,n) == size(KEv) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            sqrtKE[i,j] = sqrt(KEu[i+ep,j+1] + KEv[i+1,j])
        end
    end

    ShallowWaters.Ix!(sqrtKE_u,sqrtKE)
    ShallowWaters.Iy!(sqrtKE_v,sqrtKE)

    m,n = size(Bu)
    @boundscheck (m,n) == size(sqrtKE_u) || throw(BoundsError())
    @boundscheck (m,n) == size(h_u) || throw(BoundsError())
    @boundscheck (m+2+ep,n+2) == size(u) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            Bu[i,j] = cD*sqrtKE_u[i,j] * u[i+1+ep,j+1] / h_u[i,j]
        end
    end

    m,n = size(Bv)
    @boundscheck (m,n) == size(sqrtKE_v) || throw(BoundsError())
    @boundscheck (m,n) == size(h_v) || throw(BoundsError())
    @boundscheck (m+2,n+2) == size(v) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            Bv[i,j] = cD*sqrtKE_v[i,j] * v[i+1,j+1] / h_v[i,j]
        end
    end

    m,n = size(u) .- (2*halo,2*halo)
    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            B_u[i,j] = Bu[i+1-ep,j+1]
        end
    end

    m,n = size(v) .- (2*halo,2*halo)
    @inbounds for j ∈ 1:n
        for i ∈ 1:m 
            B_v[i,j] = Bv[i+1,j+1]
        end
    end

    return nothing

end

"""
Coriolis force computation and others -- needs to be modified but this should contain all other terms in the 
momentum equation
"""
function remaining_terms()
  # calculate layer thicknesses for initial conditions
    ShallowWaters.thickness!(S.Diag.VolumeFluxes.h, S.Prog.η, S.forcing.H)
    ShallowWaters.Ix!(S.Diag.VolumeFluxes.h_u, S.Diag.VolumeFluxes.h)
    ShallowWaters.Iy!(S.Diag.VolumeFluxes.h_v, S.Diag.VolumeFluxes.h)
    ShallowWaters.Ixy!(S.Diag.Vorticity.h_q, S.Diag.VolumeFluxes.h)

    u = S.Prog.u
    v = S.Prog.v
    halo = S.grid.halo

    # this chunk will compute ShallowWaters.UVFluxes(u,v,\eta, Diag,S) #########
    @unpack h,h_u,h_v,U,V = S.Diag.VolumeFluxes
    @unpack H = S.forcing
    @unpack ep = S.grid
    @unpack scale_inv = S.constants

    ShallowWaters.thickness!(h,S.Prog.η,H)
    ShallowWaters.Ix!(h_u,h)
    ShallowWaters.Iy!(h_v,h)

    # mass or volume flux U,V = uh,vh
    ShallowWaters.Uflux!(U,u,h_u,ep,scale_inv)
    ShallowWaters.Vflux!(V,v,h_v,scale_inv)

    # Next we want to compute ShallowWaters.advection_coriolis!(S.Prog.u, S.Prog.v, S.Prog.\eta, S.Diag, S) #########
    @unpack h = S.Diag.VolumeFluxes
    @unpack q,h_q,dvdx,dudy = S.Diag.Vorticity
    @unpack u²,v²,KEu,KEv = S.Diag.Bernoulli
    @unpack ep,f_q = S.grid

    ShallowWaters.Ixy!(h_q,h)

    # off-diagonals of stress tensor ∇(u,v)
    ShallowWaters.∂x!(dvdx,v)
    ShallowWaters.∂y!(dudy,u)

    # non-linear part of the Bernoulli potential
    ShallowWaters.speed!(u²,v²,u,v)
    ShallowWaters.Ix!(KEu,u²)
    ShallowWaters.Iy!(KEv,v²)

    # the potential vorticity computation *without* the Coriolis force
    m,n = size(q)
    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            q[i,j] = (f_q[i,j] + dvdx[i+1,j+1] - dudy[i+1+ep,j+1]) / h_q[i,j]
        end
    end

    # Now we want to compute ShallowWaters.PVadvection!(S.Diag, S) ######
    @unpack qα,qβ,qγ,qδ = S.Diag.ArakawaHsu
    ShallowWaters.AHα!(qα,q)
    ShallowWaters.AHβ!(qβ,q)
    ShallowWaters.AHγ!(qγ,q)
    ShallowWaters.AHδ!(qδ,q)

    m,n = size(S.Diag.Vorticity.qhv)
    qhv = S.Diag.Vorticity.qhv
    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            qhv[i,j] = qα[1-ep+i,j]*V[2-ep+i,j+1] + qβ[1-ep+i,j]*V[1-ep+i,j+1] + qγ[1-ep+i,j]*V[1-ep+i,j] + qδ[1-ep+i,j]*V[2-ep+i,j]
        end
    end

    m,n = size(S.Diag.Vorticity.qhu)
    qhu = S.Diag.Vorticity.qhu
    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            qhu[i,j] = qα[i,j]*U[i,j+1] + qβ[i+1,j]*U[i+1,j+1] + qγ[i+1,j+1]*U[i+1,j+2] + qδ[i,j+1]*U[i,j+2]
        end
    end

    # lastly the piece from Bernoulli needed for the advection term
    p = S.Diag.Bernoulli.p
    scale_inv = S.constants.scale_inv

    m,n = size(p)
    @boundscheck (m+ep,n+2) == size(KEu) || throw(BoundsError())
    @boundscheck (m+2,n) == size(KEv) || throw(BoundsError())
    @boundscheck (m,n) == size(S.Prog.η) || throw(BoundsError())

    one_half_scale_inv = convert(Float64,0.5)*scale_inv

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            p[i,j] = one_half_scale_inv*(KEu[i+ep,j+1] + KEv[i+1,j]) + S.constants.g*S.Prog.η[i,j]
        end
    end

    # removing the g*\eta and then computing the derivatives
    ShallowWaters.∂x!(S.Diag.Bernoulli.dpdx, p - (S.constants.g .* S.Prog.η))
    ShallowWaters.∂y!(S.Diag.Bernoulli.dpdy, p - (S.constants.g .* S.Prog.η))

end

"""
Wind-stress is found by just pulling from the model structure, it doesn't change with time
"""

"""
The following function is to compute the nonlinear advection
Milan's code hides the nonlinear advection away inside of p = 1/2 (u^2 + v^2) + gh and 
q = (f + zeta) / h. What gets added to the tendencies is ultimately
    dudt portion: qhv - partial_x p
    dvdt portion: -qhu - partial_y p
Then we have
    qhv = ((f + zeta) / h) * h * v = (f + zeta) * v = fv + zeta v = fv + v_x v - u_y v
    qhu = ((f + zeta / h)) * h * u = (f + zeta) * u = fu + zeta u = fu + (u v_x - u u_y)
    partial_x p = .5 * (2 u_x + 2 v_x) + partial_x (g h) = u u_x + v v_x + partial_x (g h)
    partial_y p = u_y + v_y + partial_y (g h)
If I want to isolate u u_x + v u_y then I need to 
    (1) find qhv = fv + v_x v - u_y v
    (2) find partial_x p = u u_x + v v_x + partial_x (g h)
    (3) Their difference is fv - v u_y - u u_x - partial_x (g h)
    (4) Then I just need to get rid of the Coriolis force and that partial term
For a given u, v the order of functions in ShallowWaters computation for qhv and qhu is
    UVFluxes
    advection_coriolis
    PV_advection
    Bernoulli
and lastly the actual rhs computation. The below function is pulling each of these pieces separately in given order

Update on 07/21/26

There was a bug in the old function, I don't know what it was, but my Bernoulli potential that I used to be 
computing was wrong. However, this has been corrected such that the difference of the 
full momentum (compute_momentum_new!()) and the euler update *without dissipation* (compute_tendencies_witheuler!()) is precisely 
the wind-stress forcing, which is what it should be. The new plan is to just force the coriolis and gravity to be zero and then compute the advection term
with the above compute_momentum!() function
"""
# function compute_advection!(adv_u, adv_v, S, t)

#     ShallowWaters.thickness!(S.Diag.VolumeFluxes.h, S.Prog.η, S.forcing.H)
#     ShallowWaters.Ix!(S.Diag.VolumeFluxes.h_u, S.Diag.VolumeFluxes.h)
#     ShallowWaters.Iy!(S.Diag.VolumeFluxes.h_v, S.Diag.VolumeFluxes.h)
#     ShallowWaters.Ixy!(S.Diag.Vorticity.h_q, S.Diag.VolumeFluxes.h)

#     # calculate PV terms for initial conditions
#     urhs = S.Diag.PrognosticVarsRHS.u .= S.Prog.u
#     vrhs = S.Diag.PrognosticVarsRHS.v .= S.Prog.v
#     ηrhs = S.Diag.PrognosticVarsRHS.η .= S.Prog.η

#     ShallowWaters.advection_coriolis!(urhs, vrhs, ηrhs, S.Diag, S)
#     ShallowWaters.PVadvection!(S.Diag, S)

#     # propagate initial conditions
#     copyto!(S.Diag.RungeKutta.u0, S.Prog.u)
#     copyto!(S.Diag.RungeKutta.v0, S.Prog.v)
#     copyto!(S.Diag.RungeKutta.η0, S.Prog.η)

#     # store initial conditions of sst for relaxation
#     copyto!(S.Diag.SemiLagrange.sst_ref, S.Prog.sst)

#     # run a single step of integration loop

#     # ghost point copy for boundary conditions
#     ShallowWaters.ghost_points!(S.Prog.u, S.Prog.v, S.Prog.η, S)
#     copyto!(S.Diag.RungeKutta.u1, S.Prog.u)
#     copyto!(S.Diag.RungeKutta.v1, S.Prog.v)
#     copyto!(S.Diag.RungeKutta.η1, S.Prog.η)

#     if S.parameters.compensated
#         fill!(S.Diag.Tendencies.du_sum, zero(S.parameters.Tprog))
#         fill!(S.Diag.Tendencies.dv_sum, zero(S.parameters.Tprog))
#         fill!(S.Diag.Tendencies.dη_sum, zero(S.parameters.Tprog))
#     end

#     if S.parameters.RKo > 2
#         error("you need to set it to use rk2")
#     end

#     for rki = 1:1
#         if rki > 1
#             ShallowWaters.ghost_points!(
#                 S.Diag.RungeKutta.u1,
#                 S.Diag.RungeKutta.v1,
#                 S.Diag.RungeKutta.η1,
#                 S
#             )
#         end

#         # type conversion for mixed precision
#         u1rhs = S.Diag.PrognosticVarsRHS.u .= S.Diag.RungeKutta.u1
#         v1rhs = S.Diag.PrognosticVarsRHS.v .= S.Diag.RungeKutta.v1
#         η1rhs = S.Diag.PrognosticVarsRHS.η .= S.Diag.RungeKutta.η1

#         ShallowWaters.rhs!(u1rhs, v1rhs, η1rhs, S.Diag, S, t)          # momentum only
#         ShallowWaters.continuity!(u1rhs, v1rhs, η1rhs, S.Diag, S, t)   # continuity equation

#         if rki < S.parameters.RKo
#             ShallowWaters.caxb!(
#                 S.Diag.RungeKutta.u1,
#                 S.Prog.u,
#                 S.constants.RKbΔt[rki],
#                 S.Diag.Tendencies.du
#             )
#             ShallowWaters.caxb!(
#                 S.Diag.RungeKutta.v1,
#                 S.Prog.v,
#                 S.constants.RKbΔt[rki],
#                 S.Diag.Tendencies.dv
#             )
#             ShallowWaters.caxb!(
#                 S.Diag.RungeKutta.η1,
#                 S.Prog.η,
#                 S.constants.RKbΔt[rki],
#                 S.Diag.Tendencies.dη
#             )
#         end

#         if S.parameters.compensated
#             ShallowWaters.axb!(S.Diag.Tendencies.du_sum, S.constants.RKaΔt[rki], S.Diag.Tendencies.du)
#             ShallowWaters.axb!(S.Diag.Tendencies.dv_sum, S.constants.RKaΔt[rki], S.Diag.Tendencies.dv)
#             ShallowWaters.axb!(S.Diag.Tendencies.dη_sum, S.constants.RKaΔt[rki], S.Diag.Tendencies.dη)
#         else
#             ShallowWaters.axb!(
#                 S.Diag.RungeKutta.u0,
#                 S.constants.RKaΔt[rki],
#                 S.Diag.Tendencies.du
#             )
#             ShallowWaters.axb!(
#                 S.Diag.RungeKutta.v0,
#                 S.constants.RKaΔt[rki],
#                 S.Diag.Tendencies.dv
#             )
#             ShallowWaters.axb!(
#                 S.Diag.RungeKutta.η0,
#                 S.constants.RKaΔt[rki],
#                 S.Diag.Tendencies.dη
#             )
#         end
#     end

#     if S.parameters.compensated
#         ShallowWaters.axb!(S.Diag.Tendencies.du_sum, -1, S.Diag.Tendencies.du_comp)
#         ShallowWaters.axb!(S.Diag.Tendencies.dv_sum, -1, S.Diag.Tendencies.dv_comp)
#         ShallowWaters.axb!(S.Diag.Tendencies.dη_sum, -1, S.Diag.Tendencies.dη_comp)

#         ShallowWaters.axb!(S.Diag.RungeKutta.u0, 1, S.Diag.Tendencies.du_sum)
#         ShallowWaters.axb!(S.Diag.RungeKutta.v0, 1, S.Diag.Tendencies.dv_sum)
#         ShallowWaters.axb!(S.Diag.RungeKutta.η0, 1, S.Diag.Tendencies.dη_sum)

#         ShallowWaters.dambmc!(
#             S.Diag.Tendencies.du_comp,
#             S.Diag.RungeKutta.u0,
#             S.Prog.u,
#             S.Diag.Tendencies.du_sum
#         )
#         ShallowWaters.dambmc!(
#             S.Diag.Tendencies.dv_comp,
#             S.Diag.RungeKutta.v0,
#             S.Prog.v,
#             S.Diag.Tendencies.dv_sum
#         )
#         ShallowWaters.dambmc!(
#             S.Diag.Tendencies.dη_comp,
#             S.Diag.RungeKutta.η0,
#             S.Prog.η,
#             S.Diag.Tendencies.dη_sum
#         )
#     end

#     ShallowWaters.ghost_points!(
#         S.Diag.RungeKutta.u0,
#         S.Diag.RungeKutta.v0,
#         S.Diag.RungeKutta.η0,
#         S
#     )

#     u0rhs = S.Diag.PrognosticVarsRHS.u .= S.Diag.RungeKutta.u0
#     v0rhs = S.Diag.PrognosticVarsRHS.v .= S.Diag.RungeKutta.v0
#     η0rhs = S.Diag.PrognosticVarsRHS.η .= S.Diag.RungeKutta.η0

#     # if S.parameters.dynamics == "nonlinear" && S.grid.nstep_advcor > 0 && (i % S.grid.nstep_advcor) == 0
#         ShallowWaters.UVfluxes!(u0rhs, v0rhs, η0rhs, S.Diag, S)
#         ShallowWaters.advection_coriolis!(u0rhs, v0rhs, η0rhs, S.Diag, S)
#     # end

#     # if (chkp.i % S.grid.nstep_diff) == 0
#         # ShallowWaters.bottom_drag!(u0rhs, v0rhs, η0rhs, S.Diag, S)
#         # ShallowWaters.diffusion!(u0rhs, v0rhs, S.Diag, S)
#         # ShallowWaters.add_drag_diff_tendencies!(
#         #     S.Diag.RungeKutta.u0,
#         #     S.Diag.RungeKutta.v0,
#         #     S.Diag,
#         #     S
#         # )
#         # ShallowWaters.ghost_points_uv!(
#         #     S.Diag.RungeKutta.u0,
#         #     S.Diag.RungeKutta.v0,
#         #     S
#         # )
#     # end

#     t += S.grid.dtint

#     u0rhs = S.Diag.PrognosticVarsRHS.u .= S.Diag.RungeKutta.u0
#     v0rhs = S.Diag.PrognosticVarsRHS.v .= S.Diag.RungeKutta.v0
#     # ShallowWaters.tracer!(i, u0rhs, v0rhs, chkp.S.Prog, chkp.S.Diag, chkp.S)

#     # copyto!(S.Prog.u, S.Diag.RungeKutta.u0)
#     # copyto!(S.Prog.v, S.Diag.RungeKutta.v0)
#     # copyto!(S.Prog.η, S.Diag.RungeKutta.η0)

#     ep = S.grid.ep
#     m,n = size(S.Diag.Tendencies.du) .- (2*S.grid.halo,2*S.grid.halo)
#     for j = 1:n
#         for i = 1:m
#         adv_u[i,j] = S.Diag.Vorticity.qhv[i,j] - S.Diag.Bernoulli.dpdx[i+1-ep,j+1]
#         end
#     end

#     m,n = size(S.Diag.Tendencies.dv) .- (2*S.grid.halo,2*S.grid.halo)
#     for j = 1:n
#         for i = 1:m
#         adv_v[i,j] = -S.Diag.Vorticity.qhu[i,j] - S.Diag.Bernoulli.dpdy[i+1,j+1]
#         end
#     end

#     adv_u .= (S.constants.scale_inv .* adv_u) ./ S.grid.Δ
#     adv_v .= (S.constants.scale_inv .* adv_v) ./ S.grid.Δ

#     return nothing

# end


"""
Just for running the above functions
"""

function save_rktendencies()

    @views u = uhrcgall;
    @views v = vhrcgall;
    @views eta = etahrcgall;

    @views u = uhrall;
    @views v = vhrall;
    @views eta = etahrall;

    # u = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/u.nc", "u");
    # v = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/v.nc", "v");
    # eta = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/eta.nc", "eta");

    # u = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/u.nc", "u");
    # v = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/v.nc", "v");
    # eta = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/eta.nc", "eta");

    # u = ncread("./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves/u.nc", "u");
    # v = ncread("./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves/v.nc", "v");
    # eta = ncread("./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves/eta.nc", "eta");

    # u = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_3years_dailysaves/u.nc", "u");
    # v = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_3years_dailysaves/v.nc", "v");
    # eta = ncread("./dissipation_constant/results/result_online_multistateweights_10dayoptimization_5-20-35-50-65-75initdays_startfrom20daystate_3years_dailysaves/eta.nc", "eta");

    Phr = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        # cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
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
        Ndays=3*365
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
        RKo=2,
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

    S = Shr
    # for S in [Slr]
        # calculate layer thicknesses for initial conditions
        ShallowWaters.thickness!(S.Diag.VolumeFluxes.h, S.Prog.η, S.forcing.H)
        ShallowWaters.Ix!(S.Diag.VolumeFluxes.h_u, S.Diag.VolumeFluxes.h)
        ShallowWaters.Iy!(S.Diag.VolumeFluxes.h_v, S.Diag.VolumeFluxes.h)
        ShallowWaters.Ixy!(S.Diag.Vorticity.h_q, S.Diag.VolumeFluxes.h)
    # end

    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86daystart_3dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-6-8-10-13-15-18-23-28-33-38-41-44-48-51-53-58-63-65-68-73-78-83-86-88daystart_2dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution;
    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_5-20-35-50-65-75daystart_10dayoptimization_initialweights20daystate_fixedcfl_15iterations_constdissipation.jld2").solution;
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

    # high-resolution model
    t = 1800 * S.grid.dtint

    # low-resolution model
    t = 225 * S.grid.dtint

    du = zeros(S.grid.nux, S.grid.nuy);
    dv = zeros(S.grid.nvx, S.grid.nvy);
    deta = zeros(S.grid.nx, S.grid.ny);
    du_all = zeros(S.grid.nux, S.grid.nuy, 1096);
    dv_all = zeros(S.grid.nvx, S.grid.nvy, 1096);
    deta_all = zeros(S.grid.nx, S.grid.ny, 1096);

    alpha = 0.1
    winu = tukey((Shr.grid.nux, Shr.grid.nuy), alpha)
    winv = tukey((Shr.grid.nvx, Shr.grid.nvy), alpha)
    # wineta = tukey((Shr.grid.nx, Shr.grid.ny), alpha)

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    for n = 1:1096

        # ufiltered = imfilter(winu.*u[:,:,n], reflect(ker))
        # vfiltered = imfilter(winv.*v[:,:,n], reflect(ker))
        # etafiltered = imfilter(eta[:,:,n], reflect(ker))

        # @views windowuhrdownsized = (ufiltered[8:8:end, 4:8:end] .+ ufiltered[8:8:end, 5:8:end]) .* 0.5
        # @views windowvhrdownsized = (vfiltered[4:8:end, 8:8:end] .+ vfiltered[5:8:end, 8:8:end]) .* 0.5
        # @views windowetahrdownsized = (etafiltered[4:8:end,4:8:end] .+ etafiltered[5:8:end,5:8:end] .+ etafiltered[4:8:end,5:8:end] .+ etafiltered[5:8:end,4:8:end]) ./ 4

        # u_, v_, eta_ = ShallowWaters.add_halo(windowuhrdownsized, windowvhrdownsized, windowetahrdownsized, S)

        u_, v_, eta_ = ShallowWaters.add_halo(u[:,:,n], v[:,:,n], eta[:,:,n], S)

        S.Prog.u = u_
        S.Prog.v = v_
        S.Prog.η = eta_

        # single_step_diff!(Bu, Bv, Mu, Mv, S, n*t)
        # compute_advection!(adv_u, adv_v, S, n*t)
        # compute_tendencies_withrk!(1, du, dv, deta, S, n*t)
        compute_tendencies_witheuler!(du, dv, deta, S, n*t)
        # saving tendencies
        @views du_all[:,:,n] .= du
        @views dv_all[:,:,n] .= dv
        @views deta_all[:,:,n] .= deta

    end

end

function save_bottomdraglr()

    @views u = uhrcgall;
    @views v = vhrcgall;
    @views eta = etahrcgall;

    # u = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/u.nc", "u");
    # v = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/v.nc", "v");
    # eta = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/eta.nc", "eta");

    # u = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/u.nc", "u");
    # v = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/v.nc", "v");
    # eta = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/eta.nc", "eta");

    # u = ncread("./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves/u.nc", "u");
    # v = ncread("./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves/v.nc", "v");
    # eta = ncread("./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves/eta.nc", "eta");

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
        RKo=2,
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

    S = Slr

    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86daystart_3dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-6-8-10-13-15-18-23-28-33-38-41-44-48-51-53-58-63-65-68-73-78-83-86-88daystart_2dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution;
    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_5-20-35-50-65-75daystart_10dayoptimization_initialweights20daystate_fixedcfl_15iterations_constdissipation.jld2").solution;
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

    # low-resolution model
    t = 225 * S.grid.dtint

    Bu = zeros(S.grid.nux, S.grid.nuy)
    Bv = zeros(S.grid.nvx, S.grid.nvy)

    Bu_all = zeros(S.grid.nux, S.grid.nuy, 1096)
    Bv_all = zeros(S.grid.nvx, S.grid.nvy, 1096)

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    for n = 1:1096

        u_, v_, eta_ = ShallowWaters.add_halo(u[:,:,n], v[:,:,n], eta[:,:,n], S)

        S.Prog.u = u_
        S.Prog.v = v_
        S.Prog.η = eta_

        compute_bottomdrag!(Bu, Bv, S)

        @views Bu_all[:,:,n] .= Bu
        @views Bv_all[:,:,n] .= Bv

    end

    return Bu_all, Bv_all

end

function save_bottomdraghr()

    @views u = uhrall;
    @views v = vhrall;
    @views eta = etahrall;

    Phr = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        # cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
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
        Ndays=3*365
    );
    Shr = ShallowWaters.model_setup(Phr);

    S = Shr

    # low-resolution model
    t = 1800 * S.grid.dtint

    Bu = zeros(S.grid.nux, S.grid.nuy)
    Bv = zeros(S.grid.nvx, S.grid.nvy)

    filteredMu = zeros(127, 128)
    filteredMv = zeros(128, 127)

    Bu_all = zeros(127, 128, 1096)
    Bv_all = zeros(128, 127, 1096)

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    for n = 1:1096

        u_, v_, eta_ = ShallowWaters.add_halo(u[:,:,n], v[:,:,n], eta[:,:,n], S)

        S.Prog.u = u_
        S.Prog.v = v_
        S.Prog.η = eta_

        compute_bottomdrag!(Bu, Bv, S)

        filteredBu = imfilter(Bu, reflect(ker))
        filteredBv = imfilter(Bv, reflect(ker))

        @views Bu_all[:,:,n] .= (filteredBu[8:8:end, 4:8:end] .+ filteredBu[8:8:end, 5:8:end]) .* 0.5
        @views Bv_all[:,:,n] .= (filteredBv[4:8:end, 8:8:end] .+ filteredBv[5:8:end, 8:8:end]) .* 0.5

    end

    return Bu_all, Bv_all

end

function save_viscositylr()

    @views u = uhrcgall;
    @views v = vhrcgall;
    @views eta = etahrcgall;

    # u = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/u.nc", "u");
    # v = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/v.nc", "v");
    # eta = ncread("./dissipation_constant/results/result_online_multistateweights_3dayoptimization_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86initdays_startfrommulti3_3years_dailysaves/eta.nc", "eta");

    # u = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/u.nc", "u");
    # v = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/v.nc", "v");
    # eta = ncread("./dissipation_constant/results/result_online_multistateweights_2dayoptimization_startfrommulti3_3years_dailysaves/eta.nc", "eta");

    # u = ncread("./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves/u.nc", "u");
    # v = ncread("./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves/v.nc", "v");
    # eta = ncread("./dissipation_constant/spinup_files/128_ZBparam_postspinup_cginitcond_3years_dailysaves/eta.nc", "eta");

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
        RKo=4,
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

    S = Slr

    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86daystart_3dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution
    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_1-4-6-8-10-13-15-18-23-28-33-38-41-44-48-51-53-58-63-65-68-73-78-83-86-88daystart_2dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution;
    # onlineweights = load_object("./dissipation_constant/tuned_weights/result_multistate_5-20-35-50-65-75daystart_10dayoptimization_initialweights20daystate_fixedcfl_15iterations_constdissipation.jld2").solution;
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

    # low-resolution model
    t = 225 * S.grid.dtint

    Mu = zeros(S.grid.nux, S.grid.nuy)
    Mv = zeros(S.grid.nvx, S.grid.nvy)

    Mu_all = zeros(S.grid.nux, S.grid.nuy, 1096)
    Mv_all = zeros(S.grid.nvx, S.grid.nvy, 1096)

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    for n = 1:1096

        u_, v_, eta_ = ShallowWaters.add_halo(u[:,:,n], v[:,:,n], eta[:,:,n], S)

        S.Prog.u = u_
        S.Prog.v = v_
        S.Prog.η = eta_

        compute_viscosity_old!(Mu, Mv, S, t*n)

        @views Mu_all[:,:,n] .= Mu
        @views Mv_all[:,:,n] .= Mv

    end

    return Mu_all, Mv_all

end

function save_viscosityhr()

    @views u = uhrall;
    @views v = vhrall;
    @views eta = etahrall;

    Phr = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        # cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
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
        Ndays=3*365
    );
    Shr = ShallowWaters.model_setup(Phr);

    S = Shr

    t = 1800 * S.grid.dtint

    Mu = zeros(S.grid.nux, S.grid.nuy)
    Mv = zeros(S.grid.nvx, S.grid.nvy)

    filteredMu = zeros(127, 128)
    filteredMv = zeros(128, 127)

    Mu_all = zeros(127, 128, 1096)
    Mv_all = zeros(128, 127, 1096)

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    for n = 1:1096

        u_, v_, eta_ = ShallowWaters.add_halo(u[:,:,n], v[:,:,n], eta[:,:,n], S)

        S.Prog.u = u_
        S.Prog.v = v_
        S.Prog.η = eta_

        # compute_viscosity!(Mu, Mv, S)
        compute_viscosity_old!(Mu, Mv, S, t*n)

        filteredMu = imfilter(Mu, reflect(ker))
        filteredMv = imfilter(Mv, reflect(ker))

        @views Mu_all[:,:,n] .= (filteredMu[8:8:end, 4:8:end] .+ filteredMu[8:8:end, 5:8:end]) .* 0.5
        @views Mv_all[:,:,n] .= (filteredMv[4:8:end, 8:8:end] .+ filteredMv[5:8:end, 8:8:end]) .* 0.5

    end

    return Mu_all, Mv_all

end

function save_eulerlr()

    @views u = uhrcgall;
    @views v = vhrcgall;
    @views eta = etahrcgall;

    Plr = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
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

    S = Slr

    # low-resolution model
    t = 225 * S.grid.dtint

    du = zeros(S.grid.nux, S.grid.nuy);
    dv = zeros(S.grid.nvx, S.grid.nvy);
    deta = zeros(S.grid.nx, S.grid.ny);
    du_all = zeros(S.grid.nux, S.grid.nuy, 1096);
    dv_all = zeros(S.grid.nvx, S.grid.nvy, 1096);
    deta_all = zeros(S.grid.nx, S.grid.ny, 1096);

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    for n = 1:1096

        u_, v_, eta_ = ShallowWaters.add_halo(u[:,:,n], v[:,:,n], eta[:,:,n], S)

        S.Prog.u = u_
        S.Prog.v = v_
        S.Prog.η = eta_

        compute_tendencies_witheuler!(du, dv, deta, S, n*t)
        # saving tendencies
        @views du_all[:,:,n] .= du
        @views dv_all[:,:,n] .= dv
        @views deta_all[:,:,n] .= deta

    end

    return du_all, dv_all, deta_all

end

function save_eulerhr()

    @views u = uhrall;
    @views v = vhrall;
    @views eta = etahrall;

    Phr = ShallowWaters.Parameter(T=Float64,
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        # cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        RKo=2,
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
        Ndays=3*365
    );
    Shr = ShallowWaters.model_setup(Phr);

    S = Shr

    ShallowWaters.thickness!(S.Diag.VolumeFluxes.h, S.Prog.η, S.forcing.H)
    ShallowWaters.Ix!(S.Diag.VolumeFluxes.h_u, S.Diag.VolumeFluxes.h)
    ShallowWaters.Iy!(S.Diag.VolumeFluxes.h_v, S.Diag.VolumeFluxes.h)
    ShallowWaters.Ixy!(S.Diag.Vorticity.h_q, S.Diag.VolumeFluxes.h)

    # high-resolution model
    t = 1800 * S.grid.dtint

    du = zeros(S.grid.nux, S.grid.nuy);
    dv = zeros(S.grid.nvx, S.grid.nvy);
    deta = zeros(S.grid.nx, S.grid.ny);

    filtereddu = zeros(127, 128)
    filtereddv = zeros(128, 127)

    du_all = zeros(127, 128, 1096);
    dv_all = zeros(128, 127, 1096);
    deta_all = zeros(S.grid.nx, S.grid.ny, 1096);

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    for n = 1:1096

        u_, v_, eta_ = ShallowWaters.add_halo(u[:,:,n], v[:,:,n], eta[:,:,n], S)

        S.Prog.u = u_
        S.Prog.v = v_
        S.Prog.η = eta_

        compute_tendencies_witheuler!(du, dv, deta, S, n*t)

        filtereddu = imfilter(du, reflect(ker))
        filtereddv = imfilter(dv, reflect(ker))

        @views du_all[:,:,n] .= (filtereddu[8:8:end, 4:8:end] .+ filtereddu[8:8:end, 5:8:end]) .* 0.5
        @views dv_all[:,:,n] .= (filtereddv[4:8:end, 8:8:end] .+ filtereddv[5:8:end, 8:8:end]) .* 0.5

        # @views deta_all[:,:,n] .= deta

    end

    return du_all, dv_all, deta_all

end

function save_adveclr()

    @views u = uhrcgall;
    @views v = vhrcgall;
    @views eta = etahrcgall;

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

    S = Slr
    t = 225 * S.grid.dtint

    adv_u = zeros(S.grid.nux, S.grid.nuy);
    adv_v = zeros(S.grid.nvx, S.grid.nvy);

    adv_u_all = zeros(S.grid.nux, S.grid.nuy, 1096);
    adv_v_all = zeros(S.grid.nvx, S.grid.nvy, 1096);

    for n = 1:1096

        u_, v_, eta_ = ShallowWaters.add_halo(u[:,:,n], v[:,:,n], eta[:,:,n], S)

        S.Prog.u = u_
        S.Prog.v = v_
        S.Prog.η = eta_

        compute_advection!(adv_u, adv_v, S, n*t)

        # saving the advection
        @views adv_u_all[:,:,n] .= adv_u
        @views adv_v_all[:,:,n] .= adv_v


    end

    return adv_u_all, adv_v_all

end

function save_advechr()

    @views u = uhrall;
    @views v = vhrall;
    @views eta = etahrall;

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
        Ndays=3*365
    );
    Shr = ShallowWaters.model_setup(Phr);

    S = Shr
    t = 1800 * S.grid.dtint

    adv_u = zeros(S.grid.nux, S.grid.nuy);
    adv_v = zeros(S.grid.nvx, S.grid.nvy);

    adv_u_all = zeros(S.grid.nux, S.grid.nuy, 1096);
    adv_v_all = zeros(S.grid.nvx, S.grid.nvy, 1096);

    for n = 1:1096

        u_, v_, eta_ = ShallowWaters.add_halo(u[:,:,n], v[:,:,n], eta[:,:,n], S)

        S.Prog.u = u_
        S.Prog.v = v_
        S.Prog.η = eta_

        compute_advection!(adv_u, adv_v, S, n*t)

        # saving the advection
        @views adv_u_all[:,:,n] .= adv_u
        @views adv_v_all[:,:,n] .= adv_v


    end

    return adv_u_all, adv_v_all

end