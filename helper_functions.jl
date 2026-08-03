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

function filter_hr(u, v)

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    filtered_u = imfilter(u, reflect(ker));
    filtered_v = imfilter(v, reflect(ker));

    return (filtered_u[8:8:end, 4:8:end] .+ filtered_u[8:8:end, 5:8:end]) .* 0.5, (filtered_v[4:8:end, 8:8:end] .+ filtered_v[5:8:end, 8:8:end]) .* 0.5;

end


function compute_transfer(u, v, Su, Sv)

        nfft = nextfastfft(size(u))

        outu, inputu, inputSu = paddingu(u, Su, nfft[1])
        fft2pow2radial!(outu, rfft(inputu), rfft(inputSu), nfft...)
        outv, inputv, inputSv = paddingv(v, Sv, nfft[1])
        fft2pow2radial!(outv, rfft(inputv), rfft(inputSv), nfft...)

    return outu, outv

end

"""
The rest of the functions here are related to computing momentum budget terms or tendencies

This is how I compute the total SGS forcing
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

    u = copy(S.Prog.u)
    v = copy(S.Prog.v)
    η = copy(S.Prog.η)

    ShallowWaters.thickness!(S.Diag.VolumeFluxes.h, η, S.forcing.H)
    ShallowWaters.Ix!(S.Diag.VolumeFluxes.h_u, S.Diag.VolumeFluxes.h)
    ShallowWaters.Iy!(S.Diag.VolumeFluxes.h_v, S.Diag.VolumeFluxes.h)
    ShallowWaters.Ixy!(S.Diag.Vorticity.h_q, S.Diag.VolumeFluxes.h)

    # calculate PV terms for initial conditions
    urhs = S.Diag.PrognosticVarsRHS.u .= u
    vrhs = S.Diag.PrognosticVarsRHS.v .= v
    ηrhs = S.Diag.PrognosticVarsRHS.η .= η

    ShallowWaters.advection_coriolis!(urhs, vrhs, ηrhs, S.Diag, S)
    ShallowWaters.PVadvection!(S.Diag, S)

    # propagate initial conditions
    copyto!(S.Diag.RungeKutta.u0, u)
    copyto!(S.Diag.RungeKutta.v0, v)
    copyto!(S.Diag.RungeKutta.η0, η)

    # store initial conditions of sst for relaxation
    copyto!(S.Diag.SemiLagrange.sst_ref, S.Prog.sst)

    # run a single step of integration loop

    # ghost point copy for boundary conditions
    ShallowWaters.ghost_points!(u, v, η, S)
    copyto!(S.Diag.RungeKutta.u1, u)
    copyto!(S.Diag.RungeKutta.v1, v)
    copyto!(S.Diag.RungeKutta.η1, η)

    if S.parameters.compensated
        fill!(S.Diag.Tendencies.du_sum, zero(S.parameters.Tprog))
        fill!(S.Diag.Tendencies.dv_sum, zero(S.parameters.Tprog))
        fill!(S.Diag.Tendencies.dη_sum, zero(S.parameters.Tprog))
    end

    if S.parameters.RKo > 2
        error("you need to set it to use rk2")
    end
    j = 1
    for rki = 1:1
        if rki > 1
            ShallowWaters.ghost_points!(
                S.Diag.RungeKutta.u1,
                S.Diag.RungeKutta.v1,
                S.Diag.RungeKutta.η1,
                S
            )
        end
        println("Number of times this print statement is hit: ", j)
        # type conversion for mixed precision
        u1rhs = S.Diag.PrognosticVarsRHS.u .= S.Diag.RungeKutta.u1
        v1rhs = S.Diag.PrognosticVarsRHS.v .= S.Diag.RungeKutta.v1
        η1rhs = S.Diag.PrognosticVarsRHS.η .= S.Diag.RungeKutta.η1

        println("Is the u being propagated still the initial u (before rhs is hit): ", norm(u1rhs .- S.Prog.u))
        println("Is the v being propagated still the initial v (before rhs is hit): ", norm(v1rhs .- S.Prog.v))
        println("Is the eta being propagated still the initial state (before rhs is hit): ", norm(η1rhs .- S.Prog.η))

        ShallowWaters.rhs!(u1rhs, v1rhs, η1rhs, S.Diag, S, t)          # momentum only
        ShallowWaters.continuity!(u1rhs, v1rhs, η1rhs, S.Diag, S, t)   # continuity equation

        if rki < S.parameters.RKo
            ShallowWaters.caxb!(
                S.Diag.RungeKutta.u1,
                u,
                S.constants.RKbΔt[rki],
                S.Diag.Tendencies.du
            )
            ShallowWaters.caxb!(
                S.Diag.RungeKutta.v1,
                v,
                S.constants.RKbΔt[rki],
                S.Diag.Tendencies.dv
            )
            ShallowWaters.caxb!(
                S.Diag.RungeKutta.η1,
                η,
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
    j += 1
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
            u,
            S.Diag.Tendencies.du_sum
        )
        ShallowWaters.dambmc!(
            S.Diag.Tendencies.dv_comp,
            S.Diag.RungeKutta.v0,
            v,
            S.Diag.Tendencies.dv_sum
        )
        ShallowWaters.dambmc!(
            S.Diag.Tendencies.dη_comp,
            S.Diag.RungeKutta.η0,
            η,
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
        # ShallowWaters.UVfluxes!(u0rhs, v0rhs, η0rhs, S.Diag, S)
        # ShallowWaters.advection_coriolis!(u0rhs, v0rhs, η0rhs, S.Diag, S)
    # end

    # if (chkp.i % S.grid.nstep_diff) == 0
        ShallowWaters.bottom_drag!(u0rhs, v0rhs, η0rhs, S.Diag, S)
        ShallowWaters.diffusion!(u0rhs, v0rhs, S.Diag, S)

    # original function call
        # ShallowWaters.add_drag_diff_tendencies!(
        #     S.Diag.RungeKutta.u0,
        #     S.Diag.RungeKutta.v0,
        #     S.Diag,
        #     S
        # )

    # modified to account for the half a timestep
        @unpack Bu,Bv = S.Diag.Bottomdrag
        @unpack LLu1,LLu2,LLv1,LLv2 = S.Diag.Smagorinsky
        @unpack halo,ep,Δt_diff = S.grid
        Tprog = S.parameters.Tprog

        m,n = size(S.Diag.RungeKutta.u0) .- (2*halo,2*halo)
        @boundscheck (m+2-ep,n+2) == size(Bu) || throw(BoundsError())
        @boundscheck (m,n+2) == size(LLu1) || throw(BoundsError())
        @boundscheck (m+2-ep,n) == size(LLu2) || throw(BoundsError())

        if S.parameters.zb_forcing_dissipation
            ShallowWaters.ZB_momentum(S.Diag.RungeKutta.u0,S.Diag.RungeKutta.v0,S,Diag)
        end

        if S.parameters.nn_forcing_dissipation
            ShallowWaters.CNN_momentum(S.Diag.RungeKutta.u0,S.Diag.RungeKutta.v0,S)
        end 

        # because the RK step only took one step with rk2 I technically did half a timestep, and need to be sure that the 
        # dissipative terms are also added with half a timestep
        # that's why here I add a multiplication by 2
        if S.parameters.zb_forcing_dissipation
        @inbounds for j ∈ 1:n
            for i ∈ 1:m
                S.Diag.RungeKutta.u0[i+2,j+2] += 0.5 * Δt_diff*(Tprog(Bu[i+1-ep,j+1]) + Tprog(LLu1[i,j+1]) + Tprog(LLu2[i+1-ep,j]) + Tprog(S.Diag.ZBVars.S_u[i,j]))
            end
        end
        elseif S.parameters.nn_forcing_dissipation
        @inbounds for j ∈ 1:n
            for i ∈ 1:m
                S.Diag.RungeKutta.u0[i+2,j+2] += 0.5 * Δt_diff*(Tprog(Bu[i+1-ep,j+1]) + Tprog(LLu1[i,j+1]) + Tprog(LLu2[i+1-ep,j]) + Tprog(S.Diag.CNNVars.S_u[i,j]))
            end
        end
        else
        @inbounds for j ∈ 1:n
            for i ∈ 1:m
                S.Diag.RungeKutta.u0[i+2,j+2] += 0.5 * Δt_diff*(Tprog(Bu[i+1-ep,j+1]) + Tprog(LLu1[i,j+1]) + Tprog(LLu2[i+1-ep,j]))
            end
        end
        end

        m,n = size(S.Diag.RungeKutta.v0) .- (2*halo,2*halo)
        @boundscheck (m+2,n+2) == size(Bv) || throw(BoundsError())
        @boundscheck (m,n+2) == size(LLv1) || throw(BoundsError())
        @boundscheck (m+2,n) == size(LLv2) || throw(BoundsError())

        if S.parameters.zb_forcing_dissipation
            @inbounds for j ∈ 1:n
                for i ∈ 1:m
                    S.Diag.RungeKutta.v0[i+2,j+2] += 0.5 * Δt_diff*(Tprog(Bv[i+1,j+1]) + Tprog(LLv1[i,j+1]) + Tprog(LLv2[i+1,j]) + Tprog(S.Diag.ZBVars.S_v[i,j]))
                end
            end
        elseif S.parameters.nn_forcing_dissipation
            @inbounds for j ∈ 1:n
                for i ∈ 1:m
                    S.Diag.RungeKutta.v0[i+2,j+2] += 0.5 * Δt_diff*(Tprog(Bv[i+1,j+1]) + Tprog(LLv1[i,j+1]) + Tprog(LLv2[i+1,j]) + Tprog(S.Diag.CNNVars.S_v[i,j]))
                end
            end
        else
            @inbounds for j ∈ 1:n
                for i ∈ 1:m
                    S.Diag.RungeKutta.v0[i+2,j+2] += 0.5 * Δt_diff*(Tprog(Bv[i+1,j+1]) + Tprog(LLv1[i,j+1]) + Tprog(LLv2[i+1,j]))
                end
            end
        end
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
                S.constants.scale_inv .* u[halo+1:end-halo, halo+1:end-halo]

    dv .= S.constants.scale_inv .* S.Diag.RungeKutta.v0[halo+1:end-halo, halo+1:end-halo] -
                S.constants.scale_inv .* v[halo+1:end-halo, halo+1:end-halo]

    deta .= S.Diag.RungeKutta.η0[haloη+1:end-haloη, haloη+1:end-haloη] -
                η[haloη+1:end-haloη, haloη+1:end-haloη]

    return nothing

end

"""
I wrote compute_advection!() to compute the nonlinear advection portion of the momentum equation. 
To use it you need to set S.parameters.g = 1e-12 and S.parameters.ω=1e-16. This will get rid of 
f_q (the Coriolis force) and g \eta in the Bernoulli potential calculaiton.

Old update:
There was a bug in the original function I'd written, I don't know what it was, but my Bernoulli potential that I was
computing was wrong. However, this has been corrected such that the difference of the 
full momentum (compute_advection!() with g and f\eta) and the euler update *without dissipation* (compute_tendencies_witheuler!() modified) is precisely 
the wind-stress forcing, which is what it should be. The new plan is to just force the coriolis and gravity to be zero and then compute the advection term
with the above compute_advection!() function

Me describing how I went about extracting the nonlinear advection
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
"""
function compute_advection!(mom_u, mom_v, S, t)

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
        # ShallowWaters.UVfluxes!(u0rhs, v0rhs, η0rhs, S.Diag, S)
        # ShallowWaters.advection_coriolis!(u0rhs, v0rhs, η0rhs, S.Diag, S)
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

"""
To get at the momentum budget terms, we need
    - nonlinear advection
    - Coriolis force
    - gravity term -g ∇η
    - wind forcing
    - bottom drag
    - viscosity
To compute these we can use the above functions.
    - The nonlinear advection is returned by compute_advection!(), just make sure to set S.parameters.g and S.parameters.ω to be sufficiently small
    - The coriolis force is added as part of the potential vorticity computation, and should be in
        S.grid.f_q
    - The gravity term is part of the Bernoulli potential calculation, and doesn't have a specific spot, but is probably accessed by
            m,n = size(p)
            one_half_scale_inv = convert(T,0.5)*scale_inv
            @inbounds for j ∈ 1:n
                for i ∈ 1:m
                    g*η[i,j] = p[i,j] - one_half_scale_inv*(KEu[i+ep,j+1] + KEv[i+1,j])
                end
            end
    and then taking derivatives of the result
    - Wind forcing is constant, just pull from model structure
    - Bottom drag is in S.Diag.Bottomdrag.Bu and S.Diag.Bottomdrag.Bv, pull from model structure
    - Viscosity is in S.Diag.Smagorinsky.LLu1[i,j+1] + S.Diag.Smagorinsky.LLu2[i+1-ep,j]
"""
function compute_mom_budget(Seuler, Sadv, t, n)

    du = zeros(Seuler.grid.nux, Seuler.grid.nuy);
    dv = zeros(Seuler.grid.nvx, Seuler.grid.nvy);
    deta = zeros(Seuler.grid.nx, Seuler.grid.ny);

    adv_u = zeros(Sadv.grid.nux, Sadv.grid.nuy);
    adv_v = zeros(Sadv.grid.nvx, Sadv.grid.nvy);

    compute_tendencies_witheuler!(du, dv, deta, Seuler, n*t)
    compute_advection!(adv_u, adv_v, Sadv, n*t)

    m,n = size(Seuler.Diag.Bernoulli.p)
    detagdx = zeros(m-1, n)
    detagdy = zeros(m, n-1)

    ShallowWaters.∂x!(detagdx, (Seuler.constants.g .* Seuler.Prog.η)./ Seuler.grid.Δ)
    ShallowWaters.∂y!(detagdy, (Seuler.constants.g .* Seuler.Prog.η)./ Seuler.grid.Δ)

    # getting the Coriolis force out of the potential vorticity computation
    m,n = size(Seuler.Diag.Vorticity.q)
    coriolisoverh = zeros(m,n)
    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            coriolisoverh[i,j] = Seuler.Diag.Vorticity.q[i,j] - (Seuler.Diag.Vorticity.dvdx[i+1,j+1] - Seuler.Diag.Vorticity.dudy[i+1,j+1]) / Seuler.Diag.Vorticity.h_q[i,j]
        end
    end

    @unpack qα,qβ,qγ,qδ = Seuler.Diag.ArakawaHsu
    ShallowWaters.AHα!(qα,coriolisoverh)
    ShallowWaters.AHβ!(qβ,coriolisoverh)
    ShallowWaters.AHγ!(qγ,coriolisoverh)
    ShallowWaters.AHδ!(qδ,coriolisoverh)

    m,n = size(Seuler.Diag.Vorticity.qhv)
    fv = zeros(m,n)
    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            fv[i,j] = qα[1+i,j]*Seuler.Diag.VolumeFluxes.V[2+i,j+1] +
                qβ[1+i,j]*Seuler.Diag.VolumeFluxes.V[1+i,j+1] +
                qγ[1+i,j]*Seuler.Diag.VolumeFluxes.V[1+i,j] +
                qδ[1+i,j]*Seuler.Diag.VolumeFluxes.V[2+i,j]
        end
    end

    m,n = size(Seuler.Diag.Vorticity.qhu)
    fu = zeros(m,n)
    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            fu[i,j] = qα[i,j]*Seuler.Diag.VolumeFluxes.U[i,j+1] +
                qβ[i+1,j]*Seuler.Diag.VolumeFluxes.U[i+1,j+1] +
                qγ[i+1,j+1]*Seuler.Diag.VolumeFluxes.U[i+1,j+2] +
                qδ[i,j+1]*Seuler.Diag.VolumeFluxes.U[i,j+2]
        end
    end

    # all of the rhs terms are computed with a *scaled* prognostic variable, so to actually 
    # see the momentum budget terms we need to remove that scaling. The division by Delta is because
    # all of the timestepping is done with dt / Delta, and we want the term that gets multiplied by dt
    # I should have corrected for the scalings in the values returned by compute_mom_budget already,
    # but for the stuff I'm pulling after (see below) I need to do it manually

    scaling = (Seuler.grid.scale * Seuler.grid.Δ)
    bottomdragu = (Seuler.Diag.Bottomdrag.Bu[2:end-1,2:end-1]) ./ scaling
    bottomdragv = (Seuler.Diag.Bottomdrag.Bv[2:end-1,2:end-1]) ./ scaling
    viscu = (Seuler.Diag.Smagorinsky.LLu1[:,2:end-1] + Seuler.Diag.Smagorinsky.LLu2[2:end-1,:]) ./ scaling;
    viscv = (Seuler.Diag.Smagorinsky.LLv1[:, 2:end-1] + Seuler.Diag.Smagorinsky.LLv2[2:end-1,:]) ./ scaling

    Fx = Seuler.forcing.Fx ./ (Seuler.grid.scale * Seuler.grid.Δ)

    return du, dv, adv_u, adv_v, fu./ (Seuler.grid.scale .* Seuler.grid.Δ), fv./ (Seuler.grid.scale .* Seuler.grid.Δ), detagdx[2:end-1, 2:end-1], detagdy[2:end-1, 2:end-1], bottomdragu, bottomdragv, viscu, viscv, Fx

end

"""
Just for running the above functions, largely not used anymore
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

    # high-resolution model
    t = 1800 * S.grid.dtint

    du = zeros(S.grid.nux, S.grid.nuy);
    dv = zeros(S.grid.nvx, S.grid.nvy);
    deta = zeros(S.grid.nx, S.grid.ny);

    filtereddu = zeros(127, 128)
    filtereddv = zeros(128, 127)

    du_all = zeros(127, 128, 1096);
    dv_all = zeros(128, 127, 1096);

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
        RKo=2,
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

    adv_u = zeros(S.grid.nux, S.grid.nuy);
    adv_v = zeros(S.grid.nvx, S.grid.nvy);

    adv_u_all = zeros(127, 128, 1096);
    adv_v_all = zeros(128, 127, 1096);

    ker = ImageFiltering.Kernel.gaussian((30e3/3750))

    for n = 1:1096

        u_, v_, eta_ = ShallowWaters.add_halo(u[:,:,n], v[:,:,n], eta[:,:,n], S)

        S.Prog.u = u_
        S.Prog.v = v_
        S.Prog.η = eta_

        compute_advection!(adv_u, adv_v, S, n*t)

        filteredadvu = imfilter(adv_u, reflect(ker))
        filteredadvv = imfilter(adv_v, reflect(ker))

        @views adv_u_all[:,:,n] .= (filteredadvu[8:8:end, 4:8:end] .+ filteredadvu[8:8:end, 5:8:end]) .* 0.5
        @views adv_v_all[:,:,n] .= (filteredadvv[4:8:end, 8:8:end] .+ filteredadvv[5:8:end, 8:8:end]) .* 0.5

    end

    return adv_u_all, adv_v_all

end