
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
